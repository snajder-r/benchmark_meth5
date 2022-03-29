import tqdm
from collections import namedtuple
from matplotlib.patches import Patch
import random
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from meth5.meth5 import MetH5File
import numpy as np
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density

from benchmark_pycometh.config import module_config
import math

matplotlib.use("Agg")

pa = PlotArchiver("segmentation", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

seg_methylkit = pd.read_csv(
    module_config.methylkit_segmentation["parents_mock_from_np_hp"],
    sep="\t",
    usecols=[0, 1, 2],
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
    skiprows=1,
)
seg_methylkit["chrom"] = seg_methylkit["chrom"].map(lambda x: x.replace("chr", ""))
seg_nanoepiseg = pd.read_csv(
    module_config.nanoepiseg_segmentation["parents"], sep="\t", names=["chrom", "start", "end"], dtype={"chrom": "str"}
)

seg_methcp = pd.read_csv(
    module_config.methcp_segmentation["parents_mock_from_np_hp"],
    sep="\t",
    usecols=[1, 2, 3],
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
    skiprows=1,
)


mfs = {sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r") for sample in ["HG003", "HG004"]}
hp_dict = {
    sample: {int(k): v for k, v in mfs[sample].h5_fp["reads/read_groups/haplotype"].attrs.items()} for sample in mfs
}
hp_ids = {sample: [k for k, v in hp_dict[sample].items() if v in {"H1", "H2"}] for sample in mfs}

seg_methylkit = seg_methylkit.groupby("chrom")
seg_nanoepiseg = seg_nanoepiseg.groupby("chrom")
seg_methcp = seg_methcp.groupby("chrom")

chrom = "21"


def compute_variance(mf, regions, wiggle=0):
    regions = list(regions)
    vars = []
    lens = []
    with tqdm.tqdm(total=len(regions)) as pbar:
        for region in regions:
            if region.start < wiggle:
                continue
            agg, _ = mf[chrom].get_values_in_range(region.start - wiggle, region.end - wiggle).get_llr_site_rate()
            var = np.nanvar(agg)
            if not np.isnan(var):
                lens.append(region.end - region.start)
                vars.append(var)
            pbar.update(1)
    vars = np.array(vars)
    lens = np.array(lens)
    return lens, vars


def permute_segments(original_segments):
    segment_lengths = original_segments.apply(lambda x: x["end"] - x["start"], axis=1).tolist()
    ends = np.array([0] + original_segments["end"].tolist())[:-1]
    starts = np.array(original_segments["start"])
    gaps = starts - ends
    offset = gaps[0]
    gaps = gaps[1:]
    random.shuffle(segment_lengths)
    random.shuffle(gaps)
    gaps = iter(gaps)
    region = namedtuple("region", ["start", "end"])
    for l in segment_lengths:
        yield region(offset, offset + l)
        try:
            offset += next(gaps) + l
        except StopIteration:
            return


segmentations = {"PycoMeth": seg_nanoepiseg, "MethylKit": seg_methylkit, "MethCP": seg_methcp}
callers = list(segmentations.keys())
variances = {sample:{} for sample in mfs}
seglengths = {sample:{} for sample in mfs}

for caller_name, segs in segmentations.items():
    chrsegs = segs.get_group(chrom)
    for sample in callers:
        seglengths[sample][caller_name], variances[sample][caller_name] = compute_variance(mfs[sample], chrsegs.itertuples())
        caller_name = f"{caller_name} randomized"
        randomsegs = list(permute_segments(chrsegs))
        seglengths[sample][caller_name], variances[sample][caller_name] = compute_variance(mfs[sample], randomsegs)


with pa.open_multipage_pdf("variance_segmentation_parents"):
    for sample in variances:
        pa.figure()
        plt.title(f"Variance of segmentation {sample} chr21")
        x = callers
        y = [variances[sample][c] for c in x]
        parts = plt.violinplot(y, positions=[0, 1, 2])
        for pc in parts["bodies"]:
            pc.set_facecolor("blue")
            pc.set_alpha(0.5)
        
        x_r = [f"{c} randomized" for c in callers]
        y = [variances[sample][c] for c in x_r]
        parts = plt.violinplot(y, positions=[0, 1, 2])
        for pc in parts["bodies"]:
            pc.set_facecolor("red")
            pc.set_alpha(0.5)
        
        legend_elements = [
            Patch(facecolor="red", edgecolor="red", label="Random"),
            Patch(facecolor="blue", edgecolor="blue", label="Segmentation"),
        ]
        plt.legend(handles=legend_elements)
        plt.xticks([0, 1, 2], x)
        pa.savefig()

        x_max = 0
        y_max = 0
        for caller_name in variances[sample]:
            y_max = max(y_max, max(variances[sample][caller_name]))
            x_max = max(x_max, max(seglengths[sample][caller_name]))
        
        for caller_name in variances[sample]:
            pa.figure()
            plt.title(f"Length vs variance {sample} chr21 {caller_name}")
            x = np.log10(seglengths[sample][caller_name])
            idx = ~(np.isnan(x) | np.isinf(x))
            plot_2d_density(x[idx], variances[sample][caller_name][idx], cmap="jet")
            plt.xlabel("Segment length (log10)")
            plt.ylabel("Variance")
            plt.ylim(0, y_max)
            plt.xlim(0, np.log10(x_max))
            pa.savefig()
