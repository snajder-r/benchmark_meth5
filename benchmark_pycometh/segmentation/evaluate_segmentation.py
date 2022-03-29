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
    module_config.methylkit_segmentation["HG003_mock_from_np"],
    sep="\t",
    usecols=[0, 1, 2],
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
    skiprows=1,
)
seg_methylkit["chrom"] = seg_methylkit["chrom"].map(lambda x: x.replace("chr", ""))
seg_nanoepiseg = pd.read_csv(
    module_config.nanoepiseg_segmentation["HG003"], sep="\t", names=["chrom", "start", "end"], dtype={"chrom": "str"}
)

seg_methcp = pd.read_csv(
    module_config.methcp_segmentation["parents_mock_from_np_hp"],
    sep="\t",
    usecols=[1, 2, 3],
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
    skiprows=1,
)



mf = MetH5File(module_config.meth5_template_file.format(sample="HG003"), "r")
hp_dict = {int(k): v for k, v in mf.h5_fp["reads/read_groups/haplotype"].attrs.items()}
hp_ids = [k for k, v in hp_dict.items() if v in {"H1", "H2"}]

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
                lens.append(region.end-region.start)
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


np_variance = compute_variance(mf, seg_nanoepiseg.get_group(chrom).itertuples())
mk_variance = compute_variance(mf, seg_methylkit.get_group(chrom).itertuples())
mcp_variance = compute_variance(mf, seg_methcp.get_group(chrom).itertuples())
np_variance_random = compute_variance(mf, list(permute_segments(seg_nanoepiseg.get_group(chrom))))
mk_variance_random = compute_variance(mf, list(permute_segments(seg_methylkit.get_group(chrom))))
mcp_variance_random = compute_variance(mf, list(permute_segments(seg_methcp.get_group(chrom))))


with pa.open_multipage_pdf("variance_segmentation"):
    pa.figure()
    plt.title("Variance of segmentation HG003 chr21")
    parts = plt.violinplot((np_variance[1], mk_variance[1], mcp_variance[1]), positions=[0, 1, 2])
    for pc in parts["bodies"]:
        pc.set_facecolor("blue")
        pc.set_alpha(0.5)
        
    parts = plt.violinplot((np_variance_random[1], mk_variance_random[1], mcp_variance_random[1]), positions=[0, 1, 2])
    for pc in parts["bodies"]:
        pc.set_facecolor("red")
        pc.set_alpha(0.5)
        

    legend_elements = [Patch(facecolor='red', edgecolor='red', label='Random'), Patch(facecolor='blue', edgecolor='blue', label='Segmentation')]
    plt.legend(handles=legend_elements)
    plt.xticks([0, 1, 2], ["PycoMeth", "MethylKit", "MethCP"])
    pa.savefig()

    for caller_name, variance in [("PycoMeth", np_variance), ("MethylKit", mk_variance), ("MethCP", mcp_variance), ("PycoMeth randomized", np_variance_random), ("MethylKit randomized", mk_variance_random), ("MethCP randomized", mcp_variance_random)]:
        pa.figure()
        plt.title(f"Length vs variance HG003 chr21 {caller_name}")
        x = np.log10(variance[0])
        idx = ~(np.isnan(x) | np.isinf(x))
        plot_2d_density(x[idx], variance[1][idx], cmap="jet")
        pa.savefig()
