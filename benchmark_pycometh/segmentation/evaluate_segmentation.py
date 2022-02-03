import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from meth5.meth5 import MetH5File
import numpy as np
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density

from benchmark_pycometh.config import module_config
import math

pa = PlotArchiver("segmentation", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})


def subsampled_segments(segment_table, subsample=1.0):
    assert (
        segment_table.columns[0] == "chrom"
        and segment_table.columns[1] == "start"
        and segment_table.columns[2] == "end"
    )
    for segment in segment_table.itertuples(index=False):
        if subsample == 1.0:
            num_chunks = 1
        else:
            num_chunks = int(math.ceil((1 / subsample) - np.random.rand()))
        subseg_len = (segment[2] - segment[1]) // num_chunks
        for i in range(num_chunks):
            yield {
                "chrom": segment[0],
                "start": int(segment[1] + subseg_len * i),
                "end": int(segment[2] + subseg_len * (i + 1)),
            }


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

subsample_rate = seg_methylkit.shape[0] / seg_nanoepiseg.shape[0]

mf = MetH5File(module_config.meth5_template_file.format(sample="HG003"), "r")
hp_dict = {int(k): v for k, v in mf.h5_fp["reads/read_groups/haplotype"].attrs.items()}
hp_ids = [k for k, v in hp_dict.items() if v in {"H1", "H2"}]




def test_idx(table):
    return table.loc[(table["chrom"] == "5") & (table["start"] <= 6051182) & (5951256 < table["end"])]

def test_idx(table):
    return table.loc[(table["chrom"] == "21")]


def compute_subseg_stdevs(*args, **kwargs):
    i = 0
    hp_dict = {int(k): v for k, v in mf.h5_fp["reads/read_groups/haplotype"].attrs.items()}
    for segment in subsampled_segments(*args, **kwargs):
        
        rg_agg = (
            mf[segment["chrom"]]
            .get_values_in_range(segment["start"], segment["end"])
            .get_llr_site_readgroup_rate("haplotype")
        )
        for hp, (bs_segment, ranges_segment) in rg_agg.items():
            if hp_dict.get(hp, "none") not in {"H1", "H2"}:
                continue
            mean_bs_segment = np.nanmean(bs_segment)
            std_bs_segment = np.nanstd(bs_segment)
            if np.isnan(mean_bs_segment):
                continue
            
            num_subsegments = len(bs_segment) // 10
            if num_subsegments < 2:
                continue
            subsegment_starts = list(range(0, len(bs_segment), len(bs_segment) // num_subsegments)) + [
                len(bs_segment) - 1
            ]
            
            diff_bs_segment = []
            for subseg_start, subseg_end in zip(subsegment_starts[:-1], subsegment_starts[1:]):
                bs_subsegment = np.nanmean(bs_segment[subseg_start:subseg_end])
                if np.isnan(bs_subsegment):
                    continue
                
                diff_bs_segment.append((bs_subsegment - mean_bs_segment) ** 2)
            yield std_bs_segment, np.sqrt(np.mean(diff_bs_segment))

stdevs_nes = list(compute_subseg_stdevs(test_idx(seg_nanoepiseg)))
stdevs_mk = list(compute_subseg_stdevs(test_idx(seg_methylkit), subsample=subsample_rate))


stdevs_nes = np.array(stdevs_nes)
stdevs_mk = np.array(stdevs_mk)

with pa.open_multipage_pdf("segmentation_purity"):
    pa.figure()
    plot_2d_density(stdevs_nes[:,0], stdevs_nes[:,1], cmap="jets")
    plt.ylim(0,0.5)
    plt.xlim(0, 0.5)
    pa.savefig()
    pa.figure()
    plot_2d_density(stdevs_mk[:, 0], stdevs_mk[:, 1], cmap="jets")
    plt.ylim(0, 0.5)
    plt.xlim(0, 0.5)
    pa.savefig()
