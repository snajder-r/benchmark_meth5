import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from benchmark_pycometh.plotter import Plotter
from nanoepitools.plotting.general_plotting import PlotArchiver
from benchmark_pycometh.bsseq.bsseq import load_metrates
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from benchmark_pycometh.config import module_config

matplotlib.use("Agg")

pa = PlotArchiver("segmentation", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

gff = GFFAnnotationsReader()
gff.read(module_config.gff_file, only_protein_coding=False)
gff.build_index()

seg_methylkit = pd.read_csv(
    module_config.methylkit_segmentation["parents_mock_from_np_diffmet"],
    sep="\t",
    usecols=[0, 1, 2],
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
    skiprows=1,
)
seg_methylkit["chrom"] = seg_methylkit["chrom"].map(lambda x: x.replace("chr", ""))

seg_nanoepiseg = pd.read_csv(
    module_config.nanoepiseg_segmentation["parents"],
    sep="\t",
    names=["chrom", "start", "end"],
    dtype={"chrom": "str"},
)

chroms = set(seg_nanoepiseg["chrom"]).intersection(set(seg_methylkit["chrom"]))
seg_methylkit = seg_methylkit.loc[seg_methylkit.apply(lambda x: x["chrom"] in chroms, axis=1)]
seg_nanoepiseg = seg_nanoepiseg.loc[seg_nanoepiseg.apply(lambda x: x["chrom"] in chroms, axis=1)]

pl = Plotter(gff, pa)




metrates = {
    "HG003_R1": load_metrates(module_config.bs_seq_files["HG003"]["R1"], True, True),
    "HG003_pseudobulk": load_metrates(module_config.mock_bsseq_template_file.format(sample="HG003"), False, False),
    "HG003_combined": load_metrates(module_config.bs_seq_files["HG003"]["combined"], False, True),
}


def plot_segmentation_comparison(chrom, start, end, title="", **kwargs):
    figure_kwargs = {"figsize": (12, 4)}
    
    def plot(seg, segtype):
        pl.plot_region(
            chrom,
            start,
            end,
            figure_kwargs=figure_kwargs,
            ws=0,
            title=f"{title} ({segtype})",
            aggregate_samples=False,
            with_no_hp=True,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=False,
            hold=False,
            vlines_in_coord_space=seg,
            vlinewidth=1,
            **kwargs
        )
    
    seg_m = seg_methylkit.loc[seg_methylkit["chrom"] == chrom]
    seg_m = seg_m.loc[(seg_m["end"] > start) & (end > seg_m["start"])]
    plot(list(seg_m["end"]), "MK")
    seg_n = seg_nanoepiseg.loc[seg_nanoepiseg["chrom"] == chrom]
    seg_n = seg_n.loc[(seg_n["end"] > start) & (end > seg_n["start"])]
    plot(list(seg_n["end"]), "NES")


with pa.open_multipage_pdf("test"):
    chrom = "9"
    start = 34990414
    end = 34990823
    plot_segmentation_comparison(chrom, start-4000, end+4000, title="Not found by NES",
                                 annotations=[{"region":(start,end), "text":"Diffmet NES", "color":"r"}])
    
    for s in metrates:
        pa.figure(figsize=(20, 4))
        idx = (metrates[s]["chrom"] == chrom) & (metrates[s]["start"] >= start) & (metrates[s]["end"] <= end)
        bs_temp = metrates[s].loc[idx]
        bs_temp = bs_temp.sort_values(["chrom", "start", "end"])
        plt.plot(bs_temp["start"], bs_temp["fracmet"], label=s)
        plt.legend()
        pa.savefig()
    
    # plot_segmentation_comparison("13", 57630000 - 50000, 57634998 + 50000, title="PCDH17")
    # plot_segmentation_comparison("9", 95498525 - 50000, 95508336 + 50000, title="PTCH1")

metrates = {sample: load_metrates(module_config.bs_seq_files[sample]["combined"]) for sample in module_config.bs_seq_files}

with pa.open_multipage_pdf("test"):
    figure_kwargs = {"figsize": (20, 4)}
    pl.plot_region(chrom, 40076150+20000, 40076150+30000, figure_kwargs=figure_kwargs, ws=0, title=f"Long diffmet region",
        aggregate_samples=True, with_no_hp=False, coordinate_space=False, marker_height=0.9, fontsize=8,
        show_has_value_indicator=False, hold=False, vlinewidth=3)
    pa.savefig()