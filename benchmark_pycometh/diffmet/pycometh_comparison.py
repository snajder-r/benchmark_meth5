from typing import Dict
import pyfaidx
import numpy as np
import tqdm
import matplotlib
import matplotlib.pyplot as plt
from nanoepitools.pycometh_result import PycomethOutput
from benchmark_pycometh.bsseq.bsseq import load_metrates
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density
from nanoepitools.reference_cpgs import ReferenceCpGs
from nanoepitools.math import fdr_from_pvals
from meth5.meth5 import MetH5File, compute_betascore

from benchmark_pycometh.config import module_config
from benchmark_pycometh.utils import unions, count_cpgs_available

matplotlib.use("Agg")
chroms = [str(i) for i in range(1, 22)]

pa = PlotArchiver("pycometh", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})


def pycometh_loader(filename, filter_singlecall=True):
    for hit in PycomethOutput(filename).read_file(
        drop_insignificant=False, min_diff=0.5, pval_threshold=0.05, use_raw_pvalue=True
    ):
        if filter_singlecall:
            sample_rates = {
                sample: sample_m5[sample][hit["chrom"]].get_values_in_range(hit["start"], hit["end"]).get_llr_site_rate()[0]
                for sample in ["HG003", "HG004"]}
            if any((~np.isnan(sample_rates[sample])).sum() <= 1 for sample in sample_rates):
                print("Skipping because based on a single call")
                continue
        yield hit


def load_methcp_result(path, min_diff=0.5):
    hits = []
    
    with open(path, "r") as f:
        header = f.readline()
        for line in f.readlines():
            line = line.strip().split("\t")
            line = [c.replace('"', "") for c in line]
            hit = dict(
                chrom=line[1], start=int(line[2]), end=int(line[3]), pval=float(line[8]), difference=-float(line[6])
            )
            hits.append(hit)
    pvals = np.array([hit["pval"] for hit in hits])
    diffs = np.array([hit["difference"] for hit in hits])
    idx = (pvals < 0.05) & (np.abs(diffs) > min_diff)
    for sig, hit in tqdm.tqdm(list(zip(idx, hits))):
        if sig:
            sample_rates = {
                sample: sample_m5[sample][hit["chrom"]]
                .get_values_in_range(hit["start"], hit["end"])
                .get_llr_site_rate()[0]
                for sample in ["HG003", "HG004"]
            }
            if any((~np.isnan(sample_rates[sample])).sum() <= 1 for sample in sample_rates):
                print("Skipping because based on a single call")
                continue
            hit["difference"] = np.nanmean(sample_rates["HG004"]) - np.nanmean(sample_rates["HG003"])
            if abs(hit["difference"]) > min_diff:
                yield hit


pm_parents = {
    "fs": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_methylkit/met_comp_adj_parents.bed",
        "name": "Fastseg (w. hp)",
        "loader": lambda x: pycometh_loader(x, False),
    },
    "nes": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_bs_diff/met_comp_adj_parents.bed",
        "name": "Nanoepiseg (BS diff)",
        "loader": lambda x: pycometh_loader(x, False),
    },
    "nes_llr": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_llr_diff/met_comp_adj_parents.bed",
        "name": "Nanoepiseg (LLR diff)",
        "loader": lambda x: pycometh_loader(x, True),
    },
    "metcp": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_methylkit/sig_diff_methcp.tsv",
        "name": "MethCP",
        "loader": load_methcp_result,
    },
}

ref = pyfaidx.Fasta(module_config.reference, "r")

ref_cpgs = ReferenceCpGs(module_config.reference)


def annotate(hits, m5s):
    for hit in tqdm.tqdm(list(hits)):
        hit["CpGs"] = count_cpgs_available(hit["chrom"], hit["start"], hit["end"], m5s, ref_cpgs)
        if len(hit["CpGs"]) == 0:
            print(hit)
        hit["diff"] = eval(hit["difference"])[0] if isinstance(hit["difference"], str) else hit["difference"]
        yield hit


def filter(hits):
    for hit in hits:
        if hit["chrom"] not in chroms:
            continue
        yield hit


sample_m5: Dict[str, MetH5File] = {
    sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r")
    for sample in ["HG003", "HG004"]  # module_config.samples
}


for tool in pm_parents:
    pm_parents[tool]["hits"] = list(
        annotate(
            filter(pm_parents[tool]["loader"](pm_parents[tool]["filename"])),
            [sample_m5["HG003"], sample_m5["HG004"]],
        )
    )


print("CpGs found using nanoepiseg: ", len(unions(h["CpGs"] for h in pm_parents["nes"]["hits"])))
print("CpGs found using nanoepiseg (LLR diff): ", len(unions(h["CpGs"] for h in pm_parents["nes_llr"]["hits"])))
print("CpGs found using fastseg: ", len(unions(h["CpGs"] for h in pm_parents["fs"]["hits"])))
print("CpGs found using MethCP: ", len(unions(h["CpGs"] for h in pm_parents["metcp"]["hits"])))

print("Segments found using nanoepiseg: ", len(pm_parents["nes"]["hits"]))
print("Segments found using nanoepiseg (LLR diff): ", len(pm_parents["nes"]["hits"]))
print("Segments found using fastseg: ", len(pm_parents["fs"]["hits"]))
print("Segments found using MethCP: ", len(pm_parents["metcp"]["hits"]))

##########################################################


pm_hg003 = {
    "mk": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/asm_fastseg/met_comp_adj_HG003.bed",
        "name": "Methylkit", "loader": lambda x: pycometh_loader(x, False),
    },
    "nes": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/asm/met_comp_adj_HG003.bed",
        "name": "Nanoepiseg", "loader": lambda x: pycometh_loader(x, False),
    },
}

for tool in pm_hg003:
    pm_hg003[tool]["pm"] = PycomethOutput(pm_hg003[tool]["filename"])
    pm_hg003[tool]["hits"] = list(
        annotate(
            filter(pm_hg003[tool]["loader"](pm_hg003[tool]["filename"])),
            [sample_m5["HG003"]],
        )
    )

print("CpGs found using nanoepiseg: ", len(unions(h["CpGs"] for h in pm_hg003["nes"]["hits"])))
print("CpGs found using methylkit: ", len(unions(h["CpGs"] for h in pm_hg003["mk"]["hits"])))

print("Segments found using nanoepiseg: ", len(pm_hg003["nes"]["hits"]))
print("Segments found using methylkit: ", len(pm_hg003["mk"]["hits"]))


###############################################################


def plot_density_effect_vs_length(pm_list, xlim=[0, 3]):
    fig, axes = plt.subplots(1, len(pm_list), figsize=(18, 5))
    for i, tool in enumerate(pm_list):
        plt.sca(axes[i])
        hits = pm_list[tool]["hits"]
        numcpgs = [len(h["CpGs"]) for h in hits]
        idx = [n>0 for n in numcpgs]
        plot_2d_density(
            np.log10([n for use, n in zip(idx, numcpgs) if use]),
            np.array([h["diff"] for use, h in zip(idx, hits) if use]),
            cmap="jet",
        )
        axes[i].title.set_text(pm_list[tool]["name"])
        axes[i].set_xlim(xlim)
        axes[i].set_ylim([-1, 1])
    fig.text(0.5, 0.04, "Number of CpGs in segment (Log10)", ha="center")
    fig.text(0.04, 0.5, "Methylation rate difference", va="center", rotation="vertical")
    pa.savefig()


with pa.open_multipage_pdf("asm_hg003_comparison"):
    figsize=(6,2)
    pa.figure(figsize=figsize)
    plt.title("Allele specific methylation in HG003")
    plt.barh(0, sum(1 for h in pm_hg003["nes"]["hits"]))
    plt.barh(1, sum(1 for h in pm_hg003["mk"]["hits"]))
    plt.barh(2, sum(1 for h in pm_hg003["fs"]["hits"]))
    plt.xlabel("Number of Segments")
    plt.yticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure(figsize=figsize)
    plt.title("Allele specific methylation in HG003")
    plt.barh(0, sum(len(h["CpGs"]) for h in pm_hg003["nes"]["hits"]))
    plt.barh(1, sum(len(h["CpGs"]) for h in pm_hg003["mk"]["hits"]))
    plt.barh(2, sum(len(h["CpGs"]) for h in pm_hg003["fs"]["hits"]))
    plt.xlabel("Number of CpGs")
    plt.yticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure(figsize=figsize)
    plt.title("Allele specific methylation in HG003")
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["nes"]["hits"]], positions=[0])
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["mk"]["hits"]], positions=[1])
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["fs"]["hits"]], positions=[2])
    plt.ylabel("Segment length")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    plot_density_effect_vs_length(pm_hg003, xlim=[0, 4])

with pa.open_multipage_pdf("diffmet_parents_comparison"):
    figsize=(6,2)
    pa.figure(figsize=figsize)
    plt.title("Differential methylation HG003 vs HG004")
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["nes"]["hits"]], positions=[0])
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["metcp"]["hits"]], positions=[1])
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["fs"]["hits"]], positions=[2])
    plt.ylabel("Segment length")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "MethCP", "Fastseg (w. HP)"])
    pa.savefig()

    figsize = (6, 2)
    pa.figure(figsize=figsize)
    plt.title("Differential methylation HG003 vs HG004")
    plt.violinplot([abs(h["diff"]) for h in pm_parents["nes"]["hits"]], positions=[0])
    plt.violinplot([abs(h["diff"]) for h in pm_parents["metcp"]["hits"]], positions=[1])
    plt.violinplot([abs(h["diff"]) for h in pm_parents["fs"]["hits"]], positions=[2])
    plt.ylabel("Effect size")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "MethCP", "Fastseg (w. HP)"])
    pa.savefig()

    pa.figure(figsize=figsize)
    plt.title("Differential methylation HG003 vs HG004")
    plt.barh(0, sum(1 for h in pm_parents["nes"]["hits"]))
    plt.barh(1, sum(1 for h in pm_parents["metcp"]["hits"]))
    plt.barh(2, sum(1 for h in pm_parents["fs"]["hits"]))
    plt.xlabel("Number of Segments")
    plt.yticks([0, 1, 2], labels=["Nanoepiseg", "MethCP", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure(figsize=figsize)
    plt.title("Differential methylation HG003 vs HG004")
    plt.barh(0, sum(len(h["CpGs"]) for h in pm_parents["nes"]["hits"]))
    plt.barh(1, sum(len(h["CpGs"]) for h in pm_parents["metcp"]["hits"]))
    plt.barh(2, sum(len(h["CpGs"]) for h in pm_parents["fs"]["hits"]))
    plt.xlabel("Number of CpGs")
    plt.yticks([0, 1, 2], labels=["Nanoepiseg", "MethCP", "Fastseg (w. HP)"])
    pa.savefig()
    
    plot_density_effect_vs_length(pm_parents, xlim=[0, 3])
