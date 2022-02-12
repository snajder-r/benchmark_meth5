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
from meth5.meth5 import MetH5File, compute_betascore

from benchmark_pycometh.config import module_config
from benchmark_pycometh.utils import unions, count_cpgs_available

matplotlib.use("Agg")
chroms = [str(i) for i in range(1, 22)]

pa = PlotArchiver("pycometh", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

pm_hg003 = {
    "fs": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/asm_fastseg/met_comp_adj_HG003.bed",
        "name": "Fastseg (w. hp)",
    },
    "mk": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/asm_fastseg/met_comp_adj_HG003.bed",
        "name": "Methylkit",
    },
    "nes": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/asm/met_comp_adj_HG003.bed",
        "name": "Nanoepiseg",
    },
}

pm_parents = {
    "fs": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_methylkit/met_comp_adj_parents.bed",
        "name": "Fastseg (w. hp)",
    },
    "nes": {
        "filename": "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents/met_comp_adj_parents.bed",
        "name": "Nanoepiseg",
    },
    "mk_diff": {
        "filename": "/omics/groups/OE0540/internal/projects/nanopore/pycometh_benchmark/diffmet_parents_methylkit/methylkit_diff_met_comp.tsv",
        "name": "Methylkit diff",
    },
}

ref = pyfaidx.Fasta(module_config.reference, "r")

ref_cpgs = ReferenceCpGs(module_config.reference)


def annotate(hits, m5s):
    for hit in tqdm.tqdm(list(hits)):
        hit["CpGs"] = count_cpgs_available(hit["chrom"], hit["start"], hit["end"], m5s, ref_cpgs)
        hit["diff"] = eval(hit["difference"])[0]
        yield hit


def filter(hits):
    for hit in hits:
        if hit["chrom"] not in chroms:
            continue
        yield hit


sample_m5: Dict[str, MetH5File] = {
    sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r") for sample in module_config.samples
}

for tool in pm_hg003:
    pm_hg003[tool]["pm"] = PycomethOutput(pm_hg003[tool]["filename"])
    pm_hg003[tool]["hits"] = list(
        annotate(
            filter(pm_hg003[tool]["pm"].read_file(drop_insignificant=False, min_diff=0.5, pval_threshold=0.05)),
            [sample_m5["HG003"]],
        )
    )

print("CpGs found using nanoepiseg: ", len(unions(h["CpGs"] for h in pm_hg003["nes"]["hits"])))
print("CpGs found using methylkit: ", len(unions(h["CpGs"] for h in pm_hg003["mk"]["hits"])))
print("CpGs found using fastseg: ", len(unions(h["CpGs"] for h in pm_hg003["fs"]["hits"])))


print("Segments found using nanoepiseg: ", len(pm_hg003["nes"]["hits"]))
print("Segments found using methylkit: ", len(pm_hg003["mk"]["hits"]))
print("Segments found using fastseg: ", len(pm_hg003["fs"]["hits"]))


def difference(hits_a, hits_b):
    not_found = 0
    found = 0
    for hit_a in tqdm.tqdm(hits_a):
        cpgs = set(hit_a["CpGs"])
        for hit_b in hits_b:
            if hit_a["chrom"] != hit_b["chrom"]:
                continue
            othercpgs = set(hit_b["CpGs"])
            cpgs = cpgs.difference(othercpgs)
        
        found += len(hit_a["CpGs"]) - len(cpgs)
        not_found += len(cpgs)
    return not_found


for tool in pm_parents:
    pm_parents[tool]["pm"] = PycomethOutput(pm_parents[tool]["filename"])
    pm_parents[tool]["hits"] = list(
        annotate(
            filter(pm_parents[tool]["pm"].read_file(drop_insignificant=False, min_diff=0.5, pval_threshold=0.05)),
            [sample_m5["HG003"], sample_m5["HG004"]],
        )
    )
    
def load_methcp_result(path):
    with open(path, "r") as f:
        header = f.readline()
        for line in f.readlines():
            line = line.strip().split("\t")
            line = [c.replace('"',"") for c in line]
            hit = dict(chrom=line[1],
                        start = int(line[2]),
                        end = int(line[3]),
                        diff = float(line[6]),
                        pval = float(line[8])
                       )
            yield hit

        

print("CpGs found using nanoepiseg: ", len(unions(h["CpGs"] for h in pm_parents["nes"]["hits"])))
print("CpGs found using fastseg: ", len(unions(h["CpGs"] for h in pm_parents["fs"]["hits"])))
print("CpGs found using methylkit diff: ", len(unions(h["CpGs"] for h in pm_parents["mk_diff"]["hits"])))

print("Segments found using nanoepiseg: ", len(pm_parents["nes"]["hits"]))
print("Segments found using fastseg: ", len(pm_parents["fs"]["hits"]))
print("Segments found using methylkit diff: ", len(pm_parents["mk_diff"]["hits"]))


def difference_diffmet(hits_a, hits_b):
    not_found = []
    for hit_a in tqdm.tqdm(hits_a):
        cpgs = set(hit_a["CpGs"])
        for hit_b in hits_b:
            if hit_a["chrom"] != hit_b["chrom"]:
                continue
            othercpgs = set(hit_b["CpGs"])
            cpgs = cpgs.difference(othercpgs)
        if len(cpgs) == len(hit_a["CpGs"]):
            # entire segment was not found
            not_found.append(hit_a)
    return not_found


def plot_density_effect_vs_length(pm_list, xlim=[0, 3]):
    fig, axes = plt.subplots(1, len(pm_list), figsize=(18, 5))
    for i, tool in enumerate(pm_list):
        plt.sca(axes[i])
        plot_2d_density(
            np.log10([len(h["CpGs"]) for h in pm_list[tool]["hits"]]),
            np.array([h["diff"] for h in pm_list[tool]["hits"]]),
            cmap="jet",
        )
        axes[i].title.set_text(pm_list[tool]["name"])
        axes[i].set_xlim(xlim)
        axes[i].set_ylim([-1, 1])
    fig.text(0.5, 0.04, "Number of CpGs in segment (Log10)", ha="center")
    fig.text(0.04, 0.5, "Methylation rate difference", va="center", rotation="vertical")
    pa.savefig()


with pa.open_multipage_pdf("asm_hg003_comparison"):
    pa.figure()
    plt.title("Allele specific methylation in HG003")
    plt.bar(0, sum(1 for h in pm_hg003["nes"]["hits"]))
    plt.bar(1, sum(1 for h in pm_hg003["mk"]["hits"]))
    plt.bar(2, sum(1 for h in pm_hg003["fs"]["hits"]))
    plt.ylabel("Number of Segments")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure()
    plt.title("Allele specific methylation in HG003")
    plt.bar(0, sum(len(h["CpGs"]) for h in pm_hg003["nes"]["hits"]))
    plt.bar(1, sum(len(h["CpGs"]) for h in pm_hg003["mk"]["hits"]))
    plt.bar(2, sum(len(h["CpGs"]) for h in pm_hg003["fs"]["hits"]))
    plt.ylabel("Number of CpGs")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure()
    plt.title("Allele specific methylation in HG003")
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["nes"]["hits"]], positions=[0])
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["mk"]["hits"]], positions=[1])
    plt.violinplot([len(h["CpGs"]) for h in pm_hg003["fs"]["hits"]], positions=[2])
    plt.ylabel("Segment length")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit", "Fastseg (w. HP)"])
    pa.savefig()
    
    plot_density_effect_vs_length(pm_hg003, xlim=[0, 4])

with pa.open_multipage_pdf("diffmet_parents_comparison"):
    pa.figure()
    plt.title("Differential methylation HG003 vs HG004")
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["nes"]["hits"]], positions=[0])
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["mk_diff"]["hits"]], positions=[1])
    plt.violinplot([len(h["CpGs"]) for h in pm_parents["fs"]["hits"]], positions=[2])
    plt.ylabel("Segment length")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit Diffmet", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure()
    plt.title("Differential methylation HG003 vs HG004")
    plt.bar(0, sum(1 for h in pm_parents["nes"]["hits"]))
    plt.bar(1, sum(1 for h in pm_parents["mk_diff"]["hits"]))
    plt.bar(2, sum(1 for h in pm_parents["fs"]["hits"]))
    plt.ylabel("Number of Segments")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit Diffmet", "Fastseg (w. HP)"])
    pa.savefig()
    
    pa.figure()
    plt.title("Differential methylation HG003 vs HG004")
    plt.bar(0, sum(len(h["CpGs"]) for h in pm_parents["nes"]["hits"]))
    plt.bar(1, sum(len(h["CpGs"]) for h in pm_parents["mk_diff"]["hits"]))
    plt.bar(2, sum(len(h["CpGs"]) for h in pm_parents["fs"]["hits"]))
    plt.ylabel("Number of CpGs")
    plt.xticks([0, 1, 2], labels=["Nanoepiseg", "Methylkit Diffmet", "Fastseg (w. HP)"])
    pa.savefig()
    
    plot_density_effect_vs_length(pm_parents, xlim=[0, 3])


metrates = {
    sample: load_metrates(module_config.bs_seq_files[sample]["combined"], has_chr=True, has_percent=False, sort=False)
    for sample in module_config.bs_seq_files
}

metrates = {sample: metrates[sample].groupby("chrom") for sample in metrates}


def pair_hits_with_diffmet_from_bsseq(hits):
    metrates_temp = {
        sample: {chrom: metrates[sample].get_group(chrom) for chrom in metrates[sample].groups} for sample in metrates
    }
    for hit in tqdm.tqdm(hits):
        mean_met = {}
        for sample in metrates_temp:
            idx_before = metrates_temp[sample][hit["chrom"]]["start"] >= hit["start"]
            idx_after = metrates_temp[sample][hit["chrom"]]["end"] < hit["end"]
            mean_met[sample] = np.nanmean(metrates_temp[sample][hit["chrom"]].loc[idx_before & idx_after]["fracmet"])
            metrates_temp[sample][hit["chrom"]] = metrates_temp[sample][hit["chrom"]].loc[idx_before]
        yield hit["diff"], mean_met["HG004"] - mean_met["HG003"]


validate_nes = np.array(list(pair_hits_with_diffmet_from_bsseq(pm_parents["nes"]["hits"])))
validate_fs = np.array(list(pair_hits_with_diffmet_from_bsseq(pm_parents["fs"]["hits"])))
validate_mk_diff = np.array(list(pair_hits_with_diffmet_from_bsseq(pm_parents["mk_diff"]["hits"])))

with pa.open_multipage_pdf("bs_validation_parents_diff"):
    pa.figure()
    plt.title("Nanoepiseg + Pycometh diffmet bisulfite validation")
    plt.scatter(validate_nes[:, 0], validate_nes[:, 1], s=1)
    plt.ylabel("BS-Seq")
    plt.ylabel("ONT")
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    pa.savefig()
    
    pa.figure()
    plt.title("Fastseg + Pycometh diffmet bisulfite validation")
    plt.scatter(validate_fs[:, 0], validate_fs[:, 1], s=1)
    plt.ylabel("BS-Seq")
    plt.ylabel("ONT")
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    pa.savefig()
    
    pa.figure()
    plt.title("Methylkit diff segmentation")
    plt.scatter(validate_mk_diff[:, 0], validate_mk_diff[:, 1], s=1)
    plt.ylabel("BS-Seq")
    plt.ylabel("ONT")
    plt.ylim(-1, 1)
    plt.xlim(-1, 1)
    pa.savefig()


from nanoepitools.annotations.annotations import GFFAnnotationsReader

gff = GFFAnnotationsReader()
gff.read(module_config.gff_file, only_protein_coding=False)
gff.build_index()

promoters_nes = pm_parents["nes"]["pm"].load_promoters_hit(gff, 2000, 500, min_diff=0.5)
promoters_mk = pm_parents["mk_diff"]["pm"].load_promoters_hit(gff, 2000, 500, min_diff=0.5)

difference(pm_parents["nes"]["hits"], pm_parents["mk_diff"]["hits"])
difference(pm_parents["mk_diff"]["hits"], pm_parents["nes"]["hits"])


"""
SANITY CHECK
"""

hp_lookup = {s: {h: int(i) for i, h in sample_m5[s].h5_fp["reads/read_groups/haplotype"].attrs.items()} for s in sample_m5}
for hit in pm_hg003["mk"]["hits"]:
    agg = (
        sample_m5["HG003"][hit["chrom"]]
        .get_values_in_range(hit["start"], hit["end"])
        .get_llr_site_readgroup_aggregate(compute_betascore, "haplotype")
    )
    rates = {hp: np.nanmean(agg[hp][0]) for hp in agg}
    if all(hp_lookup["HG003"][hp] in rates for hp in ("H1", "H2")):
        diff_diff = hit["diff"] - (rates[hp_lookup["HG003"]["H1"]] - rates[hp_lookup["HG003"]["H2"]])
        if abs(diff_diff) > 0.25:
            print(hit["diff"], (rates[hp_lookup["HG003"]["H1"]] - rates[hp_lookup["HG003"]["H2"]]))

def ratio_overmethylated(hits):
    return sum(hit["diff"] > 0 for hit in hits) / len(hits)

ratio_overmethylated(pm_hg003["nes"]["hits"])
