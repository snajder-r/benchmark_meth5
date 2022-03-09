import tqdm
import pandas as pd
import numpy as np
from meth5 import MetH5File
from benchmark_pycometh.config import module_config
from nanoepitools.pycometh_result import PycomethOutput

mf = {
    sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r") for sample in module_config.samples
}
read_groups = {sample: mf[sample].get_all_read_groups("haplotype") for sample in mf}
parents_diffmet_file = (
    "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_bs_diff/met_comp_adj_parents.bed"
)

pm = PycomethOutput(parents_diffmet_file)
hits = list(pm.read_file(drop_insignificant=False, min_diff=0.25, pval_threshold=0.05))

mapping_hps = pd.read_csv(module_config.haplotype_mapping_file, sep="\t", dtype={"chrom": str, "child_hp": str})
mapping_hps["child_hp"] = mapping_hps["child_hp"].map(lambda x: f"H{x}")
mapping_hps["parent_hp"] = mapping_hps["parent_hp"].map(lambda x: f"H{x}")
mapping_hps = mapping_hps.groupby(["chrom", "child_hp"])


def compute_hit_betascores(hit, mean=False, merge_hp=False):
    # Compute beta scores
    hp_rates = {
        sample: mf[sample][hit["chrom"]]
        .get_values_in_range(hit["start"], hit["end"])
        .get_llr_site_readgroup_rate("haplotype")
        for sample in mf
    }
    # translate haplotype names
    hp_rates = {s: {read_groups[s][hp]: r for hp, r in sr.items()} for s, sr in hp_rates.items()}
    # Remove unknown haplotypes
    hp_rates = {s: {hp: r for hp, r in sr.items() if hp in {"H1", "H2"}} for s, sr in hp_rates.items()}
    if mean:
        if merge_hp:
            hp_rates = {s: np.nanmean([ri for r in sr.values() for ri in r[0]]) for s, sr in hp_rates.items()}
        else:
            hp_rates = {s: {hp: np.nanmean(r[0]) for hp, r in sr.items()} for s, sr in hp_rates.items()}
    return hp_rates


hit = hits[1]


def get_rates_by_source(rates, coords, child_hp):
    segments = mapping_hps.get_group((hit["chrom"], child_hp))
    segments = segments.loc[(segments["start"] < coords[-1, -1]) & (coords[0, 0] <= segments["end"])]
    ret = {}
    for r, coord in zip(rates, coords):
        coord_segment = segments.loc[(segments["start"] < coord[-1]) & (coord[0] <= segments["end"])]
        key = (coord_segment.iloc[0]["parent"], coord_segment.iloc[0]["parent_hp"])
        if key not in ret:
            ret[key] = []
        ret[key].append(r)
    ret = {k: np.nanmean(v) for k, v in ret.items()}
    return ret


def generate_paired_rates(hit):
    hp_rates = compute_hit_betascores(hit)
    sample_hp_source_rates = {
        s: {
            hp: {child_hp: get_rates_by_source(*hp_rates[s][hp], child_hp) for child_hp in hp_rates["HG002"].keys()}
            for hp in hp_rates[s]
        }
        for s in hp_rates
    }
    
    for child_hp in sample_hp_source_rates["HG002"]:
        for source_key, bs in sample_hp_source_rates["HG002"][child_hp][child_hp].items():
            parent, parent_hp = source_key
            if parent not in sample_hp_source_rates:
                continue
            if parent_hp not in sample_hp_source_rates[parent]:
                continue
            parent_rates = sample_hp_source_rates[parent][parent_hp][child_hp]
            if source_key not in parent_rates:
                continue
            matched_rate = parent_rates[source_key]
            if np.isnan(bs) or np.isnan(matched_rate):
                continue
            yield bs, matched_rate


    
    
    

all_pairs = []
hg003_pairs = []
hg004_pairs = []

import scipy.stats
for hit in hits:
    for p in generate_paired_rates(hit):
        all_pairs.append(p)
    if len(all_pairs) > 3:
        print("Matched", scipy.stats.pearsonr([x[0] for x in all_pairs], [x[1] for x in all_pairs]))

    simple_bs = compute_hit_betascores(hit, mean=True, merge_hp=True)
    if "HG002" in simple_bs and "HG003" in simple_bs:
        hg003_pairs.append([simple_bs["HG002"], simple_bs["HG003"]])
    if "HG002" in simple_bs and "HG004" in simple_bs:
        hg004_pairs.append([simple_bs["HG002"], simple_bs["HG004"]])
    if len(hg003_pairs) > 3:
        print("HG003 cor", scipy.stats.pearsonr([x[0] for x in hg003_pairs], [x[1] for x in hg003_pairs]))
    if len(hg004_pairs) > 3:
        print("HG004 cor", scipy.stats.pearsonr([x[0] for x in hg004_pairs], [x[1] for x in hg004_pairs]))