import tqdm
import pandas as pd
import numpy as np
import scipy.stats
from meth5 import MetH5File
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch, Circle
from benchmark_pycometh.config import module_config
import matplotlib.pyplot as plt
from nanoepitools.pycometh_result import PycomethOutput
from nanoepitools.plotting.plot_methylation_profile import default_color_map
from nanoepitools.math import p_to_llr


def load_hp_parent_child_mapping():
    mapping_hps = pd.read_csv(module_config.haplotype_mapping_file, sep="\t", dtype={"chrom": str, "child_hp": str})
    mapping_hps["child_hp"] = mapping_hps["child_hp"].map(lambda x: f"H{x}")
    mapping_hps["parent_hp"] = mapping_hps["parent_hp"].map(lambda x: f"H{x}")
    mapping_hps = mapping_hps.groupby(["chrom", "child_hp"])
    return mapping_hps


def get_segments(mapping_hps, chrom, start, end, child_hp, **kwargs):
    key = (chrom, child_hp)
    if key not in mapping_hps.groups:
        return None
    segments = mapping_hps.get_group(key)
    segments = segments.loc[(segments["start"] < end) & (start <= segments["end"])]
    return segments


def create_faux_child(mapping_hps, hp_rates, chrom, start, end, return_sources=False, **kwargs):
    faux_rates = {k: [[np.zeros(0)], [np.zeros((0, 2))]] for k in ["H1", "H2"]}
    faux_sources = {k: {"parent": [], "parent_hp": []} for k in faux_rates}
    total_segments = 0
    for child_hp in faux_rates:
        segments = get_segments(mapping_hps, chrom, start, end, child_hp)
        if segments is not None:
            for segment in segments.itertuples():
                seg_start = max(start, segment.start)
                seg_end = min(end, segment.end)
                source_rates, source_coords = hp_rates.get(segment.parent, {}).get(
                    segment.parent_hp, [np.zeros(0), np.zeros((0, 2))]
                )
                idx = [c[0] < seg_end and seg_start < c[1] for c in source_coords]
                if sum(idx) == 0:
                    continue
                total_segments += 1
                source_rates = source_rates[idx]
                source_coords = source_coords[idx, :]
                faux_rates[child_hp][0].append(source_rates)
                faux_rates[child_hp][1].append(source_coords)
                faux_sources[child_hp]["parent"] += [segment.parent] * len(source_rates)
                faux_sources[child_hp]["parent_hp"] += [segment.parent_hp] * len(source_rates)
        faux_rates[child_hp][0] = np.concatenate(faux_rates[child_hp][0])
        faux_rates[child_hp][1] = np.concatenate(faux_rates[child_hp][1])
    
    # Decide duplicates (calls that spanned breakpoints)
    for child_hp in faux_rates:
        faux_sources[child_hp]["parent"] = np.array(faux_sources[child_hp]["parent"])
        faux_sources[child_hp]["parent_hp"] = np.array(faux_sources[child_hp]["parent_hp"])
        
        coords = faux_rates[child_hp][1][:, 0]
        if len(coords) == 0:
            continue
        idx_unique = np.array([sum(coords == p) == 1 for p in coords])
        faux_rates[child_hp][1] = faux_rates[child_hp][1][idx_unique, :]
        faux_rates[child_hp][0] = faux_rates[child_hp][0][idx_unique]
        faux_sources[child_hp]["parent"] = faux_sources[child_hp]["parent"][idx_unique]
        faux_sources[child_hp]["parent_hp"] = faux_sources[child_hp]["parent_hp"][idx_unique]
    
    if return_sources:
        return faux_rates, faux_sources
    return faux_rates


def compute_hit_betascores(mf, read_groups, chrom, start, end, mean=False, merge_hp=False, **kwargs):
    # Compute beta scores
    hp_rates = {
        sample: mf[sample][chrom].get_values_in_range(start, end).get_llr_site_readgroup_rate("haplotype")
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


def color(p):
    if np.isnan(p):
        return (255, 255, 255, 0)
    return default_color_map(p_to_llr(p))


def plot_met_patches(x, x_end, y, color):
    patches = [Rectangle((x[i], y[i]), x_end[i] - x[i], 1) for i in range(len(x))]
    patch_collection = PatchCollection(patches)
    patch_collection.set_color(color)
    patch_collection.set_edgecolor(None)
    plt.gca().add_collection(patch_collection)


def plot_faux_hg002(df):
    pa.figure(figsize=(100, 2))
    df = df.copy()
    df["chrom"] = df.index.map(lambda x: x.split("_")[0])
    df["pos"] = df.index.map(lambda x: int(x.split("_")[1]))
    df = df.sort_values(["chrom", "pos"]).groupby("chrom")
    x_off = 0
    samples = ["HG002_H1", "faux_HG002_H1", "HG002_H2", "faux_HG002_H2", "HG003_H1", "HG003_H2", "HG004_H1", "HG004_H2"]
    vlines_segments = []
    vlines_recombination = []
    for chrom in tqdm.tqdm(df.groups):
        chrom_df = df.get_group(chrom)
        x = x_off + np.arange(chrom_df.shape[0])
        for y, s in enumerate(samples):
            # plt.scatter(x, [y]*len(x), color=[color(p) for p in chrom_df[s]], s=1, marker="s")
            plot_met_patches(x, x + 0.8, [y] * len(x), [color(p) for p in chrom_df[s]])
        
        last_row = None
        
        x = x_off
        
        for row in chrom_df.itertuples():
            x += 1
            if last_row is not None:
                if row.segment != last_row.segment:
                    vlines_segments.append(x - 0.05)
                    plt.text(x - 0.1, -1, row.segment, fontsize=2)
                if last_row.H1_parent != row.H1_parent or last_row.H1_parent_hp != row.H1_parent_hp or \
                        last_row.H2_parent != row.H2_parent or last_row.H2_parent_hp != row.H2_parent_hp:
                    vlines_recombination.append(x - 0.1)
            last_row = row
        x_off = x
    
    plt.yticks(np.arange(len(samples)) + 0.5, samples)
    plt.gca().autoscale_view()
    plt.xlim(0, x_off)
    plt.vlines(vlines_recombination, plt.ylim()[0], plt.ylim()[1], linewidth=0.1, color="r")
    plt.vlines(vlines_segments, plt.ylim()[0], plt.ylim()[1], linewidth=0.2, color="k")
    pa.savefig()

if __name__ == "__main__":
    mf = {
        sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r")
        for sample in module_config.samples
    }
    
    parents_diffmet_file = (
        "/home/r933r/data/projects/nanopore/pycometh_benchmark/diffmet_parents_bs_diff/met_comp_adj_parents.bed"
    )
    
    pm = PycomethOutput(parents_diffmet_file)
    hits = list(pm.read_file(drop_insignificant=False, min_diff=0.5, pval_threshold=0.05))
    
    mapping_hps = load_hp_parent_child_mapping()
    read_groups = {sample: mf[sample].get_all_read_groups("haplotype") for sample in mf}
    hp_mean_rates_in_hits = {"HG002": [], "HG003": [], "HG004": [], "faux_HG002": []}
    for hit in tqdm.tqdm(hits):
        hp_rates = compute_hit_betascores(mf, read_groups, **hit)
        faux_child = create_faux_child(mapping_hps, hp_rates, **hit)
        hp_rates["faux_HG002"] = faux_child
        for sample in hp_mean_rates_in_hits:
            for hp in ["H1", "H2"]:
                rate = np.nanmean(hp_rates.get(sample, {}).get(hp, [np.nan])[0])
                hp_mean_rates_in_hits[sample].append(rate)
    
    hp_mean_rates_in_hits = {s: np.array(hp_mean_rates_in_hits[s]) for s in hp_mean_rates_in_hits}
    
    pa = PlotArchiver("inheritance", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})
    
    with pa.open_multipage_pdf("faux_HG002"):
        for sample in ["faux_HG002", "HG003", "HG004"]:
            idx = (~np.isnan(hp_mean_rates_in_hits[sample])) & (~np.isnan(hp_mean_rates_in_hits["HG002"]))
            print(
                sample,
                sum(idx),
                scipy.stats.spearmanr(hp_mean_rates_in_hits[sample][idx], hp_mean_rates_in_hits["HG002"][idx]),
            )
            
            pa.figure(figsize=(5, 5))
            plt.title(sample)
            plt.ylabel("HG002")
            plt.xlabel(sample)
            plot_2d_density(hp_mean_rates_in_hits[sample][idx], hp_mean_rates_in_hits["HG002"][idx], cmap="jet")
            pa.savefig()
    
    """ A very similar comparison, but this time on CpG-level """
    all_rates = {
        "HG002_H1": {},
        "HG003_H1": {},
        "HG004_H1": {},
        "faux_HG002_H1": {},
        "HG002_H2": {},
        "HG003_H2": {},
        "HG004_H2": {},
        "faux_HG002_H2": {},
    }
    all_sources = {k: {"parent": {}, "parent_hp": {}} for k in ["H1", "H2"]}
    all_segments = {}

    for segment, hit in enumerate(tqdm.tqdm(hits)):
        hp_rates = compute_hit_betascores(mf, read_groups, **hit)
        faux_child, sources = create_faux_child(mapping_hps, hp_rates, **hit, return_sources=True)
        hp_rates["faux_HG002"] = faux_child
        
        for sample in hp_rates:
            for hp in ["H1", "H2"]:
                rates_coords = hp_rates.get(sample, {}).get(hp, None)
                if rates_coords is None:
                    continue
                if len(rates_coords[1]) == 0:
                    continue
                key = f"{sample}_{hp}"
                for rate, coord, parent, parent_hp in zip(
                    *rates_coords, sources[hp]["parent"], sources[hp]["parent_hp"]
                ):
                    coord_key = f"{hit['chrom']}_{int(coord[0])}"
                    all_segments[coord_key] = segment
                    all_rates[key][coord_key] = rate
                    if sample == "faux_HG002":
                        all_sources[hp]["parent"][coord_key] = parent
                        all_sources[hp]["parent_hp"][coord_key] = parent_hp

    all_rates_df = pd.DataFrame(index=set(all_segments.keys()))
    for sample in all_rates:
        all_rates_df[sample] = pd.Series(all_rates[sample])
    
    for hp in ["H1", "H2"]:
        all_rates_df[f"{hp}_parent"] = pd.Series(all_sources[hp]["parent"])
        all_rates_df[f"{hp}_parent_hp"] = pd.Series(all_sources[hp]["parent_hp"])
    all_rates_df["segment"] = pd.Series(all_segments)
    
    with pa.open_multipage_pdf("faux_HG002_cglevel"):
        
        for sample in ["faux_HG002", "HG003", "HG004"]:
            idxh1 = (~np.isnan(all_rates_df[f"{sample}_H1"])) & (~np.isnan(all_rates_df[f"HG002_H1"]))
            idxh2 = (~np.isnan(all_rates_df[f"{sample}_H2"])) & (~np.isnan(all_rates_df[f"HG002_H2"]))
            
            vals_child = np.array(all_rates_df[f"HG002_H1"][idxh1].tolist() + all_rates_df[f"HG002_H2"][idxh2].tolist())
            vals_other = np.array(
                all_rates_df[f"{sample}_H1"][idxh1].tolist() + all_rates_df[f"{sample}_H2"][idxh2].tolist()
            )
            
            print(sample, scipy.stats.spearmanr(vals_child, vals_other))
            
            pa.figure(figsize=(5, 5))
            plt.title(sample)
            plt.ylabel("HG002")
            plt.xlabel(sample)
            plot_2d_density(vals_child, vals_other, cmap="jet")
            pa.savefig()
        
        """ More strictly, dropping all NA rows """
        all_rates_df_nonaned = all_rates_df.dropna()
        for sample in ["faux_HG002", "HG003", "HG004"]:
            vals_child = np.array(
                all_rates_df_nonaned[f"HG002_H1"].tolist() + all_rates_df_nonaned[f"HG002_H2"].tolist()
            )
            vals_other = np.array(
                all_rates_df_nonaned[f"{sample}_H1"].tolist() + all_rates_df_nonaned[f"{sample}_H2"].tolist()
            )
            
            print(sample, scipy.stats.spearmanr(vals_child, vals_other))
            
            pa.figure(figsize=(5, 5))
            plt.title(sample)
            plt.ylabel("HG002")
            plt.xlabel(sample)
            plot_2d_density(vals_child, vals_other, cmap="jet")
            pa.savefig()


    with pa.open_multipage_pdf("faux_HG002_plot"):
        plot_faux_hg002(all_rates_df_nonaned)


