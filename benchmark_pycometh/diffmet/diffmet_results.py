import numpy as np
import pandas as pd
from meth5 import MetH5File


class Results:
    def __init__(self, reference_cpgs, m5_path, include_chroms=None, read_group_key="haplotype", samples=None):
        self.reference_cpgs = reference_cpgs
        self.segments = {}
        assert (read_group_key is None or samples is None)
        
        if read_group_key is not None:
            self.m5_paths = {"asm": m5_path}
        else:
            self.m5_paths = m5_path
            
        self.include_chroms = include_chroms
        self.read_group_key = read_group_key
        self.samples = samples
    
    def annotate_cpgs(self, segments):
        segments["CpGs"] = segments.apply(
            lambda row: set(self.reference_cpgs.get_CGs(row["chrom"], row["start"], row["end"])), axis=1
        )
        return segments
    
    def compute_diffmet(self, mfs, chrom, start, end, min_calls=2, return_num_calls=False):
        if self.read_group_key is not None:
            vals = {"asm":mfs["asm"][chrom].get_values_in_range(start, end)}
            agg = vals["asm"].get_llr_site_readgroup_rate(self.read_group_key)
        else:
            vals = {s:mfs[s][chrom].get_values_in_range(start, end) for s in self.samples}
            agg = [vals[s].get_llr_site_rate() for s in self.samples]
        
        if return_num_calls:
            n_calls = sum((np.abs(v.get_llrs())>2).sum() for v in vals.values())
        
        if len(agg) < 2:
            return np.nan
        if len(agg[0][1]) < min_calls or len(agg[1][1]) < min_calls:
            # Too few calls:
            if return_num_calls:
                return np.nan, n_calls
            return np.nan
        try:
            diffmet = np.abs(np.nanmean(agg[0][0]) - np.nanmean(agg[1][0]))
            if return_num_calls:
                return diffmet, n_calls
            return diffmet
        except:
            print("Can't compute difference for ", chrom, start, end)
            if return_num_calls:
                return np.nan, n_calls
            return np.nan
    
    def load_segments(self, key, file, caller):
        print(f"Loading {key}")
        if caller == "gt":
            data = pd.read_csv(
                file, sep="\t", names=["chrom", "start", "end", "segment_type", "theta"], dtype={"chrom": str}
            )
        elif caller == "pycometh":
            data = pd.read_csv(file, sep="\t", dtype={"chromosome": str}).rename({"chromosome": "chrom"}, axis=1)
            data = data.loc[data["adj_pvalue"] < 0.05].copy()
        elif caller == "methcp":
            data = pd.read_csv(
                file,
                sep="\t",
                names=["chrom", "start", "end", "nC.valid", "nC", "diffmet_mcp", "cov", "adj_pvalue"],
                dtype={"chrom": str},
            )
            data = data.loc[data["adj_pvalue"] < 0.05].copy()
        
        if self.include_chroms is not None:
            print("Filtering chromosomes")
            data = data.loc[data["chrom"].map(lambda x: x in self.include_chroms)].copy()
        
        print(f"Annotating CpGs for {key}")
        data = self.annotate_cpgs(data)
        print(f"Loading diffmet for {key} ({data.shape[0]} segments)")
        mfs = {s: MetH5File(self.m5_paths[s], "r") for s in self.m5_paths}
        diffmet_and_ncalls = data.apply(
            lambda row: self.compute_diffmet(mfs, row["chrom"], row["start"], row["end"], return_num_calls=True), axis=1
        )
        data["n_calls"] = diffmet_and_ncalls.map(lambda x: x[1])
        data["diffmet"] = diffmet_and_ncalls.map(lambda x: x[0])
        data = data.loc[data["diffmet"].map(lambda x: not np.isnan(x))]
        print(f"Finished loading {key} ({data.shape[0]} segments)")
        self.segments[key] = data
