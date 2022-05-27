import numpy as np
import pandas as pd
from meth5 import MetH5File


class Results:
    def __init__(self, reference_cpgs, m5_path, include_chroms=None):
        self.reference_cpgs = reference_cpgs
        self.segments = {}
        self.m5_path = m5_path
        self.include_chroms = include_chroms
    
    def annotate_cpgs(self, segments):
        segments["CpGs"] = segments.apply(
            lambda row: set(self.reference_cpgs.get_CGs(row["chrom"], row["start"], row["end"])), axis=1
        )
        return segments
    
    def compute_diffmet(self, mf, chrom, start, end):
        agg = mf[chrom].get_values_in_range(start, end).get_llr_site_readgroup_rate("haplotype")
        try:
            diffmet = np.abs(np.nanmean(agg[0][0]) - np.nanmean(agg[1][0]))
            return diffmet
        except:
            print("Can't compute difference for ", chrom, start, end)
    
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
        
        if self.include_chroms is not None:
            print("Filtering chromosomes")
            data = data.loc[data["chrom"].map(lambda x: x in self.include_chroms)].copy()
        
        print(f"Annotating CpGs for {key}")
        data = self.annotate_cpgs(data)
        print(f"Loading diffmet for {key}")
        with MetH5File(self.m5_path, "r") as mf:
            data["diffmet"] = data.apply(
                lambda row: self.compute_diffmet(mf, row["chrom"], row["start"], row["end"]), axis=1
            )
        print(f"Finished loading {key}")
        self.segments[key] = data
