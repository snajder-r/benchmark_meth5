import pandas as pd
import numpy as np
from benchmark_pycometh.config import module_config

for sample, files in module_config.bs_seq_files.items():
    if sample != "HG004":
        continue
    df_r1 = pd.read_csv(files["R1"], sep="\t", names=["chrom", "start", "end", "frac", "pos", "neg"], header=None)
    df_r2 = pd.read_csv(files["R2"], sep="\t", names=["chrom", "start", "end", "frac", "pos", "neg"], header=None)
    df_r1 = df_r1.set_index(["chrom", "start", "end"])
    df_r2 = df_r2.set_index(["chrom", "start", "end"])
    
    df = df_r1.merge(df_r2, how="outer", left_index=True, right_index=True, suffixes=["_r1", "_r2"])
    
    for a in "neg", "pos":
        for r in "r1", "r2":
            df[f"{a}_{r}"] = df[f"{a}_{r}"].map(lambda x: 0 if np.isnan(x) else x)
        df[a] = df[f"{a}_r1"] + df[f"{a}_r2"]
    
    df["cov"] = (df["pos"] + df["neg"])
    df["frac"] = df["pos"] / df["cov"]
    df = df.loc[df["cov"] > 0]
    
    df = df.drop(["pos_r1", "neg_r1", "frac_r1", "pos_r2", "neg_r2", "frac_r2"], axis=1)
    
    df["strand"] = "+"
    df = df.reset_index()[["chrom", "start", "end", "frac", "cov", "pos", "strand"]]
    df.to_csv(f"{sample}_rep1_rep2_merged.bedGraph.gz", sep="\t", index=False, header=False)
