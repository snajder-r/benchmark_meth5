from pathlib import Path
import pandas as pd
import numpy as np
from benchmark_pycometh.config import module_config

basedir = Path("/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged")
files = {s: basedir.joinpath(f"{s}_mock_bsseq.bedGraph") for s in ["HG003", "HG004"]}

df = {
    s: pd.read_csv(
        files[s],
        sep="\t",
        usecols=list(range(6)),
        names=["chrom", "start", "end", "frac", "cov", "pos"],
        dtype={"chrom": str, "start": int, "end": int, "frac": float, "cov": int, "pos": int},
        header=None,
        index_col=None,
    )
    for s in files
}

df = {s: df[s].set_index(["chrom", "start", "end"]) for s in files}

df = df["HG003"].merge(df["HG004"], how="outer", left_index=True, right_index=True, suffixes=["_r1", "_r2"])

for a in "cov", "pos":
    for r in "r1", "r2":
        df[f"{a}_{r}"] = df[f"{a}_{r}"].map(lambda x: 0 if np.isnan(x) else x)
    df[a] = df[f"{a}_r1"] + df[f"{a}_r2"]

df["frac"] = df["pos"] / df["cov"]
df = df.loc[df["cov"] > 0]

df = df.drop(["pos_r1", "cov_r1", "frac_r1", "pos_r2", "cov_r2", "frac_r2"], axis=1)

df["strand"] = "+"
df = df.reset_index()[["chrom", "start", "end", "frac", "cov", "pos", "strand"]]
df["cov"] = df["cov"].astype(int)
df["pos"] = df["pos"].astype(int)
df.to_csv(basedir.joinpath("parents_mock_bsseq.bedGraph"), sep="\t", index=False, header=False)
