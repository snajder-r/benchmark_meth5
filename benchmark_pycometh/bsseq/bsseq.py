import pandas as pd

def load_metrates(filename, has_percent, has_chr, sort=True, recompute_frac = False):
    if recompute_frac:
        ret = pd.read_csv(filename, sep="\t", usecols=[0, 1, 2, 4, 5], names=["chrom", "start", "end", "pos", "neg"],
            dtype={"chrom": str})
        ret["fracmet"] = ret["pos"] / (ret["pos"] + ret["neg"])
        ret = ret.drop(["pos", "neg"], axis=1)
    else:
        ret = pd.read_csv(
            filename, sep="\t", usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "fracmet"], dtype={"chrom": str}
        )
    if has_percent:
        ret["fracmet"] = ret["fracmet"] / 100
    if has_chr:
        ret["chrom"] = ret["chrom"].map(lambda x: x.replace("chr", ""))
    if sort:
        ret = ret.sort_values(["chrom", "start", "end"])
    return ret
