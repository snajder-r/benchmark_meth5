import tqdm
import numpy as np
import pandas as pd
from meth5.meth5 import MetH5File
from benchmark_pycometh.config import module_config
import pyfaidx


def agg_fun(llrs, thres):
    num = (np.abs(llrs) > thres).sum()
    pos = (llrs > thres).sum()
    if num > 0:
        frac = pos / num
    else:
        frac = 0.5
    return {"cov": num, "pos": pos, "frac": frac}


def disentangle_combined_calls(ranges, aggs, seq):
    for agg, r in zip(aggs, ranges):
        if agg["cov"] == 0:
            continue
        
        if r[1] == r[0]:
            yield r[0], agg
        else:
            subseq = seq[r[0] : r[1]]
            for offset in range(len(subseq) - 2):
                if subseq[offset : offset + 2] == "CG":
                    yield r[0] + offset, agg


class MetH5ToBedGraph:
    def __init__(self, meth5_file):
        self.ref = pyfaidx.Fasta(module_config.reference, "r")
        
        self.f = MetH5File(meth5_file, "r")
    
    def create_dataframe(self, chrom, ranges, frequencies, sequence):
        new_bs = []
        for in_pos, in_row in disentangle_combined_calls(ranges, frequencies, sequence):
            new_bs.append({**in_row, "start": in_pos, "end": in_pos + 2, "chrom": chrom})
        
        new_bs = pd.DataFrame(new_bs)
        new_bs["strand"] = "+"
        new_bs = new_bs[["chrom", "start", "end", "frac", "cov", "pos", "strand"]]
        return new_bs
    
    def convert_all(self, out_tsv):
        virgin = True
        for chrom in tqdm.tqdm(self.f.get_chromosomes()):
            if chrom not in self.ref:
                continue
            agg, r = self.f[chrom].get_all_values().get_llr_site_aggregate(lambda x: agg_fun(x, thres=2.0))
            new_bs = self.create_dataframe(chrom, r, agg, self.ref[chrom])
            new_bs.to_csv(
                out_tsv,
                mode="a" if not virgin else "w",
                header=False,
                sep="\t",
                index=False,
            )
            virgin = False

    def convert_from_sparse_matrix(self, chrom, sparse_matrix):
        agg = [agg_fun(sparse_matrix.met_matrix[:, i], 2.0) for i in range(sparse_matrix.met_matrix.shape[1])]
        r = list(zip(sparse_matrix.genomic_coord, sparse_matrix.genomic_coord_end))
        new_bs = self.create_dataframe(chrom, r, agg, self.ref[chrom])
        return new_bs

if __name__ == "__main__":
    sample = "HG004"
    converter = MetH5ToBedGraph(module_config.meth5_template_file.format(sample=sample))
    converter.convert_all(module_config.mock_bsseq_template_file.format(sample=sample))

    
    virgin = {"H1": True, "H2": True}
    hp_labels = {int(k): v for k, v in f.h5_fp["reads/read_groups/haplotype"].attrs.items()}
    for chrom in tqdm.tqdm(f.get_chromosomes()):
        if chrom not in ref:
            continue
        agg_dict = f[chrom].get_all_values().get_llr_site_readgroup_aggregate(lambda x: agg_fun(x, thres=2.0), "haplotype")
        for hp_id in agg_dict:
            new_bs = []
            hp = hp_labels.get(hp_id, "none")
            if hp not in {"H1", "H2"}:
                continue
            agg, r = agg_dict[hp_id]
            for in_pos, in_row in disentangle_combined_calls(r, agg, ref[chrom]):
                new_bs.append({**in_row, "start": in_pos, "end": in_pos + 2, "chrom": chrom})
            new_bs = pd.DataFrame(new_bs)
            new_bs["strand"] = "+"
            new_bs = new_bs[["chrom", "start", "end", "frac", "cov", "pos", "strand"]]
            new_bs.to_csv(
                module_config.mock_bsseq_template_file.format(sample=f"{sample}_{hp}"),
                mode="a" if not virgin[hp] else "w",
                header=False,
                sep="\t",
                index=False,
            )
            virgin[hp] = False
