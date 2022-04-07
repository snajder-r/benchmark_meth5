import sys
import argparse

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
    def __init__(self, meth5_file, readgroup=None, readgroup_key=None, fasta=module_config.reference):
        self.ref = pyfaidx.Fasta(fasta, "r")
        self.f = MetH5File(meth5_file, "r")
        self.readgroup_key = readgroup_key
        if self.readgroup_key is not None:
            group_labels = {v:k for k,v in self.f.get_all_read_groups(self.readgroup_key).items()}
            self.readgroup = group_labels[readgroup] # translate string label back to rg id
        
    
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
            
            if self.readgroup_key is not None:
                agg = self.f[chrom].get_all_values().get_llr_site_readgroup_aggregate(lambda x: agg_fun(x, thres=2.0), group_key=self.readgroup_key)
                agg, r = agg[self.readgroup]
            else:
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
    parser = argparse.ArgumentParser(description='Convert m5 to pseudo-bisulfite.')
    parser.add_argument('m5file', metavar="m5file", type=str, help='Input M5 file')
    parser.add_argument('outfile', metavar="outfile", type=str, help='Output tsv file')
    parser.add_argument('fasta', metavar="fasta", type=str, help='reference genome file')
    parser.add_argument('--readgroup', metavar="readgroup", required=False, type=str, help='read group in m5 file', default=None)
    parser.add_argument('--readgroup_key', metavar="readgroup_key", required=False, type=str, help='category of read group', default=None)
    args = parser.parse_args()
    
    converter = MetH5ToBedGraph(args.m5file, readgroup=args.readgroup, readgroup_key=args.readgroup_key, fasta=args.fasta)
    converter.convert_all(args.outfile)
