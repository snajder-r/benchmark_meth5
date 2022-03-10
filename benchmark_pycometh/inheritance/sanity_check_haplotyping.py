import pysam
import pandas as pd
import numpy as np
from benchmark_pycometh.inheritance.map_child_haplotypes_to_parent import load_svs
from benchmark_pycometh.inheritance.investigate_methylation_inheritance import (
    get_segments,
    load_hp_parent_child_mapping,
)
from benchmark_pycometh.config import module_config
from typing import Optional
import pysam


def find_base_in_alignment(
    alignment: pysam.AlignedSegment, pos: int, bam_stores_revcomp: bool = False
) -> Optional[str]:
    idx_q = 0
    idx_r = pos - alignment.reference_start
    if bam_stores_revcomp:
        seq = alignment.query_sequence
    else:
        seq = alignment.get_forward_sequence()
    
    if seq is None:
        return None
    
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8}
        query_consumed = op in {0, 1, 4, 7, 8}
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r < 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]
                return base
            else:
                # position has been deleted
                return None


svs_hg002, svs_hg003, svs_hg004 = load_svs()
svs = dict(HG002=svs_hg002, HG003=svs_hg003, HG004=svs_hg004)
mapping_hps = load_hp_parent_child_mapping()


haplotags = {"HG002": "/home/r933r/data/projects/nanopore/pycometh_benchmark/whatshap/haplotags_HG002.tsv"}
for sample in haplotags:
    haplotags[sample] = pd.read_csv(haplotags[sample], sep="\t").set_index("read_name")["group"]

bams = {"HG002": "/home/r933r/data/projects/nanopore/pycometh_benchmark/mapping/HG002/100.sorted.filtered.bam"}

sample = "HG002"

with pysam.AlignmentFile(bams[sample]) as f:
    for read in f.fetch(until_eof=True):
        if read.query_sequence is None:
            continue
        hp = haplotags[sample].get(read.query_name, None)
        if hp is None or hp == "none":
            continue
        
        chrom_svs = svs[sample].loc[read.reference_name].reset_index()
        segment_svs = chrom_svs.loc[
            chrom_svs.apply(lambda x: read.reference_start < x.pos < read.reference_end, axis=1)
        ]
        
        hp_evidence = {"H1": 0, "H2": 0}
        for row in segment_svs.itertuples():
            base = find_base_in_alignment(read, row.pos, bam_stores_revcomp=True)
            if base == row.alt:
                hp_evidence[f"H{row.hp}"] += 1
            if base == row.ref:
                hp_evidence["H1" if row.hp == 2 else "H2"] += 1
        my_hp = (
            "H1"
            if hp_evidence["H1"] > hp_evidence["H2"]
            else "H2"
            if hp_evidence["H2"] > hp_evidence["H1"]
            else "unknown"
        )
        print("WhatsHap HP", hp, "my HP: ", my_hp)
