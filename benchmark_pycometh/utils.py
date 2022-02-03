import numpy as np
from typing import List, Dict
from meth5.meth5 import MetH5File

from nanoepitools.reference_cpgs import ReferenceCpGs

def merge_duplicate_diffmet_hits(oldhits, before=500, after=500):
    # Merge so we don't count them double
    newhits = []
    a = 0
    b = 0
    for i, hiti in enumerate(oldhits):
        a += 1
        duplicate = -1
        for j, hitj in enumerate(newhits):
            if hiti["start"] - before < hitj["end"] and hitj["start"] < hiti["end"] + after and hiti["chrom"] == hitj["chrom"]:
                duplicate = j
                break
        if duplicate >= 0:
            newhits[duplicate] = {
                "chrom": hiti["chrom"],
                "start": min(hiti["start"], newhits[duplicate]["start"]),
                "end": max(hiti["end"], newhits[duplicate]["end"]),
            }
        else:
            b += 1
            newhits.append(hiti)
    return newhits

def unions(iterables):
    return set(x for xl in iterables for x in xl)

def count_cpgs_available(chrom, start, end, m5s: List[MetH5File], ref_cpgs: ReferenceCpGs, llr_threshold=2.0):
    covered_coords_intersect = None
    for m5 in m5s:
        coverage, ranges = m5[chrom].get_values_in_range(start, end).get_llr_site_aggregate(lambda x: (np.abs(x) > llr_threshold).sum())
        covered_coords = unions(ref_cpgs.get_CGs(chrom, r[0]-2, r[1]+2, formatted=True) for c,r in zip(coverage, ranges) if c > 1)
        if covered_coords_intersect is None:
            covered_coords_intersect = covered_coords
        else:
            covered_coords_intersect = covered_coords_intersect.intersection(covered_coords)
    return covered_coords
    