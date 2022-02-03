from pathlib import Path
import numpy as np
import pandas as pd
import tqdm
from meth5.meth5 import MetH5File
from benchmark_pycometh.config import module_config

dir = Path("/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation/HG003/")

chroms = list({x.name.split("_")[1] for x in dir.glob("chr_*.tsv")})

def cleanup_segmentation(segment_p: np.ndarray, segments: np.ndarray, min_parameter_diff=0.1) -> np.ndarray:
    new_segments = segments.copy()
    for segment in sorted(list(set(segments))):
        if len(set(new_segments)) <= 1:
            # No need to go on if it's all just one segment
            break
        if segment == new_segments[-1]:
            candidate_replace = new_segments[new_segments != segment][-1]
        else:
            candidate_replace = segment + 1
        absdif = np.abs(segment_p[:, segment] - segment_p[:, candidate_replace]).max()
        if absdif < min_parameter_diff:
            new_segments[new_segments == segment] = candidate_replace
    
    return np.array(new_segments)

meth5 = {s:module_config.meth5_template_file.format(sample=s) for s in ["HG003"]}
meth5 = {s:MetH5File(v, "r") for s,v in meth5.items()}
haplotype_map = {s:{int(k): v for k, v in meth5[s].h5_fp["reads/read_groups"]["haplotype"].attrs.items()} for s in meth5}

def get_bs(sample, hp, chrom, start, end):
    v = meth5[sample][chrom].get_values_in_range(start, end)
    llrs = v.get_llrs()
    mask = np.array([haplotype_map[sample].get(h, "none") == hp for h in v.get_read_groups(group_key="haplotype")])
    llrs = llrs[mask]
    nums = (np.abs(llrs)>2).sum()
    if nums == 0:
        return np.nan
    return (llrs > 2).sum() / nums

def get_parameters(chrom, segments):
    chunk_start = min(segments["start"])
    chunk_end = max(segments["end"])
    
    full_matrix = []
    for sample in "HG003",:
        v = meth5[sample][chrom].get_values_in_range(chunk_start, chunk_end)
        llrs = v.get_llrs()
        ranges = v.get_ranges()
        matrix = []
        for hp in "H1", "H2":
            hp_mask = np.array([haplotype_map[sample].get(h, "none") == hp for h in v.get_read_groups(group_key="haplotype")])
            llrs_hp = llrs[hp_mask]
            ranges_hp = ranges[hp_mask]
            col = []
            for idx, row in segments.iterrows():
                seg_mask = (ranges_hp[:,0] < row["end"]) & (ranges_hp[:,1] > row["start"])
                llrs_seg = llrs_hp[seg_mask]
                nums = (np.abs(llrs_seg)>2).sum()
                if nums == 0:
                    bs = 0.5
                else:
                    bs = (llrs_seg > 2).sum() / nums
                col.append(bs)
            matrix.append(col)
        matrix = np.array(matrix)
        full_matrix.append(matrix)
    return np.concatenate(tuple(full_matrix))


def compute_new_segments(file):
    df = pd.read_csv(file, usecols=[0, 1, 2], names=["chrom", "start", "end"], sep="\t", skiprows=1)
    df = df.sort_values("start")
    segment_p = get_parameters(chrom, df)
    segments = np.array(df.index)
    clean = cleanup_segmentation(segment_p, segments, 0.35)
    for seg in set(clean):
        start = df.iloc[np.where(clean == seg)[0][0]]["start"]
        end = df.iloc[np.where(clean == seg)[0][-1]]["end"]
        yield {"chrom": chrom, "start": start, "end": end}

new_seg = []
for chrom in chroms:
    for file in tqdm.tqdm(list(dir.glob(f"chr_{chrom}_chunk_*_segmented_cpg.tsv"))):
        for seg in compute_new_segments(file):
            new_seg.append(seg)

pd.DataFrame(new_seg).to_csv("/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation/HG003_35.bed", sep="\t", index=False, header=False)