import random
from collections import namedtuple

import tqdm
import numpy as np
import pandas as pd


class SegmentsComparer:
    def __init__(self, gt_segments, predicted_segments):
        self.gt_segments = gt_segments
        self.predicted_segments = predicted_segments
        

    def compute_dist_from_a_to_b(self, a, b):
        total = len(a)
        for seg in tqdm.tqdm(a.itertuples(), total=total, disable=True):
            d = (b["start"] - seg.start).abs()
            nearest_dist = d.min()
            d = (b["end"] - seg.start).abs()
            nearest_dist = min(nearest_dist, d.min())
            d = (b["start"] - seg.end).abs()
            nearest_dist = min(nearest_dist, d.min())
            d = (b["end"] - seg.end).abs()
            nearest_dist = min(nearest_dist, d.min())
            yield nearest_dist

    def compute_gt_dist_to_predicted(self, segments):
        return self.compute_dist_from_a_to_b(self.gt_segments, segments)

    def compute_gt_dist_to_all_predicted(self):
        dists = {}
        for tool in self.predicted_segments:
            dists[tool] = np.array(list(self.compute_gt_dist_to_predicted(self.predicted_segments[tool])))
        return dists

    def compute_predicted_dist_to_gt(self, segments):
        return self.compute_dist_from_a_to_b(segments, self.gt_segments)

    def compute_predicted_dist_to_all_gt(self):
        dists = {}
        for tool in self.predicted_segments:
            dists[tool] = np.array(list(self.compute_predicted_dist_to_gt(self.predicted_segments[tool])))
        return dists

    def get_segments_around(self, pos, gt = None):
        if gt is None:
            gt = self.gt_segments
        last_seg = None
        seg_it = gt.itertuples()
        try:
            while True:
                seg = next(seg_it)
                if seg.start <= pos <= seg.end:
                    try:
                        next_seg = next(seg_it)
                    except StopIteration:
                        next_seg = None
                    return last_seg, seg, next_seg
                elif last_seg is not None:
                    if last_seg.end <= pos <= seg.start:
                        return last_seg, None, seg
                if pos < seg.start:
                    return None, None, None
                last_seg = seg
        except StopIteration:
            return None, None, None

    def find_separated_segments(self, cp, gt=None, max_dist_rel=0.05):
        """If the given changepoint marks a ground truth changepoint, this function will return a tuple
        with the two segments it separates"""
        last_seg, in_seg, next_seg = self.get_segments_around(cp, gt=gt)
        if in_seg is None and last_seg is not None and next_seg is not None:
            # If the cp is not in a sgement, it is between segments, which is a perfect changepoint
            return last_seg, next_seg
        elif in_seg is not None:
            if in_seg.end - cp < cp - in_seg.start:
                dist = in_seg.end - cp
                neighbor_segs = in_seg, next_seg
            else:
                dist = cp - in_seg.start
                neighbor_segs = last_seg, in_seg

            max_dist = (in_seg.end - in_seg.start) * max_dist_rel
            if dist <= max_dist:
                return neighbor_segs
        return None, None

    def classify_changepoint(self, cp, max_dist_rel=0.05):
        seg_a, seg_b = self.find_separated_segments(cp)
        if seg_a is None or seg_b is None:
            return 0, "unsupported"
        if seg_a.theta != 0 or seg_b.theta != 0:
            return 2, "diff_seg_cp"
        else:
            return 1, "true_cp"

    def classify_all_segments(self, tool):
        segments = self.predicted_segments[tool]
        classes = {"unsupported": 0, "diff_seg_cp": 0, "true_cp": 0}
        with tqdm.tqdm(total=len(self.predicted_segments[tool])) as pbar:
            cp_class = (0, "unsupported")
            for seg in segments.itertuples():
                # We only count one (the better) class for a pair of segment end and next segment start
                cp_class = max(cp_class, self.classify_changepoint(seg.start))
                classes[cp_class[1]] += 1
                cp_class = self.classify_changepoint(seg.end)
                pbar.update(1)
            classes[cp_class[1]] += 1
        return classes

    def count_segments_identified(self, tool):
        gt = self.gt_segments.copy()
        gt["found_right"] = False
        gt["found_left"] = False
        cps_unused = []
        with tqdm.tqdm(total=self.predicted_segments[tool].shape[0]) as pbar:
            for segment in self.predicted_segments[tool].itertuples():
                for pos in segment.start, segment.end:
                    seg_a, seg_b = self.find_separated_segments(pos, gt=gt)
                    if seg_a is None or seg_b is None:
                        cps_unused.append(pos)
                        continue
                    if gt.loc[seg_a.Index,"found_left"]:
                        gt = gt.drop(seg_a.Index)
                    else:
                        gt.loc[seg_a.Index, "found_right"] = True
                    if gt.loc[seg_b.Index, "found_right"]:
                        gt = gt.drop(seg_b.Index)
                    else:
                        gt.loc[seg_b.Index,"found_left"] = True
                pbar.update(1)
        return gt, cps_unused


def permute_segments(original_segments):
    segment_lengths = original_segments.apply(lambda x: x["end"] - x["start"], axis=1).tolist()
    ends = np.array([0] + original_segments["end"].tolist())[:-1]
    starts = np.array(original_segments["start"])
    gaps = starts - ends
    offset = gaps[0]
    gaps = gaps[1:]
    random.shuffle(segment_lengths)
    random.shuffle(gaps)
    gaps = iter(gaps)
    region = namedtuple("region", ["start", "end"])

    ret = []
    for l in segment_lengths:
        ret.append(region(offset, offset + l))
        try:
            offset += next(gaps) + l
        except StopIteration:
            continue
    return pd.DataFrame({"chrom":"1", "start": [r.start for r in ret], "end": [r.end for r in ret]})

def convert_to_cg_index(segments):
    start = 0
    cpg_indices = []
    last_end_coord = 0
    for row in segments.itertuples():
        gap_cpgs = reference_cpgs.get_CGs(row.chrom, last_end_coord, row.start-1)
        cpgs = reference_cpgs.get_CGs(row.chrom, row.start, row.end)
        last_end_coord = row.end
        start = start + len(gap_cpgs)
        end = start + len(cpgs)
        cpg_indices.append((start, end))
        start = end
    segments["real_start"] = segments["start"]
    segments["real_end"] = segments["end"]
    segments["start"] = [i[0] for i in cpg_indices]
    segments["end"] = [i[1] for i in cpg_indices]

def load_pycometh(file, merge_window="no"):
    pm = pd.read_csv(file, sep="\t", usecols=[0,1,2,4], names=["chrom", "start", "end", "type"], dtype={"chrom":str})
    if merge_window == "no":
        return pm
    rows = []
    start = 0
    for row in pm.itertuples():
        if row.type == "window_start" and start != 0:
            rows.append({"chrom": row.chrom, "end": row.end, "start": start})
        elif row.type == "window_end":
            start = row.start
        else:
            rows.append({"chrom": row.chrom, "end": row.end, "start": row.start})
    return pd.DataFrame(rows)
