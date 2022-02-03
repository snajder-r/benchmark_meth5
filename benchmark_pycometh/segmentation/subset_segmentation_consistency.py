from multiprocessing.pool import Pool

import tqdm
from nanoepiseg.segment import segment
from meth5.meth5 import MetH5File
from nanoepitools.annotations.annotations import GFFAnnotationsReader
import pandas as pd
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from nanoepitools.plotting.general_plotting import PlotArchiver
from benchmark_pycometh.config import module_config
from benchmark_pycometh.plotter import Plotter, CoordinateTranslator
from benchmark_pycometh.bsseq.create_mock_bs_from_nanopore import MetH5ToBedGraph

matplotlib.use("Agg")


def segmentation_worker(i, submatrix, max_segments):
    return (i, submatrix, segment(submatrix, max_segments))


class SegmentationConsistencyTest:
    def __init__(self, chrom, start, end):
        """INIT PLOTTING"""
        self.gff = GFFAnnotationsReader()
        self.gff.read(module_config.gff_file, only_protein_coding=False)
        self.gff.build_index()
        
        self.pa = PlotArchiver(
            "segmentation", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"}
        )
        self.pl = Plotter(self.gff, self.pa)
        
        """LOAD DATA"""
        m5_path = "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/HG003_cpg.h5"
        self.f = MetH5File(m5_path, "r", chunk_size=50000)
        self.matrix = (
            self.f[chrom].get_values_in_range(start, end).to_sparse_methylation_matrix(read_groups_key="haplotype")
        )
        self.coordinate_translator = CoordinateTranslator(self.matrix)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.converter = MetH5ToBedGraph(m5_path)
    
    def get_changepoints(self, matrix, segments):
        return matrix.genomic_coord[np.where(np.diff(segments))[0]]
    
    def plot_segmentation(self, seg, segtype, hold=False, **kwargs):
        figure_kwargs = {"figsize": (20, 4)}
        vlines = self.get_changepoints(self.matrix, seg)
        self.pl.plot_region(
            self.chrom,
            self.matrix.genomic_coord[0],
            self.matrix.genomic_coord_end[-1],
            figure_kwargs=figure_kwargs,
            ws=0,
            title=f"{segtype}",
            aggregate_samples=True,
            with_no_hp=False,
            coordinate_space=False,
            marker_height=0.9,
            fontsize=8,
            show_has_value_indicator=False,
            hold=hold,
            vlines_in_coord_space=vlines,
            vlinewidth=3,
            **kwargs,
        )
    
    def compute_changepoint_frequency_under_subsampling_input(self, rounds):
        for i in range(rounds):
            mask = np.random.random(size=self.matrix.read_names.shape) > 0.5
            submatrix = self.matrix.get_submatrix_from_read_mask(mask)
            yield i, submatrix, 30
    
    def compute_changepoint_frequency_under_subsampling(self, rounds=100):
        changepoint_frequency = np.zeros(self.matrix.genomic_coord.shape)
        with tqdm.tqdm(total=rounds) as pbar, Pool(16) as pool:
            input = self.compute_changepoint_frequency_under_subsampling_input(rounds)
            for i, submatrix, segments in pool.starmap(segmentation_worker, input):
                cps = self.get_changepoints(submatrix, segments)
                for cp in cps:
                    changepoint_frequency[self.matrix.genomic_coord == cp] += 1
                
                """ For methylkit to compute in R script methylkit_mockbsseq_consistency.R """
                output_file = module_config.mock_bsseq_template_file_consistency_testing.format(
                    sample="HG003", iteration=i
                )
                df = self.converter.convert_from_sparse_matrix(self.chrom, submatrix)
                df.to_csv(output_file, header=False, sep="\t", index=False)
                pbar.update(1)
        return changepoint_frequency
    
    def compute_changepoint_order_under_maxsegments(self, max_max_segments=30):
        changepoint_frequency = np.zeros((self.matrix.genomic_coord.shape[0], max_max_segments))
        input = ((i, self.matrix, max_segments) for i, max_segments in enumerate(range(2, max_max_segments + 1)))
        with tqdm.tqdm(total=max_max_segments - 1) as pbar, Pool(16) as pool:
            for i, _, segments in pool.starmap(segmentation_worker, input):
                cps = self.get_changepoints(self.matrix, segments)
                for cp in cps:
                    changepoint_frequency[self.matrix.genomic_coord == cp, i+1] += 1
                pbar.update(1)
        return changepoint_frequency
    
    def plot_changepoint_order(self, changepoint_order, max_max_segments=30, offset=5):
        coords = list(range(len(changepoint_order)))
        ylim_old = plt.ylim()
        
        def plot_coord(seg_num):
            return ylim_old[0] + seg_num - max_max_segments - offset
        
        for coord, min_needed in zip(coords, changepoint_order):
            if min_needed < max_max_segments - 1:
                print(min_needed, "at", plot_coord(min_needed))
                plt.vlines([coord + 0.5], plot_coord(min_needed), ylim_old[1], colors=["k"], linewidth=0.5)
        
        plt.hlines([plot_coord(x) for x in range(2, max_max_segments + 1)], 0, max(coords), linewidth=0.1, alpha=0.5)
        for i in range(2, max_max_segments + 1, 2):
            print(i, "is", plot_coord(i), "or", plot_coord(i)-0.4)
            plt.text(-6, plot_coord(i) - 0.4, str(int(i)), fontsize=8)
    
    def plot_changepoint_frequency(self, coords, cp_freq, text, offset=2, scale=2):
        plt.plot(coords, cp_freq * scale - offset)
        plt.hlines([-offset, -offset + scale], 0, max(coords), linewidth=0.1, alpha=0.5)
        plt.text(0, -offset - 0.8, text, fontsize=6)
        plt.text(-3, -offset - 0.1, "0", fontsize=6)
        plt.text(-3, -offset + 0.9 + (scale - 1), "1", fontsize=6)
    
    def plot(self, changepoint_frequency, changepoint_order):
        with self.pa.open_multipage_pdf("segmentation_consistency"):
            self.plot_segmentation([], "Nanoepiseg segmentation based on HG003", hold=True)
            if changepoint_frequency is not None:
                self.plot_changepoint_frequency(
                    list(range(len(changepoint_frequency))),
                    changepoint_frequency,
                    "Changepoint frequency under subsampling",
                    offset=6,
                    scale=2,
                )
            if changepoint_order is not None:
                self.plot_changepoint_order(changepoint_order, offset=10)
            
            methylkit_frequency = pd.read_csv(
                module_config.mock_bsseq_template_file_consistency_methylkit_result.format(sample="HG003"), sep="\t"
            )
            coord = [self.coordinate_translator(p) for p in methylkit_frequency["pos"]]
            self.plot_changepoint_frequency(
                coord,
                methylkit_frequency["changepoint_frequency"] / 100,
                "Methylkit changepoint frequency under subsampling",
                offset=10,
                scale=2,
            )
            plt.ylim(plt.ylim()[0] - max(changepoint_order) - 10, plt.ylim()[1])
            self.pa.savefig()


if __name__ == "__main__":
    test = SegmentationConsistencyTest("10", 8051635-5000, 8095709+5000)
    changepoint_frequency = test.compute_changepoint_frequency_under_subsampling()
    changepoint_order = test.compute_changepoint_order_under_maxsegments()
    test.__class__ = SegmentationConsistencyTest
    test.plot(changepoint_frequency / 100, changepoint_order)

