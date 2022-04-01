import secrets
import random
from pathlib import Path

import pandas as pd
import numpy as np

from nanoepitools.math import p_to_llr

class OmicsSimlaLoader:
    def __init__(self, simdir: Path, profile_path: Path, profile_map_path:Path):
        self.profile_path = profile_path
        self.profile_map_path = profile_map_path
        self.omics_simla_summary_file = simdir.joinpath("Methylation_summary.txt")
        self.omics_simla_file = simdir.joinpath("sim1.methy")
        self.calls = {}
        self.rates = {}
        self.segments = []
        self.segment_types = []
        self.profile_rates = []
        self.chrom = None
        self.pos = None
        self.summary = None
        self.load()

    def compute_rates(self):
        for sample in self.calls:
            self.rates[sample] = np.array([call[0]/call[1] for call in self.calls[sample]])

    def index(self):
        index = self.summary.index.map(lambda x: x.split(":"))
        self.chrom = np.array([i[0] for i in index])
        self.pos = np.array([int(i[1]) for i in index])
        self.index = {}
        for chrom in set(self.chrom):
            idx = np.where(self.chrom == chrom)[0]
            self.index[chrom] = (idx[0], idx[-1])

    def load(self):
        """ First load the summary and the map from cpg to coordinate """
        self.summary = pd.read_csv(self.omics_simla_summary_file, sep=" ").set_index("Pos")
        first_pos_str = self.summary.index[0].split(":")
        last_pos_str = self.summary.index[-1].split(":")
        for i, line in enumerate(open(self.profile_map_path, "r")):
            line = line.split(" ")
            if line[0] != first_pos_str[0] and line[0] != last_pos_str[0]:
                continue
            if line[0] == first_pos_str[0] and line[1] == first_pos_str[1]:
                first_pos = i
            if line[0] == last_pos_str[0] and line[1] == last_pos_str[1]:
                last_pos = i
                break

        """ Now load the methylation calls  """
        with open(self.omics_simla_file, "r") as f:
            for line in f.readlines():
                line = line.strip().split(" ")
                sample = line[0]
                calls = [tuple(int(a) for a in col.split(",")) for col in line[3:] ]
                self.calls[sample] = calls

        """ Finally load the mapping from CpG to segment """
        profile_read_cs = 0
        with open(self.profile_path, "r") as f:
            segment_id = 0
            for line in f.readlines():
                line = line.strip().split(" ")
                segment_type = int(line[0])
                num_cs = int(line[1])
                if profile_read_cs + num_cs < first_pos:
                    profile_read_cs += num_cs
                    continue
                else:
                    first_col_to_read = 2+max(0, first_pos - profile_read_cs)*2
                    last_col_to_read = 2+min(num_cs, last_pos - profile_read_cs+1)*2
                    self.profile_rates += [float(line[i]) for i in range(first_col_to_read, last_col_to_read, 2)]
                    self.segments += [segment_id]*(last_col_to_read-first_col_to_read)
                    self.segment_types += [segment_type]*((last_col_to_read-first_col_to_read)//2)

                segment_id+=1
                profile_read_cs += num_cs
                if profile_read_cs > last_pos:
                    break
        self.profile_rates = np.array(self.profile_rates)
        self.segment_types = np.array(self.segment_types)
        self.compute_rates()
        self.index()

    def find_index(self, chrom, start, end):
        chrstart, chrend = self.index[chrom]
        first = np.searchsorted(self.pos[chrstart:chrend], start)
        last = np.searchsorted(self.pos[chrstart:chrend], end)
        return chrstart+first, chrstart+last

    def get_region_rates(self, sample, *args):
        first, last = self.find_index(*args)
        return self.pos[first:last], self.rates[sample][first:last]


class Simulator:
    def __init__(
        self,
        *args,
        quality_alpha=0.640308238,
        quality_beta=0.208756,
        readlen_log10_means=[2.88, 3.8, 4.75],
        readlen_log10_vars=[0.38, 0.48, 0.29],
        readlen_model_weights=[0.24, 0.63, 0.13],
    ):
        # default parameters previously estimated from GIAB data
        self.omics_simla = OmicsSimlaLoader(*args)
        self.quality_alpha = quality_alpha
        self.quality_beta = quality_beta
        self.readlen_log10_means = readlen_log10_means
        self.readlen_log10_vars = readlen_log10_vars
        self.readlen_model_weights = readlen_model_weights
        self.chroms = {}
        unique_chroms = set(self.omics_simla.chrom)
        for chrom in unique_chroms:
            chrom_idx = self.omics_simla.chrom == chrom
            self.chroms[chrom] = (min(self.omics_simla.pos[chrom_idx]), max(self.omics_simla.pos[chrom_idx]))

    def simulate_read_name(self):
        read_name = secrets.token_hex(16)
        read_name = f"{read_name[:8]}-{read_name[8:12]}-{read_name[12:16]}-{read_name[16:20]}-{read_name[20:32]}"
        return read_name

    def simulate_nanopore_read(self, sample, chrom, start, length):
        pos, rates = self.omics_simla.get_region_rates(sample, chrom, start, start + length)
        binary_read = np.random.rand(len(rates)) < rates
        binary_read = (binary_read - 0.5) * 2  # 1 methylated and -1 unmethylated
        absolute_llrs = p_to_llr(np.random.beta(self.quality_alpha, self.quality_beta, len(rates)))
        llrs = binary_read * absolute_llrs
        read_name = self.simulate_read_name()
        return read_name, pos, llrs

    def get_random_read_position_and_length(self):
        chrom = random.choice(list(self.chroms.keys()))
        pos = random.randint(*self.chroms[chrom])
        length_model = np.random.choice(list(range(len(self.readlen_model_weights))), p=self.readlen_model_weights)
        length = 10**np.random.normal(self.readlen_log10_means[length_model], self.readlen_log10_vars[length_model])
        length = min(pos+length, self.chroms[chrom]-1)
        return chrom, pos, length

    def simulate_nanopore_reads(self, sample, number):
        for _ in range(number):
            yield self.simulate_nanopore_read(
                sample, *self.get_random_read_position_and_length()
            )

    def generate_meth5(self, m5_path: Path, ground_truth_path: Path, readgroup_name):
        pass