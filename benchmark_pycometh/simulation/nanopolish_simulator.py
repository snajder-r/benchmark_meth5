import secrets
import random
from typing import List
from pathlib import Path
from multiprocessing import Process, Queue, current_process
from itertools import cycle

import tqdm
import pandas as pd
import numpy as np
from meth5 import MetH5File

from nanoepitools.math import p_to_llr


class OmicsSimlaLoader:
    def __init__(self, simdir: Path, profile_path: Path, profile_map_path: Path):
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
            self.rates[sample] = np.array([call[0] / call[1] for call in self.calls[sample]])
    
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
                calls = [tuple(int(a) for a in col.split(",")) for col in line[3:]]
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
                    first_col_to_read = 2 + max(0, first_pos - profile_read_cs) * 2
                    last_col_to_read = 2 + min(num_cs, last_pos - profile_read_cs + 1) * 2
                    self.profile_rates += [float(line[i]) for i in range(first_col_to_read, last_col_to_read, 2)]
                    self.segments += [segment_id] * num_cs
                    self.segment_types += [segment_type] * num_cs
                
                segment_id += 1
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
        return chrstart + first, chrstart + last
    
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
        """
        Simulate nanopore calls from a OMICSSimla simulation
        The default parameters for quality and read length distribution are esitmated from GIAB data
        """
        
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
        self.sample_hex_ids = {sample: secrets.token_hex(2) for sample in self.omics_simla.calls}
    
    def simulate_read_name(self, sample_id="", process_id="", chunk_id=""):
        read_name = secrets.token_hex(21)
        if sample_id == "":
            sample_id = read_name[:4]
        if process_id == "":
            process_id = read_name[4:8]
        if chunk_id == "":
            chunk_id = read_name[8:12]
        read_name = f"{sample_id}{process_id}-{chunk_id}-{read_name[12:20]}-{read_name[20:28]}-{read_name[28:42]}"
        return read_name
    
    def simulate_nanopore_read(self, sample, chrom, start, length, **kwargs):
        pos, rates = self.omics_simla.get_region_rates(sample, chrom, start, start + length)
        binary_read = np.random.rand(len(rates)) < rates
        binary_read = (binary_read - 0.5) * 2  # 1 methylated and -1 unmethylated
        absolute_llrs = p_to_llr(np.random.beta(self.quality_alpha, self.quality_beta, len(rates)))
        llrs = binary_read * absolute_llrs
        read_name = self.simulate_read_name(**kwargs)
        return read_name, chrom, pos, llrs
    
    def get_random_read_position_and_length(self):
        chrom = random.choice(list(self.chroms.keys()))
        pos = random.randint(*self.chroms[chrom])
        length_model = np.random.choice(list(range(len(self.readlen_model_weights))), p=self.readlen_model_weights)
        length = int(10 ** np.random.normal(self.readlen_log10_means[length_model], self.readlen_log10_vars[length_model]))
        length = min(self.chroms[chrom][1] - 1 - pos, length)
        return chrom, pos, length
    
    def simulate_nanopore_reads(self, sample, number, **kwargs):
        for _ in range(number):
            yield self.simulate_nanopore_read(sample, *self.get_random_read_position_and_length(), **kwargs)

    def group_calls(self, pos, llrs, min_dist=10):
        if len(pos)==0:
            return
        
        it = zip(pos, llrs)
        start, llrs_sum = next(it)
        end = start
        llrs_count = 1
        for p, llr in it:
            if p - end > min_dist:
                yield start, end, llrs_sum / llrs_count
                start = p
                end = p
                llrs_sum = llr
                llrs_count = 1
            else:
                end = p
                llrs_sum += llr
                llrs_count += 1
        if llrs_count > 0:
            yield start, end, llrs_sum / llrs_count
            
    def generate_meth5(self, m5_path: Path, readgroup_name, number_per_sample, chunksize=10000, n_procs=8):
        inq = Queue()
        outq = Queue()
        
        def work():
            chunk_i = 0
            process_id = "{0:#0{1}x}".format(current_process().ident,6)[2:]
            while True:
                n = inq.get()
                if n == -1:
                    outq.put(None)
                    return
                chunk_id = "{0:#0{1}x}".format(chunk_i, 6)[2:]
                chunk_i += 1
                rows = {}
                for sample in self.omics_simla.calls:
                    rows[sample] = []
                    readname_kwargs = dict(process_id=process_id, chunk_id=chunk_id, sample_id=self.sample_hex_ids[sample])
                    for read_name, chrom, pos, llrs in self.simulate_nanopore_reads(sample, chunksize, **readname_kwargs):
                        grouped_calls = self.group_calls(pos, llrs)
                        
                        to_add = [
                            {
                                "chromosome": chrom,
                                "start": start,
                                "end": end,
                                "read_name": read_name,
                                "log_lik_ratio": llr,
                            }
                            for start, end, llr in grouped_calls
                        ]
                        rows[sample] += to_add
                chunk = {s:pd.DataFrame(rows[s]) for s in rows}
                outq.put(chunk)
        
        ps = [Process(target=work) for _ in range(n_procs)]
        for p in ps:
            p.start()
        
        with MetH5File(m5_path, "w", compression="lzf") as mf:
            number_generated = 0
            read_group_map = {}
            read_group_id_map = {i: s for i, s in enumerate(self.omics_simla.calls, 1)}
            
            while number_generated < number_per_sample:
                chunksize = min(number_per_sample - number_generated, chunksize)
                number_generated += chunksize
                inq.put(chunksize)
            
            for p in ps:
                inq.put(-1)
            
            ps_working = len(ps)
            with tqdm.tqdm(total=number_per_sample * len(read_group_id_map)) as pbar:
                while ps_working > 0:
                    chunk = outq.get()
                    if chunk is None:
                        ps_working -= 1
                        continue
                    
                    for sample_id, sample in read_group_id_map.items():
                        nreads = 0
                        for read in set(chunk[sample]["read_name"]):
                            read_group_map[read] = sample_id
                            nreads +=1
                        mf.add_to_h5_file(chunk[sample], postpone_sorting_until_close=True)
                        pbar.update(nreads)
            
            for p in ps:
                p.join()
            
            mf.annotate_read_groups(readgroup_name, read_group_map, labels=read_group_id_map)
            mf.create_chunk_index()
