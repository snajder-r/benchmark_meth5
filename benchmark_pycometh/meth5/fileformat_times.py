import time
import random
from modbampy import ModBam, pileup
import pyfaidx
import tqdm
from meth5.meth5 import MetH5File
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from nanoepitools.plotting.general_plotting import PlotArchiver
from benchmark_pycometh.config import module_config

matplotlib.use("Agg")
chroms = [str(i) for i in range(1, 22)]

pa = PlotArchiver("meth5", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

mf_lzf = MetH5File(module_config.meth5_template_file.format(sample="HG003"), "r", chunk_size=50000)
mf_gz = MetH5File(module_config.meth5_gzipped_template_file.format(sample="HG003"), "r", chunk_size=50000)


bf = module_config.modbam_template_file.format(sample="HG003")
cf = module_config.modcram_template_file.format(sample="HG003")
ff = pyfaidx.Fasta(module_config.reference)
chromlen = {chrom: len(ff[chrom]) for chrom in chroms}


def generate_random_coordinates(num, size=1000):
    chroms = [str(i) for i in range(1, 22)]
    
    for _ in range(num):
        chrom = random.choice(chroms)
        start = np.random.randint(0, chromlen[chrom] - size)
        end = start + size
        yield chrom, start, end

def benchmark_meth5(f):
    rand_coords = list(generate_random_coordinates(100))
    st = time.time()
    num = 0
    for coord in rand_coords:
        num += len(f[coord[0]].get_values_in_range(coord[1], coord[2]).get_llrs())
    return time.time() - st

def benchmark_modbam(f):
    rand_coords = list(generate_random_coordinates(100))
    st = time.time()
    num = 0
    for coord in rand_coords:
        with ModBam(f, coord[0], coord[1], coord[2]) as bam:
            for read in bam.reads():
                for pos_mod in read.mod_sites:
                    ref_pos = pos_mod[1]
                    if coord[1] <= ref_pos < coord[2]:
                        mod_base_score = pos_mod[7]
                        print(mod_base_score)
                        num += 1
    return time.time() - st


rounds = 100
times_meth5_lzf = [benchmark_meth5(mf_lzf) for _ in tqdm.tqdm(list(range(rounds)))]
times_meth5_gz = [benchmark_meth5(mf_gz) for _ in tqdm.tqdm(list(range(rounds)))]
times_modbam = [benchmark_modbam(bf) for _ in tqdm.tqdm(list(range(rounds)))]
times_modcram = [benchmark_modbam(cf) for _ in tqdm.tqdm(list(range(rounds)))]

with pa.open_multipage_pdf("meth5_modbampy_random_access"):
    pa.figure(figsize=(6, 3))
    plt.title("Seconds accessing 100 random 1000bps spans (seconds)")
    plt.barh(1, np.mean(times_meth5_lzf), xerr=np.std(times_meth5_lzf))
    plt.barh(2, np.mean(times_meth5_gz), xerr=np.std(times_meth5_gz))
    plt.barh(3, np.mean(times_modbam), xerr=np.std(times_modbam))
    plt.barh(4, np.mean(times_modcram), xerr=np.std(times_modcram))
    plt.yticks([1, 2, 3, 4], labels=["MetH5 (LZF compressed)", "MetH5 (gzip compressed)", "BAM", "CRAM"])
    plt.xlabel("Time (seconds)")
    pa.savefig()


def sequential_meth5(f, chrom):
    st = time.time()
    pos = 0
    num = 0
    for chunk in range(f[chrom].get_number_of_chunks()):
        llrs = f[chrom].get_chunk(chunk).get_llrs()
        pos += (llrs > 2).sum()
        num += (np.abs(llrs) > 2).sum()
    print(pos / num)
    return time.time() - st


def sequential_modbam(f, chrom):
    st = time.time()
    pos = 0
    count = 0
    with ModBam(f, chrom, 0, chromlen[chrom]) as bam:
        for read in bam.reads():
            for pos_mod in read.mod_sites:
                mod_base_score = pos_mod[7]
                if mod_base_score > 168:
                    pos += 1
                    count += 1
                if mod_base_score < 85:
                    count += 1
    print(pos / count)
    return time.time() - st


seconds_meth5_chrom_lzf = sequential_meth5(mf_lzf,"21")
seconds_meth5_chrom_gz = sequential_meth5(mf_gz,"21")
seconds_modbam_chrom = sequential_modbam(bf,"21")
seconds_mocram_chrom = sequential_modbam(cf,"21")


with pa.open_multipage_pdf("meth5_modbampy_sequential_access"):
    pa.figure(figsize=(6, 3))
    plt.title("Seconds computing total methylation rate of chromosome 21")
    plt.barh(1, seconds_meth5_chrom_lzf)
    plt.barh(2, seconds_meth5_chrom_gz)
    plt.barh(3, seconds_modbam_chrom)
    plt.barh(4, seconds_mocram_chrom)
    plt.yticks([1, 2, 3, 4], labels=["MetH5 (LZF compressed)", "MetH5 (gzip compressed)", "BAM", "CRAM"])
    plt.xlabel("Time (seconds)")
    pa.savefig()
