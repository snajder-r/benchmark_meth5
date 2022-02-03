import time
import random
from modbampy import ModBam, pileup
import pyfaidx
from meth5.meth5 import MetH5File
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from nanoepitools.plotting.general_plotting import PlotArchiver
from benchmark_pycometh.config import module_config

matplotlib.use("Agg")
chroms = [str(i) for i in range(1, 22)]

pa = PlotArchiver("meth5", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

mf = MetH5File(module_config.meth5_template_file.format(sample="HG003"), "r", chunk_size=50000)

bf = module_config.modcram_template_file.format(sample="HG003")
ff = pyfaidx.Fasta(module_config.reference)
chromlen = {chrom: len(ff[chrom]) for chrom in chroms}


def generate_random_coordinates(num, size=1000):
    chroms = [str(i) for i in range(1, 22)]
    
    for _ in range(num):
        chrom = random.choice(chroms)
        start = np.random.randint(0, chromlen[chrom] - size)
        end = start + size
        yield chrom, start, end

def benchmark_meth5():
    rand_coords = list(generate_random_coordinates(100))
    st = time.time()
    num = 0
    for coord in rand_coords:
        num += len(mf[coord[0]].get_values_in_range(coord[1], coord[2]).get_llrs())
    return time.time() - st

def benchmark_modbam():
    rand_coords = list(generate_random_coordinates(100))
    st = time.time()
    num = 0
    for coord in rand_coords:
        with ModBam(bf, coord[0], coord[1], coord[2]) as bam:
            for read in bam.reads():
                for pos_mod in read.mod_sites():
                    ref_pos = pos_mod[1]
                    if coord[1] <= ref_pos < coord[2]:
                        mod_base_score = pos_mod[7]
                        num += 1
    return time.time() - st


rounds = 100
times_meth5 = [benchmark_meth5() for _ in range(rounds)]
times_modbam = [benchmark_modbam() for _ in range(rounds)]

with pa.open_multipage_pdf("meth5_modbampy_random_access"):
    pa.figure()
    plt.title("MetH5 random access")
    plt.hist(times_meth5, bins=15)
    plt.xlabel("Seconds accessing 100 random 1000bps spans (seconds)")
    pa.savefig()
    
    pa.figure()
    plt.title("Modbampy random access")
    plt.hist(times_modbam, bins=15)
    plt.xlabel("Seconds accessing 100 random 1000bps spans (seconds)")
    pa.savefig()
    
    pa.figure(figsize=(8, 3))
    plt.title("Seconds accessing 100 random 1000bps spans (seconds)")
    plt.barh(1, np.mean(times_meth5), xerr=np.std(times_meth5))
    plt.barh(2, np.mean(times_modbam), xerr=np.std(times_modbam))
    plt.yticks([1, 2], labels=["MetH5", "Modbampy"])
    plt.xlabel("Time (seconds)")
    pa.savefig()


def sequential_meth5(chrom):
    st = time.time()
    pos = 0
    num = 0
    for chunk in range(mf[chrom].get_number_of_chunks()):
        llrs = mf[chrom].get_chunk(chunk).get_llrs()
        pos += (llrs > 2).sum()
        num += (np.abs(llrs) > 2).sum()
    print(pos / num)
    return time.time() - st


def sequential_modbam(chrom):
    st = time.time()
    pos = 0
    count = 0
    with ModBam(bf, chrom, 0, chromlen[chrom]) as bam:
        for read in bam.reads():
            for pos_mod in read.mod_sites():
                mod_base_score = pos_mod[7]
                if mod_base_score > 168:
                    pos += 1
                    count += 1
                if mod_base_score < 85:
                    count += 1
    print(pos / count)
    return time.time() - st


seconds_modbam_chrom = sequential_modbam("21")
seconds_meth5_chrom = sequential_meth5("21")

with pa.open_multipage_pdf("meth5_modbampy_sequential_access"):
    pa.figure(figsize=(8, 3))
    plt.title("Seconds computing total methylation rate of chromosome 21")
    plt.barh(1, np.mean(seconds_meth5_chrom))
    plt.barh(2, np.mean(seconds_modbam_chrom))
    plt.yticks([1, 2], labels=["MetH5", "Modbampy"])
    plt.xlabel("Time (seconds)")
    pa.savefig()
