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
pa = PlotArchiver("meth5", config={"plot_archive_dir": "/home/r933r/snajder/nanoepitools_plots/benchmark"})

mbs = 1024 * 1024
file_size_hg003_modbam = (88158512983) / mbs
file_size_hg003_modcram = (56391576891) / mbs
file_size_hg003_bam = (86750913912) / mbs
file_size_hg003_cram = (55561487882) / mbs

file_size_hg003_meth5 = (1580536213) / (1024 * 1024)
file_size_hg003_meth5_noreadnames = 1402323813 / (1024 * 1024)


with pa.open_multipage_pdf("meth5_modbampy_filesize"):
    pa.figure(figsize=(8, 3))
    plt.title("File size for ~30x human genome")
    plt.barh(1, file_size_hg003_meth5, color="#FC802D")
    plt.barh(3, file_size_hg003_bam, color="#2976B1")
    plt.barh(3, file_size_hg003_bam - file_size_hg003_modbam, left=file_size_hg003_bam, color="#FC802D",)
    plt.barh(2, file_size_hg003_cram, color="#2976B1")
    plt.barh(2, file_size_hg003_cram - file_size_hg003_modcram, left=file_size_hg003_cram, color="#FC802D")
    plt.yticks([1, 2, 3], labels=["MetH5", "ModCRAM", "ModBAM"])
    plt.xlabel("Megabytes")
    pa.savefig()
