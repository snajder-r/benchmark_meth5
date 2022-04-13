library(methylKit)
library(bsseq)
library(GenomicRanges)
library(MethCP)
library(BiocParallel)


bsfiles = c("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq.bedGraph",
                 "/data/r933r/giab/mock_bsseq_from_nanopore/HG004_mock_bsseq.bedGraph")

bs_object <- createBsseqObject(
    files = bsfiles, sample_names = list("HG003", "HG004"),
    chr_col = 1, pos_col = 2, m_col = 6, cov_col = 5, header = FALSE)

obj_methylKit = calcLociStat(bs_object, "HG003", "HG004", test="methylKit")


bpparams=bpparam()
bpparams$workers = 4
obj_seg_methylKit <- segmentMethCP(obj_methylKit, bs_object, region.test = "fisher", BPPARAM=bpparams)
