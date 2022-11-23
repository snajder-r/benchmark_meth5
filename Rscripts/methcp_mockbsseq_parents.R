library(methylKit)
library(bsseq)
library(GenomicRanges)
library(MethCP)
library(BiocParallel)

bsfiles = c("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq.bedGraph",
                 "/data/r933r/giab/mock_bsseq_from_nanopore/HG004_mock_bsseq.bedGraph")

out_seg_file = "/data/r933r/giab/mock_bsseq_from_nanopore/parents_mockbsseq_segments_methcp.tsv"
out_diffmet_file = "/data/r933r/giab/mock_bsseq_from_nanopore/parents_mockbsseq_segments_methcp_diffmet.tsv"

bs_object <- createBsseqObject(
    files = bsfiles, sample_names = list("HG003", "HG004"),
    chr_col = 1, pos_col = 2, m_col = 6, cov_col = 5, header = FALSE)

obj_methylKit = calcLociStat(bs_object, "HG003", "HG004", test="methylKit")

bpparams=bpparam()
bpparams$workers = 4
obj_seg_methylKit <- segmentMethCP(obj_methylKit, bs_object, region.test = "fisher", BPPARAM=bpparams)
write.table(data.frame(obj_seg_methylKit@segmentation), file=out_seg_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
region_methylKit <- getSigRegion(obj_seg_methylKit, sig.level=0.5, nC.valid=5)
write.table(region_methylKit, file=out_diffmet_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
