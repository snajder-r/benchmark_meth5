library(methylKit)
library(bsseq)
library(GenomicRanges)
library(MethCP)
library(BiocParallel)


bsfiles = c("/data/r933r/giab/mock_bsseq_from_nanopore/simulated/diffmet_merged_DCASES1.tsv",
                 "/data/r933r/giab/mock_bsseq_from_nanopore/simulated/diffmet_merged_DCONTROLS1.tsv")

bs_object <- createBsseqObject(
    files = bsfiles, sample_names = list("DCASES1", "DCONTROLS1"),
    chr_col = 1, pos_col = 2, m_col = 6, cov_col = 5, header = FALSE)

obj_methylKit = calcLociStat(bs_object, "DCASES1", "DCONTROLS1", test="methylKit")

# ========================================================================
# << << << << ABOVE IS MODIFIED VERSION OF calcLociStat << << << <<
# ========================================================================

bpparams=bpparam()
bpparams$workers = 4
obj_seg_methylKit <- segmentMethCP(obj_methylKit, bs_object, region.test = "fisher", BPPARAM=bpparams)

write.table(data.frame(obj_seg_methylKit@segmentation), file="/data/r933r/giab/mock_bsseq_from_nanopore/simulated/methcp_segments.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

region_methylKit <- getSigRegion(obj_seg_methylKit)

write.table(region_methylKit, file="/data/r933r/giab/mock_bsseq_from_nanopore/simulated/methcp_dmr.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)