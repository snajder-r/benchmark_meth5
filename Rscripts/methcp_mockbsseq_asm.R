library(methylKit)
library(bsseq)
library(GenomicRanges)
library(MethCP)
library(BiocParallel)

basedir="/data/r933r/giab/mock_bsseq_from_nanopore/"

for (sample in c("HG003", "HG004")) {
    bsfiles = c(paste0(basedir,sample,"_H1_mock_bsseq.bedGraph"), paste0(basedir, sample, "_H2_mock_bsseq.bedGraph"))

    bs_object <- createBsseqObject(files = bsfiles, sample_names = list(paste0(sample, "_H1"), paste0(sample, "_H2")), chr_col = 1, pos_col = 2, m_col = 6, cov_col = 5, header = FALSE)
    obj_methylKit = calcLociStat(bs_object, paste0(sample, "_H1"), paste0(sample, "_H2"), test="methylKit")

    bpparams=bpparam()
    bpparams$workers = 4
    obj_seg_methylKit <- segmentMethCP(obj_methylKit, bs_object, region.test = "fisher", BPPARAM=bpparams)

    write.table(data.frame(obj_seg_methylKit@segmentation), file=paste0(basedir, sample, "_asm_methcp_segments.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    region_methylKit <- getSigRegion(obj_seg_methylKit, sig.level=0.5, nC.valid=5)
    write.table(region_methylKit, file=paste0(basedir, sample,"_asm_methcp_highpval_dmr.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
