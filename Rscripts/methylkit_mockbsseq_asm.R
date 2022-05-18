library(methylKit)

basedir="/data/r933r/giab/mock_bsseq_from_nanopore/"

for (sample in c("HG002", "HG003", "HG004")) {
    bsfile = paste0(basedir, sample, "_mock_bsseq.bedGraph")
    myobj = methRead(bsfile, assembly="hg38", sample.id=sample, treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
    seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
    methSeg2bed(seg, filename=paste0(basedir, sample, "_mockbsseq_segments_methylkit.tsv"))
}