library(methylKit)

myobj = methRead("/data/r933r/giab/mock_bsseq_from_nanopore/simulated/diffmet_merged_15x_both.tsv", assembly="hg38", sample.id="simulation", treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
methSeg2bed(seg,filename="/data/r933r/giab/mock_bsseq_from_nanopore/simulated/methylkit_segmentation_15x.bed")
