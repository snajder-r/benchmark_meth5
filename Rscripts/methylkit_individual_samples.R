library(methylKit)

myobj = methRead("/data/r933r/giab/HG002_rep1_rep2_merged.bedGraph.gz", assembly="hg38", sample.id="HG002", treatment=0, header=FALSE, context="CpG",
resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
methSeg2bed(seg,filename="/data/r933r/giab/HG002_rep1_rep2_seg.bed")

myobj = methRead("/data/r933r/giab/HG003_rep1_rep2_merged.bedGraph.gz", assembly="hg38", sample.id="HG003", treatment=0, header=FALSE, context="CpG",
resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
methSeg2bed(seg,filename="/data/r933r/giab/HG003_rep1_rep2_seg.bed")

myobj = methRead("/data/r933r/giab/HG004_rep1_rep2_merged.bedGraph.gz", assembly="hg38", sample.id="HG004", treatment=0, header=FALSE, context="CpG",
resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
methSeg2bed(seg,filename="/data/r933r/giab/HG004_rep1_rep2_seg.bed")
