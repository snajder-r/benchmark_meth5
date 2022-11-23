library(methylKit)
library(fastseg)
myobj = methRead("/data/r933r/giab/mock_bsseq_from_nanopore/parents_mock_bsseq.bedGraph",
                mincov=5, assembly="hg38", sample.id="parents", treatment=0,
                header=FALSE, context="CpG", resolution="base", pipeline=
                  list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))



seg = methSeg(myobj)
methSeg2bed(seg, filename="/data/r933r/giab/mock_bsseq_from_nanopore/parents_mockbsseq_segments_methylkit.bed")
