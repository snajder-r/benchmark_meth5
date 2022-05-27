library(methylKit)
library(fastseg)
hp12 = methRead(list("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H1_mock_bsseq.bedGraph",
                     "/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H2_mock_bsseq.bedGraph",
                     "/data/r933r/giab/mock_bsseq_from_nanopore/HG004_H1_mock_bsseq.bedGraph",
                     "/data/r933r/giab/mock_bsseq_from_nanopore/HG004_H2_mock_bsseq.bedGraph"),
                mincov=5, assembly="hg38", sample.id=list("HG003", "HG003", "HG004", "HG004"), treatment=c(0,0,1,1),
                header=FALSE, context="CpG", resolution="base", pipeline=
                  list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))


united=unite(hp12)
diff = calculateDiffMeth(united)
diffseg = methSeg(diff)
methSeg2bed(diffseg,filename="/data/r933r/giab/mock_bsseq_from_nanopore/parents_mockbsseq_diff_methylkit_seg.bed")
