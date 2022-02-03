library(methylKit)

myobj = methRead("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq.bedGraph", assembly="hg38", sample.id="HG003", treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
seg = methSeg(myobj, diagnostic.plot=TRUE, maxInt=100, minSeg=10)
methSeg2bed(seg,filename="/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mockbsseq_seg.bed")

hp1 = methRead("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H1_mock_bsseq.bedGraph", assembly="hg38", sample.id="HG003", treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))
hp2 = methRead("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H2_mock_bsseq.bedGraph", assembly="hg38", sample.id="HG003", treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))


hp12 = methRead(list("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H1_mock_bsseq.bedGraph","/data/r933r/giab/mock_bsseq_from_nanopore/HG003_H2_mock_bsseq.bedGraph"), mincov=5, assembly="hg38", sample.id=list("HG003", "HG003"), treatment=c(0,0), header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4))


united=unite(hp12)
frac1 = united$numCs1/united$coverage1
frac2 = united$numCs2/united$coverage2
m = matrix(c(frac1, frac2), ncol=2)

seg = fastseg(m, cyberWeight=100)


chroms = c()
starts = c()
ends = c()
out_i = 1
for (row in 1:length(seg)){
  row = seg@ranges[row]
  range_row_start = united[row@start]
  start = range_row_start$start
  range_row_end = united[row@start + row@width -1]
  end = range_row_end$end
  chrom_a  = range_row_start$chr
  chrom_b  = range_row_end$chr
  if (chrom_a == chrom_b){
    chroms[out_i] = chrom_a
    starts[out_i] = start
    ends[out_i] = end
    out_i = out_i + 1
  }else{
    print(c(chrom_a, chrom_b))
    chroms[out_i] = chrom_a
    starts[out_i] = start
    ends[out_i] = start + 10000
    out_i = out_i + 1
    chroms[out_i] = chrom_b
    starts[out_i] = 0
    ends[out_i] = end
    out_i = out_i + 1
  }
}

seg_df = DataFrame(chrom=chroms,start=starts, end=ends)
write.table(seg_df, "/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mockbsseq_hps_fastseg.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)