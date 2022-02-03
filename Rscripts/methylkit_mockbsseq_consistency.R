library(methylKit)

parent_dir = "/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq_small/"
subset_files = list.files(parent_dir)

out_file = "/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq_small/"
subset_files = list.files(parent_dir)

coords = c()
for (file in subset_files) {
  coords = union(coords, read.table(paste0(parent_dir,file))[, 2])
}

coords = sort(coords)
frequencies = rep(0, length(coords))

for (file in subset_files) {
  myobj = methRead(paste0(parent_dir,file), assembly="hg38", sample.id="HG003", treatment=0, header=FALSE, context="CpG", resolution="base", pipeline=list(fraction=TRUE, chr.col=1, start.col=2, strand.col=7, end.col=3, coverage.col=5, freqC.col=4), mincov=1)
  seg = methSeg(myobj, diagnostic.plot=FALSE, maxInt=10, minSeg=10, alpha=30)

  changepoints = seg@ranges@start
  for (cp in changepoints){
    idx = which(coords==cp)
    frequencies[idx] = frequencies[idx] + 1
  }
}



out_df = DataFrame(chrom="10", pos=coords, changepoint_frequency=frequencies)
write.table(out_df, "/data/r933r/giab/mock_bsseq_from_nanopore/HG003_small_example_consistency.tsv", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)