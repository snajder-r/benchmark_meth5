library(methylKit)
library(bsseq)
library(GenomicRanges)
library(MethCP)
library(BiocParallel)


bsfiles = c("/data/r933r/giab/mock_bsseq_from_nanopore/HG003_mock_bsseq.bedGraph",
                 "/data/r933r/giab/mock_bsseq_from_nanopore/HG004_mock_bsseq.bedGraph")

bs_object <- createBsseqObject(
    files = bsfiles, sample_names = list("HG003", "HG004"),
    chr_col = 1, pos_col = 2, m_col = 6, cov_col = 5, header = FALSE)

obj_methylKit = calcLociStat(bs_object, "HG003", "HG004", test="methylKit")

# Unfortunately, for me the function calcLociStat would throw an obscure error, saying that a parallel job did not
# return properly. Therefore I modify it here to run sequentially, which seems to work.
# And no, I did not spend time trying to make it pretty. Already wasted enough time debugging it.
# Original source: https://github.com/boyinggong/MethCP/blob/405650ac77c6fbc9886b9f4b3330f00f0489a600/R/calcLociStat.R
# ========================================================================
# >> >> >> >> BELOW IS MODIFIED VERSION OF calcLociStat >> >> >> >>
# ========================================================================

bs.object = bs_object
group1 = "HG003"
group2 = "HG004"

test = "methylKit"

object_list <- list()
for (chr in unique(seqnames(bs.object))){
	object_list[[chr]] <- bs.object[seqnames(bs.object) == chr]
}

nsample <- length(group1) + length(group2)

statlist = list()

for (o in object_list){
	print(o)
	df <- cbind(
		as.data.frame(granges(o)),
		getCoverage(o, type = "Cov")[, c(group1, group2)],
		getCoverage(o, type = "M")[, c(group1, group2)],
		(getCoverage(o, type = "Cov") - getCoverage(
			o, type = "M"))[, c(group1, group2)])
	df$width <- NULL
	colnames(df) <- c(
		"chr", "start", "end", "strand",
		paste0("coverage", seq_len(nsample)),
		paste0("numCs", seq_len(nsample)),
		paste0("numTs", seq_len(nsample)))
	df <- df[
		, c("chr", "start", "end", "strand",
			vapply(seq_len(nsample), function(x)
				paste0(c("coverage", "numCs", "numTs"), x),
				FUN.VALUE = character(3)))]
	coverage.ind <- which(vapply(colnames(df), function(x)
		grepl("coverage", x), FUN.VALUE = logical(1)))
	names(coverage.ind) <- NULL
	obj <- new(
		"methylBase", (df),
		sample.ids = c(group1, group2),
		assembly = "-",
		context = "-",
		treatment = c(
			rep(1, length(group1)),
			rep(0, length(group2))),
		coverage.index = coverage.ind,
		numCs.index = coverage.ind+1,
		numTs.index = coverage.ind+2,
		destranded = TRUE,
		resolution = "base" )
	filter <- (rowSums(as.data.frame(getCoverage(o)[, group1])) >= 1) &
		(rowSums(as.data.frame(getCoverage(o)[, group2])) >= 1)
	errorInd <- rowSums(as.data.frame(getCoverage(o, type = "M"))) /
		rowSums(as.data.frame(getCoverage(o)))
	errorInd <- (errorInd == 0) | (errorInd == 1)
	errorInd[is.na(errorInd)] = FALSE
	obj_filtered <- obj[filter & !errorInd, ]
	invisible(capture.output(tmp <- calculateDiffMeth(obj_filtered)))
	gr <- GRanges(
		seqnames = tmp$chr,
		IRanges(start = tmp$start, end = tmp$start),
		pval = tmp$pval, methDiff = tmp$meth.diff/100)
	if(sum(errorInd) != 0){
		gr0 <- granges(o[which(errorInd)])
		gr0$pval = 1
		gr0$methDiff = 0
		gr <- unlist(GRangesList(gr, gr0))
		gr <- sortSeqlevels(gr)
		gr <- sort(gr)
	}
	statlist = append(statlist, list(gr))
}
names(statlist) <- unique(seqnames(bs.object))
statlist <- GenomicRanges::GRangesList(statlist)

obj_methylKit <- new(
	"MethCP",
	group1 = group1,
	group2 = group2,
	test = test,
	stat = statlist)

# ========================================================================
# << << << << ABOVE IS MODIFIED VERSION OF calcLociStat << << << <<
# ========================================================================

bpparams=bpparam()
bpparams$workers = 4
obj_seg_methylKit <- segmentMethCP(obj_methylKit, bs_object, region.test = "fisher", BPPARAM=bpparams)
