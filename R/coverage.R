#' Compute the number of reads per codon.
#'
#' For the specified sample(s), this function computes the codon coverage
#' defined either as the number of read footprints per codon or as the number of
#' P-sites per codon.
#'
#' @param data A list of data frames from \code{\link{psite_info}}. Data frames
#'   generated by \code{\link{bamtolist}} and \code{\link{bedtolist}} can be
#'   used only if \code{psite} is FALSE (the default).
#' @param annotation A data frame as generated by
#'   \code{\link{create_annotation}}.
#' @param sample A character string vector specifying the name of the sample(s)
#'   of interest. By default this argument is NULL, meaning that the coverage is
#'   computed for all the samples in \code{data}.
#' @param psite A logical value whether or not to return the number of P-sites
#'   per codon. Default is NULL, meaning that the number of read footprints per
#'   codon is computed instead.
#' @param min_overlap A positive integer specyfing the minimum number of
#'   overlapping positions (in nucleotides) between a reads and a codon to be
#'   considered to be overlapping. When \code{psite} is TRUE this parameter must
#'   be 1 (the default).
#' @param granges A logical value whether or not to return a GRanges object.
#'   Default is FALSE, meaning that a data frames is returned instead.
#' @details The sequence of every transcript is divided in triplets starting
#'   from the annotated translation initiation site (if any) proceeding towards
#'   the UTRs extremities, and eventually discarding the exceeding 1 or 2
#'   nucleotides at the extremities of the transcript. Please note that the
#'   transcripts not associated to any annotated \emph{5' UTR}, \emph{CDS} and
#'   \emph{3'UTR} and transcripts with coding sequence length not divisible by 3
#'   are automatically discarded.
#' @return A data frame or a GRanges object.
codon_coverage <- function(data, annotation, sample = NULL, psite = FALSE,
                         min_overlap = 1, granges = FALSE) {
  
  if((psite == TRUE || psite == T) & min_overlap != 1){
    cat("\n")
    stop("\nERROR: invalid value for min_over when psite is TRUE\n\n")
  }
  
  bin <- 3
  
  cat("1. creating codon table\n")
  rownames(annotation) <- as.character(annotation$transcript)
  l.transcripts <- rownames(annotation)[which(annotation$l_utr5 > 0 &
                                                annotation$l_cds >0 &
                                                annotation$l_cds %% 3 == 0 &
                                                annotation$l_utr3 > 0)]
  
  subanno <- subset(annotation, as.character(transcript) %in% l.transcripts)
  subanno <- subanno[order(as.character(subanno$transcript)), ]
  subanno$transcript <- factor(subanno$transcript, levels = subanno$transcript)
  subanno$start <- (subanno$l_utr5) %% bin
  
  for (coln in c("l_utr5", "l_cds", "l_utr3")) {
    subanno[, coln] <- subanno[, coln] - (subanno[, coln] %% bin)
  }
  subanno$len <- subanno$l_utr5 + subanno$l_cds + subanno$l_utr3
  subanno$stop <- subanno$start + subanno$len
  
  bin_coverage_tab <- data.frame(transcript = rep(subanno$transcript, (subanno$stop) / bin),
                                 start = as.vector(unlist(list(by(subanno,
                                                                  subanno$transcript,
                                                                  function(x) seq(from = x$start, to = x$stop - bin, by = bin))))))
  bin_coverage_tab$end = bin_coverage_tab$start + bin
  
  gr_interval <- GenomicRanges::GRanges(seqnames = bin_coverage_tab$transcript,
                                        IRanges(bin_coverage_tab$start + 1,  width = bin),
                                        strand="+")
  
  cat("2. computing distance from start/stop codon\n")
  bin_coverage_tab$start_dist <- as.vector(unlist(list(by(subanno, subanno$transcript, function(x) rep(x$start + x$l_utr5, x$len / bin)))))
  bin_coverage_tab$start_dist <- (bin_coverage_tab$start - bin_coverage_tab$start_dist) / bin
  bin_coverage_tab$stop_dist <- as.vector(unlist(list(by(subanno, subanno$transcript, function(x) rep(x$start + x$l_utr5 + x$l_cds - bin, x$len/bin)))))
  bin_coverage_tab$stop_dist <- (bin_coverage_tab$start - bin_coverage_tab$stop_dist) / bin
  
  cat("3. acquiring region information\n")
  bin_coverage_tab$region <- as.vector(unlist(list(by(subanno, subanno$transcript, function(x) rep(c("5utr", "cds", "3utr"), times = c(x$l_utr5 / bin, x$l_cds / bin, x$l_utr3 / bin))))))
  
  if(psite == T || psite == TRUE){
    cat("4. computing codon coverage based on P-site\n")
  } else {
    cat("4. computing codon coverage based on read footprints\n")
  }
  for(samp in sample){
    cat(sprintf("%s\n", samp))
    df <- subset(data[[samp]], as.character(transcript) %in% subanno$transcript)
    df$transcript <- factor(df$transcript, levels = subanno$transcript)
    
    if(psite == T || psite == TRUE){
      gr_read <- GenomicRanges::GRanges(seqnames = df$transcript,
                         IRanges(df$psite, df$psite),
                         strand="+")
    } else {
      gr_read <- GenomicRanges::GRanges(seqnames = df$transcript,
                         IRanges(df$end5, df$end3),
                         strand="+")
    }
    
    bin_coverage_tab[, samp] <- GenomicRanges::countOverlaps(gr_interval, gr_read, minoverlap=min_overlap)
  }
  
  if (granges == T || granges == TRUE) {
    bin_coverage_tab <- GenomicRanges::makeGRangesFromDataFrame(bin_coverage_tab,
                                                  keep.extra.columns=TRUE,
                                                  ignore.strand=TRUE,
                                                  seqnames.field=c("transcript"),
                                                  start.field="end5",
                                                  end.field="end3",
                                                  strand.field="strand",
                                                  starts.in.df.are.0based=FALSE)
    strand(bin_coverage_tab) <- "+"
  }
  
  return(bin_coverage_tab)
}
