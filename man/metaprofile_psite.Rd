% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaplots.R
\name{metaprofile_psite}
\alias{metaprofile_psite}
\title{Plot ribosome occupancy metaprofiles at single-nucleotide resolution.}
\usage{
metaprofile_psite(data, annotation, sample, scale_factors = NULL,
  length_range = "all", transcripts = NULL, utr5l = 25, cdsl = 50,
  utr3l = 25, plot_title = "auto")
}
\arguments{
\item{data}{A list of data frames from \code{\link{psite_info}}.}

\item{annotation}{A data frame as generated by
\code{\link{create_annotation}}.}

\item{sample}{A character string vector specifying the name of the sample (or
of its replicates) of interest. Its elements are merge together using the
scale factors specified by \code{scale_factors}.}

\item{scale_factors}{A numeric vector of scale factors for merging the
replicates (if any). The vector must contain at least a set of values, one
for each replicates, named after the strings listed in \code{sample}. No
specific order is required. Default is NULL, meaning that all the scale
factors are set to 1.}

\item{length_range}{Either "all", an integer or an integer vector. Default is
"all", meaning that all the read lengths are included in the analysis.
Otherwise, only the read lengths matching the specified value(s) are kept.}

\item{transcripts}{A character string vector specifying the name of the
transcripts to be included in the analysis. By default this argument is
NULL, meaning that all the transcripts in \code{data} are used. Note that
if either the 5' UTR, the coding sequence or the 3' UTR of a transcript is
shorther than \code{utr5l}, \eqn{2*}\code{cdsl} and \code{utr3l}
respectively, the transcript is automatically discarded.}

\item{utr5l}{A positive integer specifying the length (in nucleotides) of the
5' UTR region that in the plot flanks the start codon. The default value is
25.}

\item{cdsl}{A positive integer specifying the length (in nucleotides) of the
CDS region that in the plot will flank both the start and the stop codon.
The default value is 50.}

\item{utr3l}{A positive integer specifying the length (in nucleotides) of the
3' UTR region that in the plot flanks the start codon. The default value is
25.}

\item{plot_title}{Either "auto", NULL or any character string. The default is
"auto", meaning that the title of the plot is the name of the sample
followed by the number of the transcripts and the read lengths considered
for generating the metaprofile.}
}
\value{
A list containing a ggplot2 object, a data frame with the associated
  data and the transcripts employed for generating the plot.
}
\description{
For a specified sample this function generates a metaprofile based on the
P-site of the reads mapping around the start and the stop codon of the
annotated CDS (if any). It sums up the number of P-sites (defined by their
first nucleotide) per nucleotide computed for all the transcripts starting
from one ore more replicates.
}
\examples{
data(reads_psite_list)
data(mm81cdna)

## Generate the metaprofile for all the read lengths
metaprof_whole <- metaprofile_psite(reads_psite_list, mm81cdna, sample = "Samp1")
metaprof_whole[["plot"]]

## Generate the metaprofile employing reads of 27, 28 and 29 nucleotides and
a subset of transcripts (for example with at least one P-site mapping on
the translation initiation site)
sample_name <- "Samp1"
sub_reads_psite_list <- subset(reads_psite_list[[sample_name]], psite_from_start == 0)
transcript_names <- as.character(sub_reads_psite_list$transcript)
metaprof_sub <- metaprofile_psite(reads_psite_list, mm81cdna, sample = sample_name,
length_range = 27:29, transcripts = transcript_names)
metaprof_sub[["plot"]]
}
