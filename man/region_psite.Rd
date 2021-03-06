% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/percentage_regions.R
\name{region_psite}
\alias{region_psite}
\title{Plot the percentage of P-sites per transcript region.}
\usage{
region_psite(data, annotation, sample = NULL, transcripts = NULL,
  label = NULL, colour = c("gray70", "gray40", "gray10"))
}
\arguments{
\item{data}{A list of data frames from \code{\link{psite_info}}.}

\item{annotation}{A data frame as generated by
\code{\link{create_annotation}}.}

\item{sample}{A character string vector specifying the name of the sample(s)
of interest. By default this argument is NULL, meaning that all the samples
in \code{data} are included in the analysis.}

\item{transcripts}{A character string vector specifying the name of the
transcripts to be included in the analysis. By default this argument is
NULL, meaning that all the transcripts in \code{data} are used. Please note
that the transcripts not associated to any annotated \emph{5' UTR},
\emph{CDS} and \emph{3'UTR} are automatically discarded.}

\item{label}{A character string vector of the same length of \code{sample}
specifying the name of the samples to be displaied in the plot. By default
this argument is NULL meaning that name of the samples are used.}

\item{colour}{A character string vector of three elements specifying the
colours of the bars corresponding to the \emph{5' UTR}, the \emph{CDS} and
the \emph{3'UTR} respectively. The default is a grayscale.}
}
\value{
A list containing a ggplot2 object, and a data frame with the
  associated data.
}
\description{
For one or several samples this function computes the percentage of P-sites
falling in the three annotated regions of the transcripts (5' UTR, CDS and
3'UTR) and generates a barplot of the resulting values. The function also
calculates and plots the percentage of region length for the selected
transcripts (reported in column "RNAs").
}
\examples{
data(reads_psite_list)
data(mm81cdna)

reg_psite <- region_psite(reads_psite_list, mm81cdna, sample = "Samp1")
reg_psite[["plot"]]
}
