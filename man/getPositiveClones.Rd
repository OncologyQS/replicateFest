% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manafest_shiny_functions.r
\name{getPositiveClones}
\alias{getPositiveClones}
\title{getPositiveClones}
\usage{
getPositiveClones(
  analysisRes,
  mergedData,
  samp = names(mergedData),
  orThr = 1,
  fdrThr = 0.05,
  percentThr = 0
)
}
\arguments{
\item{analysisRes}{a list with results of Fisher's test
for pairwise comparisons}

\item{mergedData}{a list of data frames with read counts
for each sample}

\item{samp}{a vector with sample IDs to analyze}

\item{orThr}{a threshold for odds ratio}

\item{fdrThr}{a threshold for FDR}

\item{percentThr}{a threshold for the percentage of reads}

\item{nReads}{a threshold for the number of reads}
}
\value{
a vector with positive clones as names and conditions,
in which a clone is significant, as values
}
\description{
returns a vector with positive clones as names and conditions,
in which a clone is significant, as values.
For every clone, checks top two conditions with the highest number of reads
to see if a clone is uniquely expanded
}
