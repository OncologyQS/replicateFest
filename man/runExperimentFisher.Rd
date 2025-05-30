% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manafest_shiny_functions.r
\name{runExperimentFisher}
\alias{runExperimentFisher}
\title{runExperimentFisher}
\usage{
runExperimentFisher(
  files,
  refSamp,
  nReads = 50,
  fdrThr = 0.05,
  orThr = 5,
  percentThr = 0,
  excludeSamp = "",
  compareToRef = TRUE,
  outputFile = "output.xlsx",
  saveToFile = T
)
}
\arguments{
\item{files}{a list of file names with read counts}

\item{refSamp}{a reference sample ID}

\item{nReads}{a threshold for the number of reads}

\item{fdrThr}{a threshold for FDR}

\item{orThr}{a threshold for OR}

\item{percentThr}{a threshold for percentage}

\item{excludeSamp}{a sample ID to exclude from analysis}

\item{compareToRef}{a logical value indicating if the comparison to the reference sample should be performed}

\item{outputFile}{a name of the output file}

\item{saveToFile}{a logical value indicating if the output should be saved to a file}
}
\value{
a list with output tables
}
\description{
Reads in files with TCR repertoires from
a FEST experiment without replicate samples. It
performs pairwise Fisher's test find expanded clones
comparing to a reference samples.
It also compares top conditions to find unique expansions.
The results are return and saved in an Excel file.
}
