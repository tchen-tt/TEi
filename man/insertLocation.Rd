% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{insertLocation}
\alias{insertLocation}
\title{Identify the TE insertions}
\usage{
insertLocation(file, outBed, ratio = 0.1)
}
\arguments{
\item{file}{Charcater, the alignment file in bam/sam format, out from \bold{alignment}}

\item{outBed}{Character, the bed file for insertions.}

\item{ratio}{Supports the softclip reads ratio}
}
\description{
Identify the locations of TE insertions from the newly alignment result.
}
