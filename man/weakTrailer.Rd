% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassifierHelpers.R
\name{weakTrailer}
\alias{weakTrailer}
\title{Find non-active trailers}
\usage{
weakTrailer(coverage, fpkm, startRegionRelative, ORFScore)
}
\arguments{
\item{coverage}{(counts over region)}

\item{fpkm}{fragments per kilobase transcript per million}

\item{startRegionRelative}{Relative coverage of start region to upstream short region}

\item{ORFScore}{Periodicity score of triplets, > 0 if frame 0 has most reads.}
}
\value{
logical TRUE/FALSE per row
}
\description{
A filter to say if a trailer has strong indication of no translation,
used as part of negative training set
}
