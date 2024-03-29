% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassifierTrainingData.R
\name{makeTrainingData}
\alias{makeTrainingData}
\title{Get CDS and 3'UTR TrainingData of ribo seq features}
\usage{
makeTrainingData(
  tissues = "combined",
  features = c("countRFP", "disengagementScores", "entropyRFP", "floss", "fpkmRFP",
    "ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage", "startRegionCoverage",
    "startRegionRelative"),
  mode = "uORF",
  max.artificial.length,
  dataFolder = get("dataFolder", envir = .GlobalEnv),
  requiredActiveCds = 30
)
}
\arguments{
\item{tissues}{Tissue to train on, use "combined" if you want all in one,
first run of training it is a ORFik experiment.}

\item{features}{features to train model on, any of the features created
during ORFik::computeFeatures, default:
\code{c("countRFP", "disengagementScores", "entropyRFP", "floss",
"fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
"startRegionCoverage","startRegionRelative")}}

\item{mode}{character, default: "uORF". alternative "aCDS". Do you want to predict
on uORFs or artificial CDS. if "aCDS" will run twice once for whole length CDS and one for
truncated CDS to validate model works for short ORFs. "CDS" is option to predict on
whole CDS.}

\item{max.artificial.length}{integer, default: 100, only applies if mode = "aCDS",
so ignore this for most people,
when creating artificial ORFs from CDS, how large should maximum ORFs be,
this number is 1/6 of maximum size of ORFs (max size 600 if artificialLength is 100)
Will sample random size from 6 to that number, if max.artificial.length is
2, you can get artificial ORFs of size (6, 9 or 12) (6, + 6 + (3x1), 6 + (3x2))}

\item{requiredActiveCds}{numeric, default 30. How many CDSs are required to be
detected active. Size of minimum positive training set. Will abort if not
bigger than this number.}
}
\value{
invisible(NULL), saved to disc
}
\description{
Positive set is cds, negative is downstream region of 3' UTRs
}
