% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassifierHelpers.R
\name{getTissueFromFeatureTable}
\alias{getTissueFromFeatureTable}
\title{Get specific feature for specific tissues}
\usage{
getTissueFromFeatureTable(tableName, tissue)
}
\arguments{
\item{tableName}{name of table in sql database}

\item{tissue}{character, single tissue}
}
\value{
data.table of row mean values by tissue
}
\description{
Grouped by rowMeans
}
