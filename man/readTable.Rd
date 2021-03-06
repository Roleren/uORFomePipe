% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/databaseHelpers.R
\name{readTable}
\alias{readTable}
\title{Read a table from the database}
\usage{
readTable(
  tableName,
  asGR = FALSE,
  with.IDs = TRUE,
  uorfDB = get("uorfDB", envir = .GlobalEnv)
)
}
\arguments{
\item{tableName}{name of table in sql database}

\item{asGR}{convert to GRanges}

\item{with.IDs}{include ID column (remove to make calculations easier)}

\item{uorfDB}{the defined uORFome sql database, default:
get("uorfDB", envir = .GlobalEnv)}
}
\value{
the table as data.table or GRanges
}
\description{
Note: Because of SQL naming, any table you defined that has column names
starting with X, then a number directly afterwards, will lose the X.
}
