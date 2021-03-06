% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/behavior.R
\name{IfElseTimedateNA}
\alias{IfElseTimedateNA}
\title{Finds time value that is not NA}
\usage{
IfElseTimedateNA(df, col1, col2, result, default_tz, tz)
}
\arguments{
\item{df}{Dataframe of location data.}

\item{col1}{Column name, column checked to if is NA. If not NA, then col1 is
returned.}

\item{col2}{Column name, if col1 is NA, then col2 is returned.}

\item{result}{String, name of column to insert results. Default is "result".}

\item{default_tz}{String, timezone that timedate reverts to, based on OS
time.}

\item{tz}{string, timezone of original data (may be different from
default_tz).}
}
\value{
A dataframe with a new 'result' column.
}
\description{
Examines two time columns, if one is NA, the other time is returned.
}
\details{
Ensure that result is in the correct timezone. Used in
  RoostArrivalDeparture().
}
