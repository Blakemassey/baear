% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/behavior.R
\name{ConvertNestIdToNum}
\alias{ConvertNestIdToNum}
\title{Converts nest alphanumeric id to nest numeric id}
\usage{
ConvertNestIdToNum(df)
}
\arguments{
\item{df}{input dataframe with "nest_id" column}
}
\value{
A dataframe with a numberic "nest_id_num" column
}
\description{
Converts alphanumeric "nest_id" column to a numeric "nest_id_num" column,
  useful for creating RasterLayers associated with nests
}
\details{
If a "nest_id_num" column doesn't exist, the function
  automatically creates one.
}
