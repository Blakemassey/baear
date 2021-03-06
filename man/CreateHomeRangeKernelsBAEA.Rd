% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/behavior.R
\name{CreateHomeRangeKernelsBAEA}
\alias{CreateHomeRangeKernelsBAEA}
\title{Create homerange kernels for the BAEA nests}
\usage{
CreateHomeRangeKernelsBAEA(df_all, df_home, base, max_r,
  home_inflection, home_scale, avoid_inflection, avoid_scale, output_dir,
  id, write_distance, write_homerange)
}
\arguments{
\item{df_all}{dataframe of all homerange centroids}

\item{df_home}{dataframe of homerange centroids to calculate homerange
kernels (these must be a subset of the df_all dataframe), Default is to
use 'df_all' dataframe}

\item{base}{base Raster that sets the projection, extent, and dimensions of
the study area}

\item{max_r}{maximum radius to calculate the homerange raster from each
df_home centroid}

\item{home_inflection}{inflection point of the Logistic function that
governs the home kernel}

\item{home_scale}{scale parameter of the Logistic function that governs the
scale parameter}

\item{avoid_inflection}{inflection point of the Logistic function that
governs the conspecific avoidance kernel}

\item{avoid_scale}{scale parameter of the Logistic function that governs the
conspecific avoidance kernel}

\item{min_prob}{numeric, minimum probablity threshold; all probabilities
below this value will be rounded to zero. Default is 0.}

\item{id}{column name of df_home that identifies the homerange. Default is
NULL, which sets the names to df_home row number.}

\item{name}{column name of df_home that identifies the name of the homerange
 and is used to write the file name of the .tif Default is 'id', which sets
the names to df_home row number if 'id' is NULL.}

\item{output_dir}{directory for output files (distance, homerange)}

\item{write_distance}{logical, write distance raster to file. Default is
FALSE.}

\item{write_homerange}{logical, write home range raster to file. Default is
FALSE.}
}
\value{
A list containing homerange kernel Rasters for all the df_home
  centroids.
}
\description{
Creates RasterLayers of homerange kernels based on home centroids
}
