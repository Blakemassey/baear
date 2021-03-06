% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/behavior.R
\name{CreateHomeRangeKernelsParetoGamma}
\alias{CreateHomeRangeKernelsParetoGamma}
\title{Create home range kernels based on Pareto and Gamma functions}
\usage{
CreateHomeRangeKernelsParetoGamma(
  df_all,
  df_home = df_all,
  base = base,
  max_r,
  pareto_shape,
  pareto_scale,
  gamma_shape,
  gamma_scale,
  id = "nest_id",
  output_dir = getwd(),
  write_distance = TRUE,
  write_homerange = TRUE
)
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

\item{pareto_shape}{shape parameter of the Logistic function that governs
the home kernel}

\item{pareto_scale}{scale parameter of the Logistic function that governs
the home kernel}

\item{gamma_shape}{shape parameter of the Gamma function that governs the
conspecific avoidance kernel}

\item{gamma_scale}{scale parameter of the Gamma function that governs the
conspecific avoidance kernel}

\item{id}{column name of df_home that identifies the homerange. Default is
NULL, which sets the names to df_home row number.}

\item{output_dir}{directory for output files (distance, homerange)}

\item{write_distance}{logical, write distance raster to file. Default is
FALSE.}

\item{write_homerange}{logical, write home range raster to file. Default is
FALSE}
}
\value{
A list containing homerange kernel Rasters for all the df_home
  centroids
}
\description{
Creates RasterLayers of homerange kernels based on home centroids
}
