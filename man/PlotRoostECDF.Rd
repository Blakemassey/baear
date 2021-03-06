% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/behavior.R
\name{PlotRoostECDF}
\alias{PlotRoostECDF}
\title{Plot roost behavior ECDF}
\usage{
PlotRoostECDF(df, pars)
}
\arguments{
\item{df}{Dataframe of location data.}

\item{pars}{List of simulation parameters with male/female, arrive/depart,
shape/scale.}
}
\value{
Facetted plot with roost data and fitted Weibull distributions
}
\description{
Plots roost empirical distribution funtion and fitted Weibull cumulative
  distribution function
}
\details{
Empirical distribution function extends to 15 min past the end of
  last time.
}
