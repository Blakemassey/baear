% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{ExtractRoostPars}
\alias{ExtractRoostPars}
\title{Extracts roost parameters}
\usage{
ExtractRoostPars(pars = sim_pars)
}
\arguments{
\item{pars}{List of simulation parameters with male/female, arrive/depart,
shape/scale}
}
\value{
Dataframe of roost parameters for male/female, arrive/depart,
  shape/scale
}
\description{
Returns df of roost parameters for male/female, arrive/depart, shape/scale
}
\details{
Used in PlotRoostECDF
}
