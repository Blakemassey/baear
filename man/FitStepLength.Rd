% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{FitStepLength}
\alias{FitStepLength}
\title{Fit step length data}
\usage{
FitStepLength(data, location, scale, shape)
}
\arguments{
\item{df}{Dataframe with a "step_length" column}

\item{location}{Numeric, starting value of location, default is:
min(data$step_length)-1}

\item{scale}{Numeric, starting value of scale, default is 1.}

\item{shape}{Numeric, starting value of shape, default is 0.}

\item{pars}{}
}
\value{
List of (step_location, step_scale, step_shape)
}
\description{
Fits a generalized Pareto distribution to step-length data
}
