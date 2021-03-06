% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longCorreloSingle.R
\name{longCorreloSingle}
\alias{longCorreloSingle}
\title{Plots the correlogram given the taxon index.}
\usage{
longCorreloSingle(ps, main_factor, time_var, taxon, taxlevel = "Species")
}
\arguments{
\item{ps}{(Required). A \code{\link{phyloseq-class}}.}

\item{main_factor}{Character string. The name of the variable to map main factor in the plot (factor vector).}

\item{time_var}{Character string. The name of the variable to map time in the plot (integer vector).}

\item{taxon}{Numeric, index of taxon to get the correlogram}

\item{taxlevel}{Character string. The taxonomy level for the plot title}
}
\value{
\code{ggplot2} object of correlogram for one taxon according to the taxon index.
}
\description{
Plots the correlogram given the taxon index.
}
