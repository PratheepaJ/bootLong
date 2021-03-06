% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longPACFMultiple.R
\name{longPACFMultiple}
\alias{longPACFMultiple}
\title{Plots the PACF given the multiple taxa indices.}
\usage{
longPACFMultiple(pstr, main_factor, time_var, starttaxa = 1,
  endtaxa = 4, taxlevel = "Species", lag.max = 10)
}
\arguments{
\item{pstr}{(Required).  A \code{\link{phyloseq-class}}. The arcsinh transformation of the otu table.}

\item{main_factor}{Character string. The name of the variable to map main factor in the plot (factor vector).}

\item{time_var}{Character string. The name of the variable to map time in the plot (integer vector).}

\item{starttaxa}{Numeric. The first index of taxa ordered by total reads}

\item{endtaxa}{Numeric. The last index of taxa ordered by total reads}

\item{taxlevel}{Character string. The taxonomy level for the plot title.}

\item{lag.max}{Numeric. The maximum lags that PACF to be computed.}
}
\value{
\code{\link{ggplot2-class}} object of PACF for multiple taxa.
}
\description{
Plots the PACF given the multiple taxa indices.
}
