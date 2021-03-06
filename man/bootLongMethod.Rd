% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootLongMethod.R
\name{bootLongMethod}
\alias{bootLongMethod}
\title{bootLongMethod}
\usage{
bootLongMethod(ps, main_factor, time_var, subjectID_var, sampleID_var, b,
  R, RR, FDR = 0.05, compStatParallel = FALSE)
}
\arguments{
\item{ps}{phyloseq object.}

\item{main_factor}{Character string. The name of the covariate variable.}

\item{time_var}{Character string. The name of the time variable.}

\item{subjectID_var}{Character string. The name of the subject ID variable.}

\item{sampleID_var}{Character string. The name of the sample ID variable.}

\item{b}{A numeric. The block size to account for the dependence within-subject.}

\item{R}{A numeric. The number of block bootstrap realizations.}

\item{RR}{A numeric. The number of double block bootstrap realizations.}

\item{FDR}{A numeric. False discovery rate.}

\item{compStatParallel}{A logical. True is used to use parallel in \code{computeStat} and disable parallel in \code{bootLongPsi} or \code{bootLongMethod}.}
}
\value{
list of dataframe with ASV, observed stat, pvalues, adjusted pvalues, lcl, ucl, observed pivotal quantity; stat.obs; stat.star; stat.star.star; T.star_obs.
}
\description{
bootLongMethod
}
