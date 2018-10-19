## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width=7.5, fig.height=5.5)

## ----install_packages, eval = FALSE--------------------------------------
#  pkgs <- c("ggplot2","dplyr","tidyr",
#            "phyloseq", "limma","ashr",
#            "gridExtra", "geepack","MASS",
#            "geeM", "R.utils", "BiocParallel",
#            "doParallel", "parallel","magrittr",
#            "joineR", "DESeq2")
#  
#  #   installed packages that were not installed already
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(setdiff(pkgs,installed.packages()), suppressUpdates = TRUE)
#  
#  devtools::install_github("PratheepaJ/bootLong")

## ----load_packages-------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(limma)
library(ashr)
library(gridExtra)
library(MASS)
library(geeM)
library(R.utils)
library(BiocParallel)
library(doParallel)
library(parallel)
library(magrittr)
library(joineR)
library(DESeq2)
#library(bootLong)

## ------------------------------------------------------------------------
ncores = as.integer(Sys.getenv("SLURM_NTASKS"))
if(is.na(ncores)) ncores <- parallel::detectCores()
ncores

## ------------------------------------------------------------------------
ps <- pssim

## ----eval=FALSE----------------------------------------------------------
#  names(sample_data(ps))
#  table(sample_data(ps)$SubjectID, sample_data(ps)$Time)

