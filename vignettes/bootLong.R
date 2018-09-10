## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----install_packages----------------------------------------------------
pkgs <- c("ggplot2","dplyr","tidyr",
          "phyloseq","DESeq2","edgeR",
          "limma","ashr","gridExtra",
          "geepack","MASS","geeM",
          "R.utils", "BiocParallel","doParallel",
          "parallel","magrittr")

#   installed packages that were not installed already
source("http://bioconductor.org/biocLite.R")
biocLite(setdiff(pkgs,installed.packages()), suppressUpdates = TRUE)

#devtools::install_github("PratheepaJ/bootLong")

## ----load_packages-------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(DESeq2)
library(edgeR)
library(limma)
library(ashr)
library(gridExtra)
library(geepack)
library(MASS)
library(geeM)
library(R.utils)
library(BiocParallel)
library(doParallel)
library(parallel)
library(magrittr)
devtools::load_all(".")

## ------------------------------------------------------------------------
ncores = as.integer(Sys.getenv("SLURM_NTASKS"))
if(is.na(ncores)) ncores <- parallel::detectCores()

## ------------------------------------------------------------------------
ps <- pssim

## ----eval=FALSE----------------------------------------------------------
#  psm <- merge_samples(ps, "SubjectID")
#  sample_data(psm)$SubjectID <- rownames(sample_data(psm))
#  sample_data(psm)$Preterm <- factor(sample_data(psm)$Preterm)
#  omat <- as(otu_table(psm), "matrix")
#  nsam <- table(sample_data(ps)$SubjectID)
#  omn <- omat/as.vector(nsam)
#  otu_table(psm) <- otu_table(omn, taxa_are_rows=FALSE)
#  if(dim(otu_table(psm))[2]!=nsamples(psm)){otu_table(psm)=t(otu_table(psm))}
#  

## ----eval=FALSE----------------------------------------------------------
#  table(sample_data(ps)$SubjectID,sample_data(ps)$Time)

## ------------------------------------------------------------------------
theme_set(theme_bw())
set.theme <- theme_update(panel.border = element_blank(),
                    panel.grid = element_line(size = .8),
                    axis.ticks = element_blank(),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 6),
                    axis.text = element_text(size = 8),
                    axis.title = element_text(size = 8),
                    strip.background = element_blank(),
                    strip.text = element_text(size = 8),
                    legend.key = element_blank())


sample_data(ps)$Preterm <- as.factor(sample_data(ps)$Preterm)
p <- plot_sampling_schedule(ps, time_var = "Time", subjectID_var = "SubjectID", main_factor = "Preterm", theme_manual = set.theme)
p <- p + scale_color_discrete(name  ="Group",breaks=c("FALSE", "TRUE"),labels=c("Term", "Preterm"))
#ggsave("./sampling_schedule.pdf",plot=p,width = 8,height = 5.5)

## ------------------------------------------------------------------------
threshold <- .01
keep_asv <- apply(otu_table(ps),1,function(x){sum(x>0)}) > threshold*nsamples(ps)
ps <- prune_taxa(keep_asv,ps)

## ------------------------------------------------------------------------
plot_common_legend <- function(p){
    ggplot_feature = ggplot_gtable(ggplot_build(p))
    le <- which(sapply(ggplot_feature$grobs, function(x) x$name) == "guide-box")
    l <- ggplot_feature$grobs[[le]]
    return(l)
}

