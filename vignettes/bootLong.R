## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width=7.5, fig.height=5.5)

## ----install_packages----------------------------------------------------
pkgs <- c("ggplot2","dplyr","tidyr",
          "phyloseq", "limma","ashr",
          "gridExtra", "geepack","MASS",
          "geeM", "R.utils", "BiocParallel",
          "doParallel", "parallel","magrittr",
          "joineR", "DESeq2")

#   installed packages that were not installed already
source("http://bioconductor.org/biocLite.R")
biocLite(setdiff(pkgs,installed.packages()), suppressUpdates = TRUE)

devtools::install_github("PratheepaJ/bootLong")

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
library(bootLong)

## ------------------------------------------------------------------------
ncores = as.integer(Sys.getenv("SLURM_NTASKS"))
if(is.na(ncores)) ncores <- parallel::detectCores()
ncores

## ------------------------------------------------------------------------
ps <- pssim

## ----eval=FALSE----------------------------------------------------------
#  names(sample_data(ps))
#  table(sample_data(ps)$SubjectID, sample_data(ps)$Time)

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

p <- plotSamplingSchedule(ps, time_var = "Time", subjectID_var = "SubjectID", main_factor = "Preterm", theme_manual = set.theme)

p <- p + scale_color_discrete(name  ="Group",breaks=c("FALSE", "TRUE"),labels=c("Term", "Preterm"))
p
#ggsave("./sampling_schedule.eps", plot=p, width = 8, height = 5.5)

## ----eval=FALSE----------------------------------------------------------
#  if (dim(otu_table(ps))[2] != nsamples(ps)) {
#          otu_table(ps) <- t(otu_table(ps))
#      }
#  threshold <- .1
#  keep_asv <- apply(otu_table(ps),1,function(x){sum(x>0)}) > threshold*nsamples(ps)
#  ps <- prune_taxa(keep_asv,ps)

## ----common-legend-------------------------------------------------------
plot_common_legend <- function(p){
    ggplot_feature = ggplot_gtable(ggplot_build(p))
    le <- which(sapply(ggplot_feature$grobs, function(x) x$name) == "guide-box")
    l <- ggplot_feature$grobs[[le]]
    return(l)
}

## ----pacf, message=FALSE, warning=FALSE----------------------------------
ps.tr <- psTransform(ps, 
                     main_factor = "Preterm") 

p.all <- longPACFMultiple(ps.tr[[1]],
                             ps.tr[[2]], 
                             main_factor = "Preterm", 
                             time_var = "Time", 
                             starttaxa = 1, 
                             endtaxa = 6,
                             taxlevel = "Genus")

#   Change the legend labels
p.all <- lapply(p.all, function(x){
    x + scale_fill_discrete(name  ="Group", breaks=c("FALSE", "TRUE"), labels = c("Term", "Preterm"))
    })

#   extract the common legend for all taxa
leg <- plot_common_legend(p.all[[1]])

plist <- lapply(p.all,function(x){
    x+theme(legend.position="none")
    })

# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# ggsave("./core_Sim.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
             leg,
             ncol=2,
             widths=c(10,2))

## ----corr, message=FALSE, warning=FALSE----------------------------------
p.all <- longCorreloMultiple(ps.tr[[1]],
                             ps.tr[[2]], 
                             main_factor = "Preterm", 
                             time_var = "Time", 
                             starttaxa = 1, 
                             endtaxa = 6,
                             taxlevel = "Genus")

#   Change the legend labels
p.all <- lapply(p.all, function(x){
    x + scale_fill_discrete(name  ="Group", breaks=c("FALSE", "TRUE"), labels = c("Term", "Preterm"))
    })

#   extract the common legend for all taxa
leg <- plot_common_legend(p.all[[1]])

plist <- lapply(p.all,function(x){
    x+theme(legend.position="none")
    })

# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# ggsave("./core_Sim.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
             leg,
             ncol=2,
             widths=c(10,2))

## ----lag-plots, message=FALSE,warning=FALSE------------------------------
lags <- as.list(seq(1,8))

p.lags <- lapply(lags,function(x){
    longLagPlot(ps.tr[[2]], 
        main_factor="Preterm", 
        time_var = "Time", 
        taxon = 1, 
        x, 
        taxlevel="Genus")})

#   Change the legend labels
p.lags <- lapply(p.lags, function(x){
    x+scale_color_discrete(name  ="Group", breaks=c("FALSE", "TRUE"), labels=c("Term", "Preterm"))
    })

leg <- plot_common_legend(p.lags[[1]])

plist <- lapply(p.lags,function(x){
    x+theme(legend.position="none")
    })


# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# 
# ggsave("./lag_sim.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3,3)), 
             leg, 
             ncol=2, 
             widths=c(10,2))


## ----message=FALSE, warning=FALSE----------------------------------------
p.all <- longVarioMultiple(ps.tr[[1]],
                           ps.tr[[2]],
                           main_factor = "Preterm",
                           time_var="Time",
                           subjectID_var = "SubjectID",
                           starttaxa = 1,
                           endtaxa = 6,
                           point = TRUE,
                           taxlevel = "Species")

#   Change the legend labels
p.all <- lapply(p.all, function(x){
    x + scale_color_discrete(name  ="Group", breaks=c("FALSE", "TRUE"), labels=c("Term", "Preterm "))
    })

leg <- plot_common_legend(p.all[[1]])

plist <- lapply(p.all,function(x){
    x + theme(legend.position="none")
    })

# p <- grid.arrange(arrangeGrob(grobs=plist, nrow=2, widths=c(3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# 
# ggsave("./vario_sim.eps", plot=p, width = 8, height = 5.5)

grid.arrange(arrangeGrob(grobs=plist, nrow=2, widths=c(3,3,3)),
             leg,
             ncol=2,
             widths=c(10,2))

