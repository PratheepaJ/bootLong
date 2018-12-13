params <-
list(lC1 = 2L, lC2 = 2L)

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 7.5, fig.height = 5.5)

## ------------------------------------------------------------------------
pkgs <- c("ggplot2","dplyr","tidyr",
          "phyloseq", "limma","ashr",
          "gridExtra","MASS",
          "geeM", "BiocParallel",
          "doParallel", "parallel","magrittr",
          "joineR", "DESeq2", "grid", "devtools")

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
#library(R.utils)
library(BiocParallel)
library(doParallel)
library(parallel)
library(magrittr)
library(joineR)
library(DESeq2)
library(grid)
library(bootLong)
#devtools::load_all(".")

## ------------------------------------------------------------------------
ncores = as.integer(Sys.getenv("SLURM_NTASKS"))
if(is.na(ncores)) ncores <- parallel::detectCores()
ncores

setting_name <- "Sim"

## ----read_arg------------------------------------------------------------
lC1 <- params$lC1
lC1
lC2 <- params$lC2
lC2

## ------------------------------------------------------------------------
ps <- psSim

## ------------------------------------------------------------------------
names(sample_data(ps))
table(sample_data(ps)$SubjectID, sample_data(ps)$Time)

## ------------------------------------------------------------------------
theme_set(theme_bw())
set.theme <- theme_update(panel.border = element_blank(),
                    panel.grid = element_line(size = .8),
                    axis.ticks = element_blank(),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 8),
                    axis.text = element_text(size = 8),
                    axis.title = element_text(size = 8),
                    strip.background = element_blank(),
                    strip.text = element_text(size = 8),
                    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5))


sample_data(ps)$Group <- as.factor(sample_data(ps)$Group)

p <- plotSamplingSchedule(ps, time_var = "Time", subjectID_var = "SubjectID", main_factor = "Group", theme_manual = set.theme)

p <- p + scale_color_discrete(name  ="Group",breaks=c("Case", "Control"),labels=c("Case", "Control"))
p
#ggsave("./Example_sampling_schedule.eps", plot = p, width = 8, height = 5.5)

## ------------------------------------------------------------------------
if (dim(otu_table(ps))[2] != nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps))
    }
threshold <- .1
keep_asv <- apply(otu_table(ps), 1, function(x){
    sum(x > 0)
    }) > threshold*nsamples(ps)

ps <- prune_taxa(keep_asv, ps)

### Better to save the ps after all filtering
# fileN <- paste0("ps", setting_name, ".rds")
# saveRDS(ps, fileN)

## ----common-legend-------------------------------------------------------
plot_common_legend <- function(p){
    ggplot_feature = ggplot_gtable(ggplot_build(p))
    le <- which(sapply(ggplot_feature$grobs, function(x) x$name) == "guide-box")
    l <- ggplot_feature$grobs[[le]]
    return(l)
}

## ----eval=FALSE----------------------------------------------------------
#  tax_table(ps)[, "Genus"] <- paste0(tax_table(ps)[, "Genus"],"_st_" ,seq(1, ntaxa(ps)))
#  

## ----pacf----------------------------------------------------------------
ps.tr <- psTransform(ps, main_factor = "Group") 

p.all <- longPACFMultiple(pstr  = ps.tr[[1]],
                             psres = ps.tr[[1]], 
                             main_factor = "Group", 
                             time_var = "Time", 
                             starttaxa = 1, 
                             endtaxa = 6,
                             taxlevel = "Genus")

#   Change the legend labels
p.all <- lapply(p.all, function(x){
    x + scale_fill_discrete(name  ="Group", breaks=c("Case", "Control"), labels=c("Case", "Control"))
    })

#   extract the common legend for all taxa
leg <- plot_common_legend(p.all[[1]])

plist <- lapply(p.all,function(x){
    x + theme(legend.position = "none")
    })

# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
                   
# ggsave("./Example_PAC.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs = plist, nrow=2, widths=c(3,3,3)),
             leg,
             ncol = 2,
             widths = c(10, 2))

## ----lag-plots, message=FALSE,warning=FALSE------------------------------
lags <- as.list(seq(1, 8))

p.lags <- lapply(lags,function(x){
    longLagPlot(ps.tr[[2]], 
        main_factor = "Group", 
        time_var = "Time", 
        taxon = 30, 
        x, 
        taxlevel = "Genus")})

#   Change the legend labels
p.lags <- lapply(p.lags, function(x){
    x + scale_color_discrete(name  = "Group", breaks = c("Case", "Control"), labels = c("Case", "Control"))
    })

leg <- plot_common_legend(p.lags[[1]])

plist <- lapply(p.lags,function(x){
    x+theme(legend.position = "none")
    })


# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# 
# ggsave("./Example_lag_plot.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs = plist, nrow = 2, widths = c(3, 3, 3, 3)), 
             leg, 
             ncol = 2, 
             widths = c(10, 2), top = textGrob("ASV 30", gp = gpar(fontsize = 8, font = 3)))


## ----corr, message=FALSE, warning=FALSE----------------------------------
p.all <- longCorreloMultiple(ps.tr[[1]],
                             ps.tr[[1]], 
                             main_factor = "Group", 
                             time_var = "Time", 
                             starttaxa = 1, 
                             endtaxa = 6,
                             taxlevel = "Genus")

#   Change the legend labels
p.all <- lapply(p.all, function(x){
    x + scale_fill_discrete(name  ="Group",breaks=c("Case", "Control"),labels=c("Case", "Control"))
    })

#   extract the common legend for all taxa
leg <- plot_common_legend(p.all[[1]])

plist <- lapply(p.all,function(x){
    x + theme(legend.position = "none")
    })

# p <- grid.arrange(arrangeGrob(grobs=plist,nrow=2,widths=c(3,3,3)),
#                   leg,
#                   ncol=2,
#                   widths=c(10,2))
# ggsave("./Example_correlogram.eps",plot=p,width = 8,height = 5.5)

grid.arrange(arrangeGrob(grobs = plist, nrow = 2, widths = c(3, 3, 3)),
             leg,
             ncol = 2,
             widths=c(10, 2))

## ----eval=FALSE----------------------------------------------------------
#  R <- 100
#  RR <- 50
#  main_factor <- "Group"
#  time_var <- "Time"
#  subjectID_var = "SubjectID"
#  sampleID_var = "SampleID"
#  lI <- 5
#  omega <- .6
#  
#  system.time(
#      mse.results <- bootLongSubsampling(ps = ps,
#          main_factor = main_factor,
#          time_var = time_var,
#          subjectID_var = subjectID_var,
#          sampleID_var = sampleID_var,
#          lI = lI,
#          R = R,
#          RR = RR,
#          omega = omega,
#          lC1 = lC1, lC2 = lC2,
#          ncores = ncores, psi.hat.lI = FALSE, psi.hat.lI.val = NULL, compStatParallel = FALSE)
#  )
#  
#  fileN <- paste0("./psi.hat.lI_", setting_name,".rds")
#  saveRDS(mse.results, fileN)

## ----eval=FALSE----------------------------------------------------------
#  # fileN <- paste0("ps", setting_name, ".rds")
#  # ps <- readRDS(fileN)
#  
#  fileN <- paste0("./psi.hat.lI_", setting_name,".rds")
#  psi.hat.lI.val <- readRDS(fileN)
#  
#  R <- 100
#  RR <- 50
#  main_factor <- "Group"
#  time_var <- "Time"
#  subjectID_var = "SubjectID"
#  sampleID_var = "SampleID"
#  lI <- 5
#  omega <- .6
#  
#  system.time(
#      mse.results <- bootLongSubsampling(ps = ps,
#          main_factor = main_factor,
#          time_var = time_var,
#          subjectID_var = subjectID_var,
#          sampleID_var = sampleID_var,
#          lI = lI,
#          R = R,
#          RR = RR,
#          omega = omega,
#          lC1 = lC1, lC2 = lC2,
#          ncores = ncores, psi.hat.lI = TRUE, psi.hat.lI.val = psi.hat.lI.val, compStatParallel = FALSE)
#  )
#  
#  fileN <- paste0("MSE_", setting_name, "_lC1", lC1,".rds")
#  saveRDS(mse.results, fileN)

## ----eval = FALSE--------------------------------------------------------
#  MSE.all <- list()
#  for(lC1 in 2:4){
#    filename <- paste0("MSE_",setting_name,"_lC1", lC1,".rds")
#    MSE.all[(lC1-1)] <- readRDS(filename)
#  }
#  
#  fileN <- paste0("MSE_",setting_name,".rds")
#  saveRDS(MSE.all, fileN)

## ----eval=FALSE----------------------------------------------------------
#  lC1 <- 2
#  omega <- .6
#  fileN <- paste0("MSE_",setting_name,".rds")
#  mse.results <- readRDS(fileN)
#  
#  blks <- length(mse.results)
#  mse <- list()
#  
#  for(i in 1:blks){
#      mse[[i]] <- mse.results[[i]]$MSE_i
#      }
#  
#  mse.avg <- lapply(mse,function(x){mean(x, na.rm=T)})
#  mse.sum <- lapply(mse,function(x){sum(x, na.rm=T)})
#  
#  mse <- unlist(mse.sum)
#  lblk <- seq(lC1, (length(mse)+1), by=1)
#  lblk.f <- as.factor(lblk)
#  dfp <- data.frame(mse = mse,lblk = lblk.f)
#  
#  l.omega <- lblk[mse==min(mse)]
#  l.omega
#  l.opt <- round((100/(omega*100))^(1/5)*l.omega, digits = 0)
#  l.opt
#  
#  p.mse <- ggplot(dfp, aes(x=lblk, y=mse, group=1))+
#      geom_point()+
#      geom_line()+
#      xlab("block size")+
#      ylab("Mean squared error")+
#      ggtitle(paste0("MSE for  ", setting_name, " ", omega*100,"%"," subsample")) +
#    theme(plot.title = element_text(hjust = 0.5))
#  
#  p.mse
#  
#  fileN <- paste0("MSE_", setting_name, ".eps")
#  ggsave(fileN, plot = p.mse, width = 8, height = 5.5)

## ----eval=FALSE----------------------------------------------------------
#  R <- 200
#  RR <- 50
#  main_factor <- "Group"
#  time_var <- "Time"
#  subjectID_var = "SubjectID"
#  sampleID_var = "SampleID"
#  
#  sample_data(ps)$Group <- factor(sample_data(ps)$Group)
#  sample_data(ps)$Group <- relevel(sample_data(ps)$Group, ref = "Control")
#  
#  system.time(
#      MBB.all <- bootLongMethod(ps,
#          main_factor = main_factor,
#          time_var = time_var,
#          sampleID_var = sampleID_var,
#          subjectID_var = subjectID_var,
#          b = l.opt,
#          R = R,
#          RR = RR, FDR = .05, compStatParallel = FALSE)
#      )
#  
#  fileN <- paste0("MBB_",setting_name,".rds")
#  saveRDS(MBB.all, fileN)

## ----eval=FALSE----------------------------------------------------------
#  fileN <- paste0("ps", setting_name, ".rds")
#  ps <- readRDS(fileN)
#  
#  FDR <- .05
#  taxalevel <- "Genus"
#  
#  fileN <- paste0("MBB_",setting_name,".rds")
#  boot.res.all <- readRDS(fileN) # list("summary", "beta.hat", "beta.hat.star", "beta.hat.star.star", "T.obs")
#  
#  summ <- boot.res.all$summary
#  
#  ### add strain number to interested taxonomy same as in EDA
#  tax_table(ps)[, taxalevel] <- paste0(tax_table(ps)[, taxalevel],"_st_",seq(1, ntaxa(ps)))
#  
#  ### if we change FDR, wee need to adjust the confidence interval
#  T.stars.and.obs <- boot.res.all$beta.hat.star
#  
#  lcl <- apply(T.stars.and.obs, 1, FUN=function(x){
#    quantile(x, probs=FDR/2, na.rm=TRUE)
#    })
#  ucl <- apply(T.stars.and.obs, 1, FUN=function(x){
#    quantile(x, probs=(1-FDR/2), na.rm=TRUE)
#    })
#  
#  summ$lcl <- lcl
#  summ$ucl <- ucl
#  
#  summ.filt.FDR <- dplyr::filter(summ, pvalue.adj <= FDR)
#  summ.filt.FDR <- dplyr::select(summ.filt.FDR, one_of("ASV","stat","pvalue.adj","lcl","ucl"))
#  summ.filt.FDR <- dplyr::filter(summ.filt.FDR, !is.na(pvalue.adj))
#  summ.filt.FDR <- dplyr::arrange(summ.filt.FDR, desc(stat))
#  
#  ### To save the results
#  summ.filt.FDR$stat <- round(summ.filt.FDR$stat, digits = 2)
#  summ.filt.FDR$pvalue.adj <- round(summ.filt.FDR$pvalue.adj, digits = 4)
#  summ.filt.FDR$pvalue.adj[which(summ.filt.FDR$pvalue.adj==0)] <- "<.0001"
#  summ.filt.FDR$lcl <- round(summ.filt.FDR$lcl, digits = 2)
#  summ.filt.FDR$ucl <- round(summ.filt.FDR$ucl, digits = 2)
#  
#  df.mbb <- data.frame(Taxa = summ.filt.FDR$ASV, betahat = summ.filt.FDR$stat, lcl = summ.filt.FDR$lcl, ucl = summ.filt.FDR$ucl, p.adj = summ.filt.FDR$pvalue.adj)
#  
#  df.mbb.write <- df.mbb
#  df.mbb.write$Taxa <- tax_table(ps)[as.character(df.mbb.write$Taxa), taxalevel] %>% as.character
#  
#  ####if you would like to avoid completely the ASV with strain number only.
#  
#  # is.it.st <- lapply(df.mbb.write$Taxa, function(x){(strsplit(x, "_") %>% unlist)[1]}) %>% unlist
#  # is.it.ind <- which(is.it.st == "st")
#  # df.mbb.write <- df.mbb.write[-is.it.ind, ]
#  
#  library(xtable)
#  fileN <- paste0("MBB_",setting_name,".tex")
#  print(xtable(df.mbb.write, type = "latex", digits = c(0,0,2,2,2,4)), file = fileN)

## ----eval=FALSE----------------------------------------------------------
#  ####if you would like to avoid completely the ASV with strain number only.
#  #df.mbb <- df.mbb[-is.it.ind, ]
#  #
#  df.mbb <- dplyr::filter(df.mbb, abs(beta) > 1)
#  ps.tr <- psTransform(ps, main_factor = "Group")[[1]]
#  
#  taxaName <-  tax_table(ps.tr)[, taxalevel] %>% as.character
#  ind <- which(taxa_names(ps) %in% as.character(df.mbb$Taxa))
#  tax.plot <- taxa_names(ps.tr)[ind]
#  
#  
#  otu.df <- otu_table(ps.tr) %>% t %>% data.frame
#  samdf <-  sample_data(ps.tr) %>% data.frame
#  otu.df <- bind_cols(otu.df, samdf)
#  colnames(otu.df)[1:ntaxa(ps.tr)] <- taxa_names(ps.tr)
#  
#  otu.df.l <- gather(otu.df, key=taxa, value=abund, 1:ntaxa(ps.tr))
#  otu.df.l$taxa <- otu.df.l$taxa %>% factor
#  otu.df.l$Time <- otu.df.l$Time %>% factor
#  
#  otu.df.l <- filter(otu.df.l,taxa %in% tax.plot)
#  otu.df.l$taxa <- droplevels(otu.df.l$taxa)
#  
#  txname <- taxaName[ind]
#  
#  match.tax.name <- data.frame(taxa = tax.plot, nam = txname)
#  rownames(match.tax.name) <- as.character(match.tax.name$taxa)
#  repl.tax.n <- match.tax.name[levels(otu.df.l$taxa),"nam"]
#  
#  otu.df.l$taxa <- factor(otu.df.l$taxa,labels  = as.character(repl.tax.n))
#  
#  otu.df.l$Time <- as.numeric(otu.df.l$Time)
#  
#  p <- ggplot(otu.df.l)+
#      geom_line(aes(x = Time, y = abund, linetype = Group, col = SubjectID, group = SubjectID))+
#      facet_wrap(~taxa,scales = "fixed")+
#      ylab("asinh")+
#    xlab("Gestational Week") +
#    ggtitle(paste0("Differentially Abundant Taxa (MBB) in ", setting_name)) +
#      theme(strip.text.x = element_text(size = 10, face="italic", margin = margin(.1, 0, .05, 0, "cm")), axis.text=element_text(size=8, face = "bold"), plot.title = element_text(hjust = 0.5)) + scale_color_discrete(guide = FALSE) +
#    scale_linetype_discrete(name  ="",breaks=c("Control", "Case"),labels=c("Control", "Case")) +
#    stat_summary(aes(x = Time,y = abund, linetype= Group),fun.y = mean, geom = "line")
#  
#  p
#  
#  fileN <- paste0("Inspect_diff_abun_", setting_name, ".eps")
#  ggsave(fileN, plot=p, width = 12, height = 8)

