## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Install-Pkgs,message=FALSE,warning=FALSE,results='hide'-------------
pkgs <- c("ggplot2","doParallel","foreach","reshape2","gridExtra","DESeq2","BiocParallel","ggfortify")
# git packages : "username/pkgname"
# source("http://bioconductor.org/biocLite.R")
# biocLite(setdiff(pkgs,installed.packages()))
#install_github("PratheepaJ/bbootLong")

## ----Load-pkgs,message=FALSE,warning=FALSE,results='hide'----------------
lpckges<-pkgs%in%(.packages())
if(any(!lpckges)){
    sapply(pkgs[!lpckges],require,character.only=TRUE)
}
#library(bbootLong)

## ----eval=FALSE----------------------------------------------------------
#  doParallel::registerDoParallel(parallel::detectCores())
#  BiocParallel::register(BiocParallel::DoparParam())

## ----eval=FALSE----------------------------------------------------------
#  ps <- pssim

## ----eval=FALSE----------------------------------------------------------
#  ps.tr <- ps_trans(ps,factors="Preterm")
#  p.all <- LMic_vario_multaxa(ps.tr, factors="Preterm",time="Time",1,6,taxlevel = "Species")
#  do.call("grid.arrange", c(p.all, ncol = 3))

## ----eval=FALSE----------------------------------------------------------
#  R <- 3
#  RR <- 2
#  factors <- "Preterm"
#  time <- "Time"
#  lI <- 5
#  omega <- .6
#  system.time(mse_results <- bboot_optblocksize(ps,R=R,RR=RR,factors=factors,time=time,lI=lI,omega=omega))
#  #saveRDS(mse_results,"./Results/mse_results.rds")

## ----eval=FALSE----------------------------------------------------------
#  #mse_results <- readRDS("./Results/mse_results.rds")
#  blks <- length(mse_results)
#  mse <- list();for(i in 1:blks){mse[[i]] <- mse_results[[i]]$MSE_i}
#  
#  mse.avg <- lapply(mse,function(x){mean(x,na.rm=T)})
#  mse.sum <- lapply(mse,function(x){sum(x,na.rm=T)})
#  
#  mse <- unlist(mse.sum)
#  lblk <- seq(1,length(mse),by=1)
#  lblk.f <- as.factor(lblk)
#  dfp <- data.frame(mse=mse,lblk=lblk.f)
#  
#  ggplot(dfp,aes(x=lblk,y=mse))+
#      geom_point()+
#      xlab("block size")+
#      ylab("Mean square error")
#  
#  l.M <- lblk[mse==min(mse)]
#  l.M
#  l.opt <- ceiling((100/(omega*100))^(1/5)*l.M)
#  l.opt

## ----eval=FALSE----------------------------------------------------------
#  R <- 3
#  RR <- 2
#  factors <- "Preterm"
#  time <- "Time"
#  
#  system.time(bboot_res <- bboot_LMic(ps,b=l.opt,R=R,RR=RR,factors=factors,time=time))
#  #saveRDS(bboot_res,"./Results/bboot_res.rds")

## ----eval=FALSE----------------------------------------------------------
#  #bboot_res <- readRDS("./Results/bboot_res.rds")
#  FDR <- .1
#  taxalevel <- "Species"
#  out <- bboot_res[[1]]
#  
#  T.star_obs <- bboot_res[[5]]
#  
#  #  Filter by FDR
#  out <- dplyr::filter(out,pvalue.adj <= FDR)
#  
#  #   Choose the taxalevel
#  taxaName <- as.character(tax_table(ps)[as.character(out$ASV),taxalevel])
#  
#  # #   if you want to append Species names to Genus taxalevel
#  # specName <- as.character(tax_table(ps)[as.character(out$ASV),"Species"])
#  # tog <- character()
#  #
#  # for(i in 1:length(taxaName)){
#  #     if(!is.na(specName[i])){
#  #     tog[i] <- paste(taxaName[i],specName[i])
#  #     }else{
#  #         tog[i] <- taxaName[i]
#  #     }
#  # }
#  # out$ASV <- tog
#  
#  out$ASV <- taxaName
#  #   remove ASV with NA
#  out <- dplyr::filter(out.res,!is.na(ASV))
#  out <- dplyr::arrange(out,desc(stat))
#  
#  
#  #   Write the results table in latex
#  library(xtable)
#  print(xtable(out.res, type = "latex",digits = 3), file = "./Results/MBB_Sim.tex")

