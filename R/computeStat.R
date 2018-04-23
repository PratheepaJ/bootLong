#' computeStat
#'
#' This function computes size factors for each samples using \code{estimateSizeFactors}.
#' Then, it transforms count data to inverse hyperbolic sine transformation.
#' Using the transformed counts, this function estimates the mean-variance relationship and
#' appropriate observation-level weights. These weights and the transformed data is used to fit a linear model
#' for each taxa. The empirical Bayes is applied to compute the shirnkage estimates for each factor.
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param time Time variable at repeated observations.
#' @param b numeric, optimal block size to account for dependence within-subject.
#'
#' @return dataframe with the last column corresponds to taxa names
#' and other columns are for estimated coefficients for each factors using GEE, including Intercept.
#'
#' @import "edgeR"
#' @import "phyloseq"
#' @import "DESeq2"
#' @import "limma"
#' @import "ashr"
#' @export
computeStat <- function(ps,factors,time,b){
    ot <- as.matrix(otu_table(ps))
    samdf <- data.frame(sample_data(ps))
    anno <- data.frame(tax_table(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samdf)

    #   to estimate the size factors
    # #   median-of-ratios
    #sizefa <- EBSeq::MedianNorm(ot)
    geo.mean.raw <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
    sj.com <- 0
    for(i in 1:dim(ot)[2]){
        sj.com[i] <- median(ot[,i]/geo.mean.raw)
    }

    pse <- ps
    otu_table(pse) <- otu_table(ps)+1
    pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
    #rm(pse)
    pDE <- DESeq2::estimateSizeFactors(pDE)
    sj <- DESeq2::sizeFactors(pDE)

    # View(data.frame(sj.com,sj))

    v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)

    df.beta.hat <- list()
    for(i in 1:ntaxa(ps)){
        # df.per.taxa <- data.frame(samdf,otu=as.numeric(t(ot[1,])),sj=sj,weights=v$weights[1,])
        df.per.taxa <- data.frame(samdf,otu=as.numeric(t(v$E[i,])),sj=sj,weights=v$weights[i,])

        #   need numeric values for Subject ID

        df.per.taxa$idvar <- as.numeric(as.factor(df.per.taxa$SubjectID))
        df.per.taxa <- arrange(df.per.taxa,SubjectID)

        #   compute residuals
        # gee1 <- geeglm(otu~Preterm, data=df.per.taxa,id=idvar,corstr = "independence",weights = df.per.taxa$weights,waves = df.per.taxa$Time,offset = df.per.taxa$sj,family = poisson(link = log))
        gee1 <- geeglm(otu~Preterm, data=df.per.taxa,id=idvar,corstr = "independence",weights = df.per.taxa$weights,waves = df.per.taxa$Time,offset = df.per.taxa$sj,family = gaussian)
        res <- df.per.taxa$otu-as.vector(residuals(gee1))

        #   compute sample correlation
        df.per.taxa$res <- res
        df.per.taxa$Time <- as.factor(df.per.taxa$Time)
        df.per.s <- df.per.taxa %>% group_by(Time) %>% summarise(meanr=mean(res))
        acf.res <- as.numeric(acf(df.per.s$meanr,plot = F)$acf)
        acf.res[(b+1):length(acf.res)] <- 0

        workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
        for(rw in 1:length(acf.res)){
            for(nc in 1:length(acf.res)){
                workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
            }
        }

        zcor <- fixed2Zcor(workCorr,id=df.per.taxa$idvar,waves = df.per.taxa$Time)
        #   gee with the working correlation
        # gee2 <- geeglm(otu~Preterm, data=df.per.taxa,id=idvar,corstr = "fixed",weights = v$weights[1,],waves = df.per.taxa$Time,offset = sj,family = poisson(link = log),zcor=zcor)
        gee2 <- geeglm(otu~Preterm, data=df.per.taxa,id=idvar,corstr = "fixed",weights = df.per.taxa$weights,waves = df.per.taxa$Time,offset = df.per.taxa$sj,family = gaussian,zcor=zcor)

        df.beta.hat[[i]] <- t(gee2$coefficients)
    }


    df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))

    ASV <- taxa_names(ps)
    res <- bind_cols(df.beta.hat,ASV=ASV)
    # res <- res[,c(2,3)]
    rm(ps)
    return(res)
}

### 1st attempt
# computeStat <- function(ps,factors){
#     ot <- as.matrix(otu_table(ps))
#     anno <- data.frame(tax_table(ps))
#     samdf <- data.frame(sample_data(ps))
#
#     #   setting up the model
#     des <- as.formula(paste("~", paste(factors, collapse="+")))
#     mm <- model.matrix(des,data=samdf)
#
#     #   to estimate the size factors
#     pse <- ps
#     otu_table(pse) <- otu_table(ps)+1
#     pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
#     rm(pse)
#     pDE <- DESeq2::estimateSizeFactors(pDE)
#     sj <- DESeq2::sizeFactors(pDE)
#
#     dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)
#
#     v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
#
#     df.mm <- data.frame(mm)
#     #   Estimate correlation
#     dupcor <- duplicateCorrelation(v,design = df.mm,block=samdf$SubjectID)
#
#     #   estimate coefficients
#     #   we have asinh transformed data
#     #   fit linear model using computed weights
#     #   fit linear model for each taxa
#
#     #fit <- limma::lmFit(v, design = mm,weights = v$weights,correlation = dupcor$atanh.correlations)
#
#     #   extract coefficients
#     df.beta.hat <- data.frame(fit$coefficients)
#     df.sd.beta <- data.frame(fit$stdev.unscaled)
#     dfreedom <- fit$df.residual[1]
#
#     #   hist(data.frame(fit$coefficients)$PretermTRUE)
#     #   empirical Bayes shrinkage of the estimators: ashr::ash()
#     res <- df.beta.hat
#     for(i in 1:dim(df.beta.hat)[2]){
#         shrinkage_fit <- suppressMessages(ash(df.beta.hat[,i],sebetahat = df.sd.beta[,i],mixcompdist = "normal"))
#         res[,i] <- shrinkage_fit$result$PosteriorMean
#     }
#
#     ASV <- taxa_names(ps)
#     res <- bind_cols(res,ASV=ASV)
#     res <- res[,c(2,3)]
#     rm(ps)
#     return(res)
# }
#



###### 2nd attempt

# computeStat <- function(ps,factors,time,b){
#     ot <- as.matrix(otu_table(ps))
#     anno <- data.frame(tax_table(ps))
#     samdf <- data.frame(sample_data(ps))
#
#     #   setting up the model
#     des <- as.formula(paste("~", paste(factors, collapse="+")))
#     mm <- model.matrix(des,data=samdf)
#
#     #   to estimate the size factors
#     pse <- ps
#     otu_table(pse) <- otu_table(ps)+1
#     pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
#     rm(pse)
#     pDE <- DESeq2::estimateSizeFactors(pDE)
#     sj <- DESeq2::sizeFactors(pDE)
#
#     dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)
#
#     blk <- apply(data.frame(samdf$SubjectID),1,function(x){
#     substr(x,start=2,stop = nchar(x))
#     })
#     blk <- factor(blk)
#     df.mm <- data.frame(mm)
#
#
#     df.taxa <- data.frame(sample_data(ps),otu=as.numeric(t(otu_table(ps)[1,])))
#     names(df.taxa)[names(df.taxa)==factors] <- "Group"
#     names(df.taxa)[names(df.taxa)==time] <- "Time"
#
#     df.taxa$Time <- as.factor(df.taxa$Time)
#     df.taxa <- df.taxa %>% group_by(Time) %>% summarise(meant=mean(otu))
#     acf.res <- as.numeric(acf(df.taxa$meant,plot=F)$acf)
#     acf.res[(b+1):length(acf.res)] <- 0
#
#     workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
#     for(rw in 1:length(acf.res)){
#         for(nc in 1:length(acf.res)){
#             workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
#         }
#     }
#
#     geese(as.numeric(ot[1,]+1)~samdf$Preterm,weights = v$weights[1,],id=blk,family = poisson(link = log),waves = samdf$Time,offset = sj)
#
#
#     v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
#
#     # #   estimate coefficients
#     # #   we have voom transformed data
#     # #   Use limma::lmFit to fit linear model using computed weights
#     # #   fit linear model for each taxa
#     # if correlation is given then, gls is used, default is weighted least squares
#     fit <- limma::lmFit(v, design = mm,weights = v$weights)
#
#     #fit <- limma::lmFit(v, design = mm,weights = v$weights,offset=sj,method = "robust")
#
#     # #   Testing
#     # fit4 <- lapply(as.list(1:dim(ot)[1]),function(x){
#     #     #glm(as.numeric(ot[x,])~samdf$Preterm,offset=sj,weights = v$weights[x,],family = quasipoisson(link = asinhLink))
#     #     geeglm(as.numeric(ot[x,])~samdf$Preterm,weights = v$weights[x,],id=samdf$Time,family = poisson)
#     # })
#     # fit4 <- lapply(fit4,function(x)coefficients(x))
#     # fit4 <- do.call("rbind",fit4)
#     ##########
#     #   estimated coefficients
#     df.beta.hat <- data.frame(fit$coefficients)
#     df.sd.beta <- data.frame(fit$stdev.unscaled)
#     dfreedom <- fit$df.residual[1]
#
#     #   hist(data.frame(fit$coefficients)$PretermTRUE)
#     #   empirical Bayes shrinkage of the estimators: ashr::ash()
#     res <- df.beta.hat
#     for(i in 1:dim(df.beta.hat)[2]){
#         shrinkage_fit <- suppressMessages(ash(df.beta.hat[,i],sebetahat = df.sd.beta[,i],mixcompdist = "normal"))
#         res[,i] <- shrinkage_fit$result$PosteriorMean
#     }
#
#     ASV <- taxa_names(ps)
#     res <- bind_cols(res,ASV=ASV)
#     res <- res[,c(2,3)]
#     rm(ps)
#     return(res)
# }

