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
#'
#' @return dataframe with the last column corresponds to taxa names
#' and other columns are for moderated coefficients for each factors, including Intercept.
#'
#' @import "edgeR"
#' @import "phyloseq"
#' @import "DESeq2"
#' @import "limma"
#' @import "ashr"
#' @export
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





computeStat <- function(ps,factors){
    ot <- as.matrix(otu_table(ps))
    anno <- data.frame(tax_table(ps))
    samdf <- data.frame(sample_data(ps))

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samdf)

    #   to estimate the size factors
    pse <- ps
    otu_table(pse) <- otu_table(ps)+1
    pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
    rm(pse)
    pDE <- DESeq2::estimateSizeFactors(pDE)
    sj <- DESeq2::sizeFactors(pDE)

    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)

blk <- apply(data.frame(samdf$SubjectID),1,function(x){
substr(x,start=2,stop = nchar(x))
})
blk <- factor(blk)
df.mm <- data.frame(mm)

geeglm(as.numeric(ot[1,]+1)~samdf$Preterm,weights = v$weights[1,],id=blk,family = poisson(link = log),waves = samdf$Time,offset = sj)


    v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)

    # #   estimate coefficients
    # #   we have voom transformed data
    # #   Use limma::lmFit to fit linear model using computed weights
    # #   fit linear model for each taxa
    # if correlation is given then, gls is used, default is weighted least squares
    fit <- limma::lmFit(v, design = mm,weights = v$weights)

    #fit <- limma::lmFit(v, design = mm,weights = v$weights,offset=sj,method = "robust")

    # #   Testing
    # fit4 <- lapply(as.list(1:dim(ot)[1]),function(x){
    #     #glm(as.numeric(ot[x,])~samdf$Preterm,offset=sj,weights = v$weights[x,],family = quasipoisson(link = asinhLink))
    #     geeglm(as.numeric(ot[x,])~samdf$Preterm,weights = v$weights[x,],id=samdf$Time,family = poisson)
    # })
    # fit4 <- lapply(fit4,function(x)coefficients(x))
    # fit4 <- do.call("rbind",fit4)
    ##########
    #   estimated coefficients
    df.beta.hat <- data.frame(fit$coefficients)
    df.sd.beta <- data.frame(fit$stdev.unscaled)
    dfreedom <- fit$df.residual[1]

    #   hist(data.frame(fit$coefficients)$PretermTRUE)
    #   empirical Bayes shrinkage of the estimators: ashr::ash()
    res <- df.beta.hat
    for(i in 1:dim(df.beta.hat)[2]){
        shrinkage_fit <- suppressMessages(ash(df.beta.hat[,i],sebetahat = df.sd.beta[,i],mixcompdist = "normal"))
        res[,i] <- shrinkage_fit$result$PosteriorMean
    }

    ASV <- taxa_names(ps)
    res <- bind_cols(res,ASV=ASV)
    res <- res[,c(2,3)]
    rm(ps)
    return(res)
}


