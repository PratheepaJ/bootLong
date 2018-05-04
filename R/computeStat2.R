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
#' @export
computeStat2 <- function(ps,factors,time,b){
    ot <- as.matrix(round(otu_table(ps),digits = 0))
    samd <- data.frame(sample_data(ps))
    names(samd)[names(samd)==time] <- "Time"
    anno <- data.frame(tax_table(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samd)
    des2 <- as.formula(paste("otu","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))
    #   estimate the size factors
    geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
    sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})

    v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
    we <- v$weights

    nt <- as.list(seq(1,ntaxa(ps)))
    #names(samd)[names(samd)==factors] <- "Group"

    # com.beta <- list()
    # for(ind in 1:ntaxa(ps)){
    #     otu <- as.numeric(ot[ind,])
    #     if(sum(otu==0)>10){otu <- otu+1}
    #     sj <- as.numeric(sj)
    #     we.ind <- as.numeric(we[ind,])
    #     dff <- samd
    #     dff <- cbind(samd,otu=otu,sj=sj,we=we.ind)
    #     #   need numeric values for Subject ID
    #     dff$idvar <- as.numeric(as.factor(dff$SubjectID))
    #     dff <- arrange(dff,SubjectID)
    # 
    #     #       negative binomial family with arcsinh link
    #     glmft.tx <- glm.nb(des2,data = dff,weights = we,method = "glm.fit",link = arcsinhLink())
    # 
    #     #   residuals
    #     rese <- as.vector(residuals(glmft.tx))
    # 
    #     #   compute sample correlation
    #     dff$res <- rese
    #     dfsub <- dff
    #     if(!is.factor(dfsub$Time)){dfsub$Time <- as.factor(dfsub$Time)}
    #     dfsub <- dfsub %>% group_by(Time) %>% summarise(meanr=mean(res))
    #     acf.res <- as.numeric(acf(dfsub$meanr,plot = F,lag.max = dim(dfsub)[1])$acf)
    #     acf.res[(b+1):length(acf.res)] <- 0
    # 
    #     workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
    #     for(rw in 1:length(acf.res)){
    #         for(nc in 1:length(acf.res)){
    #             workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
    #         }
    #     }
    # 
    #     if(!is.numeric(dff$Time)){dff$Time <- as.numeric(as.character(dff$Time))}
    #     wa <- dff$Time
    #     theta <- glmft.tx$theta
    #     init.beta <- as.numeric(glmft.tx$coefficients)
    #     fit <- geeM::geem(des2,id=idvar,waves = Time,data = dff,family=arcsinhlstLink(),weights = we.ind,corstr = "fixed",corr.mat = workCorr,nodummy=TRUE,init.beta = init.beta)
    # 
    #     com.beta[[ind]] <- fit$beta
    # 
    # }



    com.beta <- function(ind,samd,ot,sj,we,des2,b){
        otu <- as.numeric(ot[ind,])
        sj <- as.numeric(sj)
        we.ind <- as.numeric(we[ind,])
        dff <- samd
        dff <- cbind(samd,otu=otu,sj=sj,we=we.ind)
        #   need numeric values for Subject ID
        dff$idvar <- as.numeric(as.factor(dff$SubjectID))
        dff <- arrange(dff,SubjectID)

        #       negative binomial family with arcsinh link
        glmft.tx <- glm.nb(des2,data = dff,weights = we,method = "glm.fit",link = arcsinhLink())

        #   residuals
        rese <- as.vector(residuals(glmft.tx))

        #   compute sample correlation
        dff$res <- rese
        dfsub <- dff
        if(!is.factor(dfsub$Time)){dfsub$Time <- as.factor(dfsub$Time)}
        dfsub <- dfsub %>% group_by(Time) %>% summarise(meanr=mean(res))
        acf.res <- as.numeric(acf(dfsub$meanr,plot = F,lag.max = dim(dfsub)[1])$acf)
        acf.res[(b+1):length(acf.res)] <- 0

        workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
        for(rw in 1:length(acf.res)){
            for(nc in 1:length(acf.res)){
                workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
            }
        }

        if(!is.numeric(dff$Time)){dff$Time <- as.numeric(as.character(dff$Time))}
        wa <- dff$Time
        theta <- glmft.tx$theta
        init.beta <- as.numeric(glmft.tx$coefficients)
        fit <- geeM::geem(des2,id=idvar,waves = Time,data = dff,family=arcsinhlstLink(),weights = we,corstr = "fixed",corr.mat = workCorr,nodummy=TRUE,init.beta = init.beta,scale.fix = TRUE)

        return(fit$beta)

    }

    df.beta.hat <- lapply(nt,function(x){com.beta(x,samd,(ot+1),sj,we,des2,b)})

    df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))

    ASV <- taxa_names(ps)
    res <- bind_cols(df.beta.hat,ASV=ASV)
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


#### 3rd attempt
# computeStat2 <- function(ps,factors,time,b){
#     ot <- as.matrix(round(otu_table(ps),digits = 0))
#     samd <- data.frame(sample_data(ps))
#     names(samd)[names(samd)==time] <- "Time"
#     anno <- data.frame(tax_table(ps))
#     dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)
#
#     #   setting up the model
#     des <- as.formula(paste("~", paste(factors, collapse="+")))
#     mm <- model.matrix(des,data=samd)
#     des2 <- as.formula(paste("otu","~", paste(factors, collapse="+")))
#     #   estimate the size factors
#     geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
#     sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})
#
#     v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
#     v.E <- v$E
#     we <- v$weights
#
#     nt <- as.list(seq(1,ntaxa(ps)))
#     #names(samd)[names(samd)==factors] <- "Group"
#     com.beta <- function(ind,samd,ot,sj,we,des2,b){
#         otu <- as.numeric(ot[ind,])
#         sj <- as.numeric(sj)
#         we <- as.numeric(we[ind,])
#         dff <- samd
#         dff <- cbind(samd,otu=otu,sj=sj,we=we)
#         #   need numeric values for Subject ID
#         dff$idvar <- as.numeric(as.factor(dff$SubjectID))
#         dff <- arrange(dff,SubjectID)
#
#         #   compute residuals
#         gee1 <- geeglm(des2,data=dff,id=idvar,corstr = "independence",weights = we,offset = sj,family = poisson)
#         #gee1 <- glm(des2,data=dff,offset = sj,family = quasipoisson(link=arcsinhLink()))
#         rese <- as.vector(residuals(gee1))
#
#         #   compute sample correlation
#         dff$res <- rese
#         dfsub <- dff
#         dfsub$Time <- as.factor(dfsub$Time)
#         dfsub <- dfsub %>% group_by(Time) %>% summarise(meanr=mean(res))
#         acf.res <- as.numeric(acf(dfsub$meanr,plot = F)$acf)
#         acf.res[(b+1):length(acf.res)] <- 0
#
#         workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
#         for(rw in 1:length(acf.res)){
#             for(nc in 1:length(acf.res)){
#                 workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
#             }
#         }
#
#         id.zcor <- dff$idvar
#         wa <- dff$Time
#         zcor <- fixed2Zcor(workCorr,id=id.zcor,waves = wa)
#         #   gee with the working correlation
#         gee2 <- geeglm(des2, data=dff,id=idvar,corstr = "userdefined",weights = we,waves = Time,offset = sj,family = poisson,zcor=zcor)
#         #gee2 <- glm(des2,data=dff,offset = sj,family = quasipoisson(link=arcsinhLink()))
#         return(t(gee2$coefficients))
#         #return(gee2)
#
#     }
#
#     df.beta.hat <- lapply(nt,function(x){com.beta(x,samd,ot,sj,we,des2,5)})
#
#     # df.beta.hat <- list()
#     # for(i in 1:length(nt)){
#     #     df.beta.hat[[i]] <- com.beta(nt[[i]],v=v,sj=sj)
#     # }
#
#     df.beta.hat4 <- data.frame(do.call("rbind",df.beta.hat))
#
#     ASV <- taxa_names(ps)
#     res <- bind_cols(df.beta.hat,ASV=ASV)
#     rm(ps)
#     return(res)
# }
