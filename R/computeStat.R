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
#' @import "geepack"
#' @export
computeStat <- function(ps,factors,time,b){
    ot <- as.matrix(round(otu_table(ps),digits = 0))
    samd <- data.frame(sample_data(ps))
    anno <- data.frame(tax_table(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samd)
    des2 <- as.formula(paste("otu","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))
    #des2 <- as.formula(paste("otu","~", paste(factors, collapse="+")))
    #   estimate the size factors
    geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
    sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})

    v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
    we <- v$weights

    nt <- as.list(seq(1,ntaxa(ps)))
    #names(samd)[names(samd)==factors] <- "Group"
    com.beta <- function(ind,samd,ot,sj,we,des2,b){
            otu <- as.numeric(ot[ind,])
            sj <- as.numeric(sj)
            we <- as.numeric(we[ind,])
            dff <- samd
            dff <- cbind(samd,otu=otu,sj=sj,we=we)
            #       negative binomial family with arcsinh link
            glmft.tx <- glm.nb(des2,data = dff,weights = we,method = "glm.fit",link = arcsinhLink())
            return(t(glmft.tx$coefficients))
    }

    df.beta.hat <- lapply(nt,function(x){com.beta(x,samd,(ot+1),sj,we,des2,b)})

    df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))

    ASV <- taxa_names(ps)
    res <- bind_cols(df.beta.hat,ASV=ASV)
    rm(ps)
    return(res)
}


# computeStat <- function(ps,factors,time,b){
#         ot <- as.matrix(round(otu_table(ps),digits = 0))
#         samd <- data.frame(sample_data(ps))
#         names(samd)[names(samd)==time] <- "Time"
#         anno <- data.frame(tax_table(ps))
#         dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)
#
#         #   setting up the model
#         des <- as.formula(paste("~", paste(factors, collapse="+")))
#         mm <- model.matrix(des,data=samd)
#         des2 <- as.formula(paste("otu","~", paste(factors, collapse="+")))
#         #   estimate the size factors
#         geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
#         sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})
#
#         v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
#         v.E <- v$E
#         we <- v$weights
#
#         nt <- as.list(seq(1,ntaxa(ps)))
#         #names(samd)[names(samd)==factors] <- "Group"
#         com.beta <- function(ind,samd,ot,sj,we,des2,b){
#                 otu <- as.numeric(ot[ind,])
#                 sj <- as.numeric(sj)
#                 we <- as.numeric(we[ind,])
#                 dff <- samd
#                 dff <- cbind(samd,otu=otu,sj=sj,we=we)
#                 #   need numeric values for Subject ID
#                 dff$idvar <- as.numeric(as.factor(dff$SubjectID))
#                 dff <- arrange(dff,SubjectID)
#
#                 #   compute residuals
#                 gee1 <- geeglm(des2,data=dff,id=idvar,corstr = "independence",weights = we,offset = sj,family = poisson)
#                 #gee1 <- glm(des2,data=dff,offset = sj,family = quasipoisson(link=arcsinhLink()))
#                 rese <- as.vector(residuals(gee1))
#
#                 #   compute sample correlation
#                 dff$res <- rese
#                 dfsub <- dff
#                 if(!is.factor(dfsub$Time)){dfsub$Time <- as.factor(dfsub$Time)}
#                 dfsub <- dfsub %>% group_by(Time) %>% summarise(meanr=mean(res))
#                 if(!is.numeric(dfsub$Time)){dfsub$Time <- as.numeric(as.character(dfsub$Time))}
#                 dfsub <- arrange(dfsub,Time)
#                 acf.res <- as.numeric(acf(dfsub$meanr,plot = F,lag.max = dim(dfsub)[1])$acf)
#                 acf.res[(b+1):length(acf.res)] <- 0
#
#                 workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
#                 for(rw in 1:length(acf.res)){
#                         for(nc in 1:length(acf.res)){
#                                 workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
#                         }
#                 }
#
#                 id.zcor <- dff$idvar
#                 if(!is.numeric(dff$Time)){dff$Time <- as.numeric(as.character(dff$Time))}
#                 wa <- dff$Time
#                 zcor <- fixed2Zcor(workCorr,id=id.zcor,waves = wa)
#                 #   gee with the working correlation
#                 gee2 <- geeglm(des2, data=dff,id=idvar,corstr = "fixed",weights = we,waves = Time,offset = sj,family = poisson,zcor=zcor)
#                 #gee2 <- glm(des2,data=dff,offset = sj,family = quasipoisson(link=arcsinhLink()))
#                 return(t(gee2$coefficients))
#                 #return(gee2)
#
#         }
#
#         df.beta.hat <- lapply(nt,function(x){com.beta(x,samd,ot,sj,we,des2,b)})
#
#         df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))
#
#         ASV <- taxa_names(ps)
#         res <- bind_cols(df.beta.hat,ASV=ASV)
#         rm(ps)
#         return(res)
# }
