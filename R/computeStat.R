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
#' @param time numeric, time variable at repeated observations.
#' @param b numeric, optimal block size to account for dependence within-subject.
#'
#' @return dataframe with the last column corresponds to taxa names
#' and other columns are for estimated coefficients for each factors using GEE, including Intercept.
#'
#' @import "edgeR"
#' @import "phyloseq"
#' @import "MASS"
#' @import "geeM"
#' @export
computeStat <- function(ps,factors,time,b,SubjectID_n="SubjectID"){
        
        otu_full <- as.matrix(round(otu_table(ps),digits = 0))+1
        samdf <- data.frame(sample_data(ps))
        if(!is.numeric(samdf[,time])){samdf[,time] <- as.numeric(samdf[,time])}
        if(!is.factor(samdf[,SubjectID_n])){samdf[,SubjectID_n] <- as.factor(samdf[,SubjectID_n])}
        
        anno <- data.frame(tax_table(ps))
        
        dgeList <- edgeR::DGEList(counts=otu_full, genes=anno, samples = samdf)

        #   setting up the model
        des <- as.formula(paste("~", paste(factors, collapse="+")))
        mm <- model.matrix(des,data=samdf)
        
        des2 <- as.formula(paste("otuT","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(allSj))"))
        
        #   estimate the size factors
        geo.mean.row <- apply((otu_full+1),1,function(x){exp(sum(log(x))/length(x))})
        sj <- apply((otu_full+1),2,function(x){median(x/geo.mean.row)})

        v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj)
        weights.cal <- v$weights

        taxaLst <- as.list(seq(1,ntaxa(ps)))

        com_beta <- function(taxIndex,sampleDf,otuDf,allSj,weightDf,desingGEE,b,SubjectID_n,time){
                otuT <- as.numeric(otuDf[taxIndex,])
                allSj <- as.numeric(allSj)
                weightT <- as.numeric(weightDf[taxIndex,])
                dffT <- cbind(sampleDf,otuT =otuT,allSj=allSj,weightT=weightT)
                
                #   need numeric values for Subject ID
                dffT$idvar <- as.numeric(as.factor(dffT[,SubjectID_n]))
                idvar <- "idvar"
                dffT <- arrange_(dffT,idvar,time)
                
                #       negative binomial family with arcsinh link
                glmft.tx <- glm.nb(formula=desingGEE,data = dffT,weights = weightT,method = "glm.fit",link = arcsinhLink())
                #   residuals
                rese <- as.vector(residuals(glmft.tx))
                
                #   compute sample correlation
                dffT$res <- rese
                dfsub <- dffT
                if(!is.factor(dfsub[,time])){dfsub[,time] <- as.factor(dfsub[,time])}
                meanT <- data.frame(dfsub %>% group_by_(time) %>% summarise(meanr=mean(res)))
                
                if(!is.numeric(meanT[,time])){meanT[,time] <- as.numeric(as.character(meanT[,time]))}
                meanT <- arrange_(meanT,time)
                
                meanr <- ts(meanT$meanr,start = min(meanT[,time]),end = max(meanT[,time]),frequency = 1)
                
                acf.res <- as.numeric(acf(meanr,plot = F,lag.max = length(meanr))$acf)
                acf.res[(b+1):length(acf.res)] <- 0
                
                
                workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
                for(rw in 1:length(acf.res)){
                        for(nc in 1:length(acf.res)){
                                workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
                        }
                }
                
                
                if(!is.numeric(dffT[,time])){dffT[,time] <- as.numeric(as.character(dffT[,time]))}
                wavesTime <- dffT[,time]
                idvarV <- dffT[,"idvar"]
                theta <- glmft.tx$theta
                 
                init.beta <- as.numeric(glmft.tx$coefficients)
                
                fit <-  tryCatch(geeM::geem(formula=desingGEE,id=idvarV,waves = wavesTime,data = dffT,family=arcsinhlstLink(),corstr = "fixed",weights = weightT,corr.mat = workCorr,init.beta = init.beta,nodummy=TRUE)$beta,error=function(e){t(glmft.tx$coefficients)})
                
                return(fit)

        }
        
        df.beta.hat <- lapply(taxaLst,function(x){com_beta(x,sampleDf=samdf,otuDf=otu_full,allSj=sj,weightDf=weights.cal,desingGEE=des2,b=b,SubjectID_n=SubjectID_n,time=time)})

        df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))

        ASV <- taxa_names(ps)
        res <- bind_cols(df.beta.hat,ASV=ASV)
        rm(ps)
        return(res)
}



### without correlation and used for oral data latest analysis
# computeStat <- function(ps,factors,time,b){
#     ot <- as.matrix(round(otu_table(ps),digits = 0))
#     samd <- data.frame(sample_data(ps))
#     anno <- data.frame(tax_table(ps))
#     dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)
#
#     #   setting up the model
#     des <- as.formula(paste("~", paste(factors, collapse="+")))
#     mm <- model.matrix(des,data=samd)
#     des2 <- as.formula(paste("otu","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))
#     #des2 <- as.formula(paste("otu","~", paste(factors, collapse="+")))
#     #   estimate the size factors
#     geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
#     sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})
#
#     v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
#     we <- v$weights
#
#     nt <- as.list(seq(1,ntaxa(ps)))
#     #names(samd)[names(samd)==factors] <- "Group"
#     com.beta <- function(ind,samd,ot,sj,we,des2,b){
#             otu <- as.numeric(ot[ind,])
#             sj <- as.numeric(sj)
#             we <- as.numeric(we[ind,])
#             dff <- samd
#             dff <- cbind(samd,otu=otu,sj=sj,we=we)
#             #       negative binomial family with arcsinh link
#             glmft.tx <- glm.nb(des2,data = dff,weights = we,method = "glm.fit",link = arcsinhLink())
#             return(t(glmft.tx$coefficients))
#     }
#
#     df.beta.hat <- lapply(nt,function(x){com.beta(x,samd,(ot+1),sj,we,des2,b)})
#
#     df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))
#
#     ASV <- taxa_names(ps)
#     res <- bind_cols(df.beta.hat,ASV=ASV)
#     rm(ps)
#     return(res)
# }


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
