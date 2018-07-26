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

        des <- as.formula(paste("~", paste(factors, collapse="+")))
        mm <- model.matrix(des,data=samdf)

        des2 <- as.formula(paste("otuT","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(allSj))"))

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

                dffT$idvar <- as.numeric(as.factor(dffT[,SubjectID_n]))
                idvar <- "idvar"
                dffT <- arrange_(dffT,idvar,time)

                glmft.tx <- glm.nb(formula=desingGEE,data = dffT,weights = weightT,method = "glm.fit",link = arcsinhLink())

                rese <- as.vector(residuals(glmft.tx))

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



