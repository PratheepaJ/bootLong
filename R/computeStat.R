#' computeStat
#'
#' This function computes size factors for each samples using \code{estimateSizeFactors}.
#' Then, it transforms count data to inverse hyperbolic sine transformation.
#' Using the transformed counts, this function estimates the mean-variance relationship and
#' appropriate observation-level weights. These weights and the transformed data is used to fit a linear model
#' for each taxa. The empirical Bayes is applied to compute the shirnkage estimates for each factor.
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param main_factor vector of factor variable(s) in the sample data of ps.
#' @param time_var numeric, time_var variable at repeated observations.
#' @param b numeric, optimal block size to account for dependence within-subject.
#'
#' @return dataframe with the last column corresponds to taxa names
#' and other columns are for estimated coefficients for each factors using GEE, including Intercept.
#'
#' @import "edgeR"
#' @import "phyloseq"
#' @export
computeStat <- function(ps,
                        main_factor,
                        time_var,
                        subjectID_var,
                        b){

    if(dim(otu_table(ps))[2] != nsamples(ps)){
        otu_table(ps) = t(otu_table(ps))
    }

    ot = otu_table(ps) %>% data.frame %>% as.matrix

    if(!isTRUE(all(ot == floor(ot)))){
        stop("otu_table entries must be integer")
    }

    #   variance stabilization for negative binomial distribution accounting for library sizes.
    geo_mean = function(x){
        val = exp(mean(log(x[x>0])))
    }

    geom_mean_row = apply(ot, 1, FUN = geo_mean)

    median_ratios = function(x, geom_mean_row){
        rat = x/geom_mean_row
        median_rat <- median(rat)
        return(median_rat)
    }

    # library size normalization factor
    sj = apply(ot, 2, FUN = median_ratios, geom_mean_row = geom_mean_row)

    # ot_trans = t(asinh(t(ot)*sj))

    ##   compute weights

    # samdf = sample_data(ps) %>% data.frame
    des = as.formula(paste("otuT","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))
    #
    # response_fitted <- function(ind, samdf, ot, sj, des){
    #     otu = as.numeric(ot[ind,])
    #     sj = as.numeric(sj)
    #     dff = mutate(samdf, otu=otu, sj=sj)
    #     glmft = MASS::glm.nb(des, data = dff, method = "glm.fit", link = arcsinhLink())
    #     res_fitted = fitted(glmft, "response")
    #     return(res_fitted)
    # }
    #
    # fited = lapply(seq_len(ntaxa(ps)),function(x){
    #     response_fitted(x, samdf = samdf, ot = (ot+1), sj = sj, des = des)
    # })
    #
    # fited = data.frame(do.call("rbind",fited))
    # colnames(fited) = sample_names(ps)
    # rownames(fited) = taxa_names(ps)
    #
    # #   compute weight for each observation using the mean-variance relationship of asinh transformed observations accounting for library size
    # mean_asv = apply(ot_trans, 1, FUN = mean)
    # sd_asv = sqrt(apply(ot_trans, 1, FUN = sd))
    #
    # allzero = (rowSums(ot) == 0)
    #
    # if (any(allzero)) {
    #     mean_asv = mean_asv[!allzero]
    #     sd_asv = sd_asv[!allzero]
    # }
    #
    # l = lowess(mean_asv, sd_asv, f = span)
    #
    # f = approxfun(l, rule = 2)
    #
    # fited_trans = t(asinh(t(fited)/sj))
    #
    # w = 1/f(fited_trans)^4
    # dim(w) = dim(fited_trans)
    # weights.cal = w

    samdf = sample_data(ps) %>% data.frame
    des_v <- as.formula(paste("~", paste(main_factor, collapse="+")))
    mm <- model.matrix(des_v,data=samdf)
    v <- asinh_voom(counts=ot, design=mm, lib.size=sj)
    weights.cal = v$weights
    #   Estimate regression coefficients
    com_beta <- function(taxIndex,
                         sampleDf,
                         otuDf,
                         allSj,
                         weightDf,
                         desingGEE,
                         b,
                         subjectID_var,
                         time_var){

        otuT = as.numeric(otuDf[taxIndex,])
        allSj = as.numeric(allSj)
        weightT = as.numeric(weightDf[taxIndex,])
        dffT = cbind(sampleDf, otuT =otuT, allSj=allSj, weightT=weightT)

        dffT$idvar = as.numeric(as.factor(dffT[,subjectID_var]))
        idvar = "idvar"
        dffT = arrange_(dffT, idvar, time_var)

        glmft.tx = glm.nb(formula = desingGEE,
                          data = dffT,
                          weights = weightT,
                          method = "glm.fit",
                          link = arcsinhLink())

        #rese = as.vector(residuals(glmft.tx))
        rese = resid(glmft.tx,"response")
        rese = t(asinh(t(rese)*sj))

        dffT$res = rese
        dfsub = dffT
        if(!is.factor(dfsub[,time_var])){
            dfsub[,time_var] = as.factor(dfsub[,time_var])
        }

        meanT = data.frame(dfsub %>% group_by_(time_var) %>% summarise(meanr=mean(res)))

        if(!is.numeric(meanT[,time_var])){
            meanT[,time_var] <- as.numeric(as.character(meanT[,time_var]))
        }

        meanT = arrange_(meanT, time_var)

        meanr = ts(meanT$meanr, start = min(meanT[,time_var]), end = max(meanT[,time_var]), frequency = 1)

        acf.res = as.numeric(acf(meanr, plot = F, lag.max = length(meanr))$acf)
        acf.res[(b+1):length(acf.res)] = 0


        workCorr = matrix(nrow=length(acf.res), ncol = length(acf.res))
        for(rw in 1:length(acf.res)){
            for(nc in 1:length(acf.res)){
                workCorr[rw,nc] = acf.res[(abs(rw-nc)+1)]
            }
        }


        if(!is.numeric(dffT[,time_var])){
            dffT[,time_var] = as.numeric(as.character(dffT[,time_var]))
            }
        wavesTime = dffT[,time_var]
        idvarV = dffT[,"idvar"]
        theta = glmft.tx$theta

        init.beta = as.numeric(glmft.tx$coefficients)

        fit <-  tryCatch(geeM::geem(formula = desingGEE,
                                    id = idvarV,
                                    waves = wavesTime,
                                    data = dffT,
                                    family=arcsinhlstLink(),
                                    corstr = "fixed",
                                    weights = weightT,
                                    corr.mat = workCorr,
                                    init.beta = init.beta,
                                    nodummy=TRUE)$beta,
                         error = function(e){
                             t(glmft.tx$coefficients)
                             })

        return(fit)

    }

    df.beta.hat <- lapply(seq_len(ntaxa(ps)), function(x){
        com_beta(x,
                 sampleDf = samdf,
                 otuDf = (ot+1),
                 allSj = sj,
                 weightDf = weights.cal,
                 desingGEE = des,
                 b=b,
                 subjectID_var=subjectID_var,
                 time_var=time_var)
        })

    df.beta.hat = data.frame(do.call("rbind", df.beta.hat))



}


# otu_full <- as.matrix(round(otu_table(ps),digits = 0))+1
# samdf <- data.frame(sample_data(ps))
# if(!is.numeric(samdf[,time_var])){samdf[,time_var] <- as.numeric(samdf[,time_var])}
# if(!is.factor(samdf[,subjectID_var])){samdf[,subjectID_var] <- as.factor(samdf[,subjectID_var])}
#
# anno <- data.frame(tax_table(ps))
#
# dgeList <- edgeR::DGEList(counts=otu_full, genes=anno, samples = samdf)
#
# des <- as.formula(paste("~", paste(main_factor, collapse="+")))
# mm <- model.matrix(des,data=samdf)
#
# des2 <- as.formula(paste("otuT","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(allSj))"))
#
# geo.mean.row <- apply((otu_full+1),1,function(x){exp(sum(log(x))/length(x))})
# sj <- apply((otu_full+1),2,function(x){median(x/geo.mean.row)})
#
# v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj)
# weights.cal <- v$weights
#
# taxaLst <- as.list(seq(1,ntaxa(ps)))
#
# com_beta <- function(taxIndex, sampleDf, otuDf, allSj, weightDf, desingGEE, b, subjectID_var, time_var){
#     otuT <- as.numeric(otuDf[taxIndex,])
#     allSj <- as.numeric(allSj)
#     weightT <- as.numeric(weightDf[taxIndex,])
#     dffT <- cbind(sampleDf,otuT =otuT,allSj=allSj,weightT=weightT)
#
#     dffT$idvar <- as.numeric(as.factor(dffT[,subjectID_var]))
#     idvar <- "idvar"
#     dffT <- arrange_(dffT,idvar,time_var)
#
#     glmft.tx <- glm.nb(formula=desingGEE,data = dffT,weights = weightT,method = "glm.fit",link = arcsinhLink())
#
#     rese <- as.vector(residuals(glmft.tx))
#
#     dffT$res <- rese
#     dfsub <- dffT
#     if(!is.factor(dfsub[,time_var])){dfsub[,time_var] <- as.factor(dfsub[,time_var])}
#     meanT <- data.frame(dfsub %>% group_by_(time_var) %>% summarise(meanr=mean(res)))
#
#     if(!is.numeric(meanT[,time_var])){meanT[,time_var] <- as.numeric(as.character(meanT[,time_var]))}
#     meanT <- arrange_(meanT,time_var)
#
#     meanr <- ts(meanT$meanr,start = min(meanT[,time_var]),end = max(meanT[,time_var]),frequency = 1)
#
#     acf.res <- as.numeric(acf(meanr,plot = F,lag.max = length(meanr))$acf)
#     acf.res[(b+1):length(acf.res)] <- 0
#
#
#     workCorr <- matrix(nrow=length(acf.res),ncol = length(acf.res))
#     for(rw in 1:length(acf.res)){
#         for(nc in 1:length(acf.res)){
#             workCorr[rw,nc] <- acf.res[(abs(rw-nc)+1)]
#         }
#     }
#
#
#     if(!is.numeric(dffT[,time_var])){dffT[,time_var] <- as.numeric(as.character(dffT[,time_var]))}
#     wavesTime <- dffT[,time_var]
#     idvarV <- dffT[,"idvar"]
#     theta <- glmft.tx$theta
#
#     init.beta <- as.numeric(glmft.tx$coefficients)
#
#     fit <-  tryCatch(geeM::geem(formula=desingGEE,id=idvarV,waves = wavesTime,data = dffT,family=arcsinhlstLink(),corstr = "fixed",weights = weightT,corr.mat = workCorr,init.beta = init.beta,nodummy=TRUE)$beta,error=function(e){t(glmft.tx$coefficients)})
#
#     return(fit)
#
# }
#
# df.beta.hat <- lapply(taxaLst,function(x){com_beta(x, sampleDf=samdf, otuDf=otu_full, allSj=sj, weightDf=weights.cal, desingGEE=des2, b=b, subjectID_var=subjectID_var, time_var=time_var)})
#
# df.beta.hat <- data.frame(do.call("rbind",df.beta.hat))
#
# ASV <- taxa_names(ps)
# res <- bind_cols(df.beta.hat,ASV=ASV)
# rm(ps)
# return(res)
