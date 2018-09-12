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

    #   library size normalization factor.
    geo_mean = function(x){
        val = exp(mean(log(x[x>0])))
    }

    geom_mean_row = apply(ot, 1, FUN = geo_mean)

    median_ratios = function(x, geom_mean_row){
        rat = x/geom_mean_row
        median_rat <- median(rat)
        return(median_rat)
    }

    sj = apply(ot, 2, FUN = median_ratios, geom_mean_row = geom_mean_row)


    des = as.formula(paste("otuT","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))


    ##   compute weights
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
        rese = resid(glmft.tx, "response")
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

    ASV <- taxa_names(ps)
    res <- bind_cols(df.beta.hat,ASV=ASV)
    rm(ps)
    return(res)

}


