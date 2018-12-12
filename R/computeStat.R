#' computeStat
#'
#' Computes the regression coefficients.
#'
#' Computes the library normalization factors using the median-ratio method, computes the observation-level weight using the mean-variance relationship (modification to \code{\link[limma]{voom}}), computes the residuals using glm \code{\link[MASS]{glm.nb}}, computes the working correlation based on the given block size, computes the regression coefficient using the generalized estimating equation \code{\link[geeM]{geem}}.
#'
#' @param ps phyloseq object.
#' @param main_factor Character string. The name of the covariate variable.
#' @param time_var Character string. The name of the time variable.
#' @param subjectID_var Character string. The name of the subject ID variable.
#' @param b A numeric. The block size to account for the dependence within-subject.
#' @param compStatParallel A logical. True is used to use parallel in \code{computeStat} and disable parallel in \code{bootLongPsi} or \code{bootLongMethod}.
#' @return dataframe with the last column corresponds to taxa names
#' and other columns are for estimated coefficients for each covariate using GEE, including Intercept.
#' @export
computeStat <- function(ps, main_factor, time_var, subjectID_var, b, compStatParallel = FALSE) {

    # doParallel::registerDoParallel(parallel::detectCores())
    # BiocParallel::register(BiocParallel::DoparParam())

    if (dim(otu_table(ps))[2] != nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps))
    }

    ot <- otu_table(ps) %>% data.frame %>% as.matrix

    if (!isTRUE(all(ot == floor(ot)))) {
        # regression for negative binomial data so the
        stop("otu_table entries must be integer")
    }

    # library size normalization factor.
    geo_mean <- function(x) {
        if(all(x == 0)){
            val <- 0
        }else{
            val <- exp(sum(log(x[x > 0]))/length(x))
        }
        return(val)
    }

    geom_mean_row <- apply(ot, 1, FUN = geo_mean)


    sj <- estimateSizeFactorsForMatrix(ot, median, geoMeans = geom_mean_row)

    des <- as.formula(paste("otuT", "~", paste(main_factor, collapse = "+"), "+", "offset(asinh(allSjT))"))

    ## computing weights
    samdf <- sample_data(ps) %>% data.frame
    if(!is.numeric(samdf[, time_var])){
        samdf[, time_var] <- as.numeric(samdf[, time_var])
        }
    if(!is.factor(samdf[, subjectID_var])){
        samdf[, subjectID_var] <- as.factor(samdf[, subjectID_var])
        }

    des_v <- as.formula(paste("~", paste(main_factor, collapse = "+")))
    mm <- model.matrix(des_v, data = samdf)
    v <- asinhVoom(counts = ot, design = mm, sj = sj)
    weights.cal <- v$weights

    ot_trans <- t(asinh(t(ot)/sj))
    ## Estimating regression coefficients (using generalized estimating equation)

    com_beta <- function(taxIndex, sampleDf, otuDf, allSj, weightDf, desingGEE, b, subjectID_var, time_var, ot_trans) {

        otuT <- as.numeric(otuDf[taxIndex, ])
        allSj <- as.numeric(allSj)
        weightT <- as.numeric(weightDf[taxIndex, ])
        ot_transT <- as.numeric(ot_trans[taxIndex, ])
        dffT <- cbind(sampleDf, otuT = otuT, allSjT = allSj, weightT = weightT, ot_transT = ot_transT)
        dffT$idvar <- as.numeric(as.factor(dffT[, subjectID_var]))
        idvar <- "idvar"
        dffT <- arrange_(dffT, idvar, time_var)


        doGLMnbGEE <- function(){
            glmft.tx <- MASS::glm.nb(desingGEE, data = dffT, weights = weightT, method = "glm.fit", link = arcsinhLink())

            dfsub <- dffT

            dfsub[, subjectID_var] <- factor(dfsub[, subjectID_var], levels = unique(dfsub[, subjectID_var]))

            g <- dfsub[, subjectID_var]
            dfsub.sp <- split(dfsub, g)


            bootCorr <- lapply(dfsub.sp, function(x){
                bootLongWorkingCor(x$ot_transT, b)
            })

            workCorr <- bdiag(bootCorr) %>% as.matrix

            if (!is.numeric(dffT[, time_var])) {
                dffT[, time_var] <- as.numeric(as.character(dffT[, time_var]))
            }

            wavesTime <- dffT[, time_var]
            idvarV <- dffT[, "idvar"]

            theta <- glmft.tx$theta

            init.beta <- as.numeric(glmft.tx$coefficients)

            theta <- glmft.tx$theta

            LinkFun <- function(y){
                log(y + sqrt(y^2 + 1))
            }

            VarFun <- function(y){
                y * (1 + y/theta)
            }

            InvLink <- function(eta){
                pmax((0.5 * exp(-eta) * (exp(2 * eta) - 1)), .Machine$double.eps)
            }

            InvLinkDeriv <- function(eta){
                pmax((0.5 * (exp(eta) + exp(-eta))), .Machine$double.eps)
            }

            FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)

            geeM::geem(formula = desingGEE, id = idvarV, waves = wavesTime, data = dffT, family = FunList, corstr = "fixed", weights = weightT, corr.mat = workCorr, init.beta = init.beta, nodummy = TRUE)
        }

        doGLMPoiss <- function(){
            glm(desingGEE, data = dffT, weights = weightT, method = "glm.fit", family = poisson())
        }

        fit.m <- tryCatch(doGLMnbGEE(),
            error=function(e){
            doGLMPoiss()
        })


        if(class(fit.m) == "geem"){
            fit <- fit.m$beta
        }else{
            fit <- t(fit.m$coefficients)
        }

        return(fit)

    }

    ind <- as.list(c(1:ntaxa(ps)))

    if(compStatParallel){
        doParallel::registerDoParallel(parallel::detectCores())
        BiocParallel::register(BiocParallel::DoparParam())

        df.beta.hat <- BiocParallel::bplapply(ind, function(x) {
            rt <- com_beta(x, sampleDf = samdf, otuDf = ot, allSj = sj, weightDf = weights.cal, desingGEE = des, b = b, subjectID_var = subjectID_var, time_var = time_var, ot_trans = ot_trans)
            return(rt)
        })
    }else{
        df.beta.hat <- lapply(ind, function(x) {
            rt <- com_beta(x, sampleDf = samdf, otuDf = ot, allSj = sj, weightDf = weights.cal, desingGEE = des, b = b, subjectID_var = subjectID_var, time_var = time_var, ot_trans = ot_trans)
            return(rt)
        })
    }


    df.beta.hat <- data.frame(do.call("rbind", df.beta.hat))

    ASV <- taxa_names(ps)
    res <- bind_cols(df.beta.hat, ASV = ASV)
    rm(ps)
    return(res)

}


