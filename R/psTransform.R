#'psTransform
#'
#' Variance-stablized values of otu table and residuals of glm.nb fit.
#'
#' @inheritParams plotSamplingSchedule
#'
#' @return A list of phyloseq objects. The first element is the asinh transformed of otu_table and the second element is the residuals of \code{\link[MASS]{glm.nb}} fit.
#' @export
#' @importFrom MASS glm.nb
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
psTransform <- function(ps, main_factor) {

    if (dim(otu_table(ps))[2] != nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps))
    }

    ot <- otu_table(ps) %>% data.frame %>% as.matrix

    if (!isTRUE(all(ot == floor(ot)))) {
        stop("otu_table entries must be integer")
    }

    # variance stabilization for negative binomial distribution accounting for
    # library sizes.
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

    ot_trans <- t(asinh(t(ot)/sj))
    colnames(ot_trans) <- sample_names(ps)
    rownames(ot_trans) <- taxa_names(ps)

    # computing weights and residuals
    des <- as.formula(paste("otu", "~", paste(main_factor, collapse = "+"), "+", "offset(asinh(sj))"))

    ## computing weights
    samdf <- sample_data(ps) %>% data.frame
    des_v <- as.formula(paste("~", paste(main_factor, collapse = "+")))
    mm <- model.matrix(des_v, data = samdf)
    v <- asinhVoom(counts = ot, design = mm, sj = sj)
    weights.cal <- v$weights

    ## computing residuals
    response_residulas_fitted <- function(ind, samdf, ot, sj, des, weights.cal) {
        otu <- as.numeric(ot[ind, ])
        sj <- as.numeric(sj)
        weightT <- as.numeric(weights.cal[ind, ])
        dff <- mutate(samdf, otu = otu, sj = sj, weightT = weightT)

        glmft <- tryCatch(MASS::glm.nb(des, data = dff, weights = weightT, method = "glm.fit", link = arcsinhLink()),
            error = function(e){
                dff$otu <- dff$otu + 1;glm(des, data = dff, weights = weightT, method = "glm.fit", family = poisson())#when count is very small
            })

        res_residuals <- resid(glmft)
        rt <- list(res_residuals)
        names(rt) <- c("response_residuals")
        return(rt)
    }


    resi_fitted <- lapply(seq_len(ntaxa(ps)), function(x) {
        response_residulas_fitted(x, samdf = samdf, ot = ot, sj = sj,
            des = des, weights.cal = weights.cal)
    })

    resi <- lapply(resi_fitted, "[[", 1)
    resi <- do.call("rbind", resi) %>% data.frame
    colnames(resi) <- sample_names(ps)
    rownames(resi) <- taxa_names(ps)

    ps_resid_asinh <- phyloseq(otu_table(resi, taxa_are_rows = T), sample_data(ps),
        tax_table(ps))

    ps_ot_asinh <- phyloseq(otu_table(ot_trans, taxa_are_rows = TRUE), sample_data(ps),
        tax_table(ps))

    rt <- list(ps_ot_asinh, ps_resid_asinh)
    names(rt) <- c("asinh_transformed_counts", "asinh_transformed_residulas")
    return(rt)
}

