#'psTransform
#'
#' Variance-stablized values of otu table and residuals.
#'
#' @inheritParams plotSamplingSchedule
#'
#' @return A list of phyloseq objects. The first element is the asinh transfomed of the otu_table and the second element is the asinh transformed of the response residuals of the \code{\link[MASS]{glm.nb}} fit.
#' @export
#' @importFrom MASS glm.nb
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
        val <- exp(mean(log(x[x > 0])))
    }

    geom_mean_row <- apply(ot, 1, FUN = geo_mean)

    median_ratios <- function(x, geom_mean_row) {
        rat <- x/geom_mean_row
        median_rat <- median(rat [rat >0])
        return(median_rat)
    }

    # library size normalization factor
    sj <- apply(ot, 2, FUN = median_ratios, geom_mean_row = geom_mean_row)
    # matrix/vector or matrix*vector - rowwise division or rowwise
    # multiplication
    ot_trans <- t(asinh(t(ot) * sj))

    ## compute residuals and weights
    samdf <- sample_data(ps) %>% data.frame
    des <- as.formula(paste("otu", "~", paste(main_factor, collapse = "+"),
        "+", "offset(arcsinhLink()$linkfun(sj))"))

    ## compute weights
    des_v <- as.formula(paste("~", paste(main_factor, collapse = "+")))
    mm <- model.matrix(des_v, data = samdf)
    v <- asinhVoom(counts = ot, design = mm, sj = sj)
    weights.cal <- v$weights

    response_residulas_fitted <- function(ind, samdf, ot, sj, des, weights.cal) {
        otu <- as.numeric(ot[ind, ])
        sj <- as.numeric(sj)
        weightT <- as.numeric(weights.cal[ind, ])
        dff <- mutate(samdf, otu = otu, sj = sj, weightT = weightT)
        dff <- mutate(dff, weightT = ifelse(otu==0, 0, weightT))
        glmft <- MASS::glm.nb(des, data = dff, weights = weightT, method = "glm.fit",
            link = arcsinhLink())
        res_residuals <- resid(glmft, "response")
        rt <- list(res_residuals)
        names(rt) <- c("response_residuals")
        return(rt)
    }


    resi_fitted <- lapply(seq_len(ntaxa(ps)), function(x) {
        response_residulas_fitted(x, samdf = samdf, ot = ot, sj = sj,
            des = des, weights.cal = weights.cal)
    })

    resi <- lapply(resi_fitted, "[[", 1)
    resi <- data.frame(do.call("rbind", resi))
    colnames(resi) <- sample_names(ps)
    rownames(resi) <- taxa_names(ps)

    resi_trans <- t(asinh(t(resi) * sj))

    ps_resid_asinh <- phyloseq(otu_table(resi_trans, taxa_are_rows = T), sample_data(ps),
        tax_table(ps))

    ps_ot_asinh <- phyloseq(otu_table(ot_trans, taxa_are_rows = TRUE), sample_data(ps),
        tax_table(ps))

    rt <- list(ps_ot_asinh, ps_resid_asinh)
    names(rt) <- c("asinh_transformed_counts", "asinh_transformed_residulas")
    return(rt)
}

