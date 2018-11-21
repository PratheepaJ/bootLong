#'psTransform
#'
#' Variance-stablized values of otu table and residuals of glm.nb fit.
#'
#' @inheritParams plotSamplingSchedule
#'
#' @return A list of phyloseq objects. The first phyloseq element contains the arcsinh transformed of otu_table and the second phyloseq element contains the residuals of \code{\link[MASS]{glm.nb}} fit in otu table.
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

    geoMean <- function(x) {
        if(all(x == 0)){
            val <- 0
        }else{
            val <- exp(sum(log(x[x > 0]))/length(x))
        }
        return(val)
    }

    geom.mean.row <- apply(ot, 1, FUN = geoMean)

    sj <- estimateSizeFactorsForMatrix(ot, median, geoMeans = geom.mean.row)

    ot.trans <- t(asinh(t(ot)/sj))
    colnames(ot.trans) <- sample_names(ps)
    rownames(ot.trans) <- taxa_names(ps)

    ps.ot.asinh <- phyloseq(otu_table(ot.trans, taxa_are_rows = TRUE), sample_data(ps),tax_table(ps))

    # computing weights and residuals
    des <- as.formula(paste("otu", "~", paste(main_factor, collapse = "+"), "+", "offset(asinh(sj))"))

    ## computing weights
    samdf <- sample_data(ps) %>% data.frame
    des.v <- as.formula(paste("~", paste(main_factor, collapse = "+")))
    mm <- model.matrix(des.v, data = samdf)
    v <- asinhVoom(counts = ot, design = mm, sj = sj)
    weights.cal <- v$weights

    ## computing residuals
    residulasFitted <- function(ind.I, samdf.I, ot.I, sj.I, des.I, weights.I) {
        otu.T <- ot.I[ind.I, ] %>% as.numeric
        sj.T <- sj.I %>% as.numeric
        weight.T <- weights.I[ind.I, ] %>% as.numeric
        dff <- mutate(select_(samdf.I, main_factor), otu = otu.T, sj = sj.T, weightT = weight.T)

        dff[, main_factor] <- as.factor(dff[, main_factor])

        glmft <- suppressWarnings(MASS::glm.nb(des.I, data = dff, weights = weightT, method = "glm.fit", link = arcsinhLink()))
        res.residuals <- resid(glmft)
        rt <- list(res.residuals)
        names(rt) <- c("response_residuals")
        return(rt)
    }


    resi.fitted <- lapply(seq_len(ntaxa(ps)), function(x) {
        residulasFitted(ind.I = x, samdf.I = samdf, ot.I = ot, sj.I = sj,
            des.I = des, weights.I = weights.cal)
    })


    resi <- lapply(resi.fitted, "[[", 1)
    resi <- do.call("rbind", resi) %>% data.frame
    colnames(resi) <- sample_names(ps)
    rownames(resi) <- taxa_names(ps)

    ps.resid.asinh <- phyloseq(otu_table(resi, taxa_are_rows = T), sample_data(ps), tax_table(ps))

    rt <- list(ps.ot.asinh, ps.resid.asinh)
    names(rt) <- c("asinh_transformed_counts", "asinh_transformed_residulas")
    return(rt)
}

