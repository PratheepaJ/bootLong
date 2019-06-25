#'arcsinh transformation of counts
#'
#' Variance-stabilizing transformation of otu table.
#'
#' @inheritParams plotSamplingSchedule
#'
#' @return A \code{\link{phyloseq-class}} with arcsinh transformed of otu_table.
#' @export
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

    sj <- estimateSizeFactorsForMatrix(ot, median,
        geoMeans = geom.mean.row)

    ot.trans <- t(asinh(t(ot)/sj))
    colnames(ot.trans) <- sample_names(ps)
    rownames(ot.trans) <- taxa_names(ps)

    ps.ot.asinh <- phyloseq(otu_table(ot.trans,
        taxa_are_rows = TRUE),
        sample_data(ps),
        tax_table(ps))

    rt <- ps.ot.asinh
    return(rt)
}

