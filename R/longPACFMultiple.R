#' Plots the PACF given the multiple taxa indices.
#'
#' @param pstr (Required).  A \code{\link{phyloseq-class}}. The arcsinh transformation of the otu table.
#' @param starttaxa Numeric. The first index of taxa ordered by total reads
#' @param endtaxa Numeric. The last index of taxa ordered by total reads
#' @inheritParams longPACFSingle
#'
#' @return \code{\link{ggplot2-class}} object of PACF for multiple taxa.
#' @export
#'
longPACFMultiple <- function(pstr,
    main_factor,
    time_var,
    starttaxa = 1,
    endtaxa = 4,
    taxlevel = "Species",
    lag.max = 10) {

    taxa_order <- sort(taxa_sums(pstr), decreasing = T)
    ind <- which(taxa_names(pstr) %in% names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa, function(x) {
        longPACFSingle(ps = pstr,
            main_factor = main_factor,
            time_var = time_var,
            taxon = x,
            taxlevel = taxlevel,
            lag.max = lag.max)
    })
    return(p.all)
}

