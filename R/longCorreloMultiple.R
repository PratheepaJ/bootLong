#' Plots the correlogram given the multiple taxa indices.
#'
#' @param pstr phyloseq object. The asinh transfomed of the otu_table.
#' @param psres phyloseq object. The asinh transformed of the response residuals.
#' @param starttaxa Numeric. The first index of taxa ordered by total reads
#' @param endtaxa Numeric. The last index of taxa ordered by total reads
#' @inheritParams longCorreloSingle
#'
#' @return \code{ggplot2} object of correlogram for multiple taxa.
#' @export
#'
longCorreloMultiple <- function(pstr, psres, main_factor, time_var, starttaxa = 1,
    endtaxa = 4, taxlevel = "Species") {

    taxa_order <- sort(taxa_sums(pstr), decreasing = T)
    ind <- which(taxa_names(pstr) %in% names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa, function(x) {
        longCorreloSingle(ps = psres, main_factor = main_factor, time_var = time_var,
            taxon = x, taxlevel = taxlevel)
    })
    return(p.all)
}

