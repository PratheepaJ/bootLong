#' Plots the variogram for multiple taxa given the indices.
#'
#' @inheritParams longCorreloMultiple
#' @inheritParams longVarioSingle
#'
#' @return \code{ggplot} object of variogram for for multiple taxa.
#' @export

longVarioMultiple <- function(pstr, psres, main_factor, time_var, subjectID_var,
    starttaxa = 1, endtaxa = 4, point = FALSE, taxlevel = "Species") {

    taxa_order <- sort(taxa_sums(pstr), decreasing = T)
    ind <- which(taxa_names(pstr) %in% names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa, function(x) {
        longVarioSingle(ps = psres, main_factor = main_factor, time_var = time_var,
            subjectID_var = subjectID_var, taxon = x, point = point, taxlevel = taxlevel)
    })
    return(p.all)
}
