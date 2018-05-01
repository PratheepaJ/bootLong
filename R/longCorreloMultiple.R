#' longCorreloMultiple
#'
#' This function plots the correlogram given the taxa indices.
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param time character, time variable at repeated observations.
#' @param starttaxa numeric, first index of taxa ordered by total reads
#' @param endtaxa numeric, the last index of taxa ordered by total reads
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot2} object of correlogram for multiple taxa.
#' @export
#'
longCorreloMultiple <- function(pstr,psres,factors,time,starttaxa=1,endtaxa=4,taxlevel="Species"){
    taxa_order <- sort(taxa_sums(pstr),decreasing = T)
    ind <- which(taxa_names(pstr)%in%names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa,function(x){longCorreloSingle(ps=psres,factors=factors,time=time,taxon=x,taxlevel = taxlevel)})
    return(p.all)
}

