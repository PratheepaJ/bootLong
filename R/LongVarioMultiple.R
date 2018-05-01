#' longVarioMultiple
#'
#' This function plots the variogram for multiple taxa given the indeices.
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param starttaxa numeric, starting index of taxon: taxa will be ordered inside the function
#' @param endtaxa numeric, last index of taxon: taxa will be ordered inside the function
#' @param time character, time variable at repeated observations.
#' @param taxlevel character, taxonomy level to make the title
#'
#' @return \code{ggplot} object of variogram for starttaxa:endtaxa taxa
#' @export
#'
longVarioMultiple <- function(pstr,psres,factors,time,starttaxa=1,endtaxa=4,point=FALSE,taxlevel="Species"){
    # ps.tr <- ps_trans(ps,factors=factors)
    #   order tax by taxa sum
    #taxa_order <- names(sort(taxa_sums(ps),decreasing = TRUE))
    taxa_order <- sort(taxa_sums(pstr),decreasing = T)
    ind <- which(taxa_names(pstr)%in%names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa,function(x){longVarioSingle(ps=psres,factors=factors,time=time,taxon=x,point=point,taxlevel = taxlevel)})
    return(p.all)
}
