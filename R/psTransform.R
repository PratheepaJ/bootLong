#' psTransform
#' This function takes the \code{phyloseq} object and apply the arcsinh transformation.
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#'
#' @return \code{phyloseq} object with transformed \code{otu_table}
#' @export
#' @import "joineR"
psTransform <- function(ps,factors){
    ot <- as.matrix(otu_table(ps))
    anno <- data.frame(tax_table(ps))
    samdf <- data.frame(sample_data(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samdf)
    pse <- ps
    otu_table(pse) <- otu_table(ps)+1
    pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
    rm(pse)
    pDE <- DESeq2::estimateSizeFactors(pDE)
    sj <- DESeq2::sizeFactors(pDE)
    #   NOTE: if we want to use gene specific size factor, then
    #sij <- normalizationFactors(pDE) # matrix

    v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)

    transformed.ot <- data.frame(v$E)
    colnames(transformed.ot) <- sample_names(ps)

    pstr <- merge_phyloseq(otu_table(transformed.ot,taxa_are_rows = TRUE),sample_data(ps),tax_table(ps))
    return(pstr)
}
