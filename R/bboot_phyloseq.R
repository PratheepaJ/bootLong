#' bboot_phyloseq
#'
#' Creates block bootstrap realizations as \code{phyloseq} objects.
#'
#' @param ps a \code{phyloseq} object.
#' @param b a numeric. Block size.
#' @param time a character. Name of the time variable.
#'
#' @return a \code{phyloseq} object. This will be a block bootstrap realization of input \code{ps}.
#' @export
bboot_phyloseq <- function(ps,b,time){

    #   add ``Index``
    sample_data(ps)$Index <- seq(1,nsamples(ps),by=1)

    #   sample data
    samdf <- data.frame(sample_data(ps))

    #   otu table
    otu.tab <- data.frame(otu_table(ps))

    #   tax table
    tax.tab <- tax_table(ps)

    #   split data frame by subjects ---
    samdf.split.by.subjects <- split(samdf,samdf$SubjectID)

    #   number of repeated biological samples per subject
    num.of.rep.obs <- lapply(samdf.split.by.subjects,function(x){dim(x)[1]})

    #   number of ovelapping blocks per subject
    num.of.blks <- lapply(num.of.rep.obs,function(x){x-b+1})
    L <- max(unlist(num.of.blks))
    blks.first.index <- sample(1:L,L,replace = T)

    #   block bootstrap samples' indices
    sampling.blks.within.subject.indices <- lapply(samdf.split.by.subjects,FUN=function(q){bboot_indices(q,b,time,L,blks_first_index =blks.first.index)})

    #   indices are first element of the list
    sampling.blks.within.subject.indices <- lapply(sampling.blks.within.subject.indices,"[[",1)

    #   all sample indices
    boot.sample.indices <- unlist(sampling.blks.within.subject.indices)

    #   block bootstrap realization
    blk.boot.otu.tab <- otu.tab[,boot.sample.indices]
    blk.boot.samdf <- samdf[boot.sample.indices,]
    colnames(blk.boot.otu.tab) <- rownames(blk.boot.samdf)
    rownames(blk.boot.otu.tab) <- taxa_names(ps)
    ps.boot <- merge_phyloseq(X=otu_table(blk.boot.otu.tab,taxa_are_rows = T),sample_data(blk.boot.samdf),tax_table(tax.tab))
    return(list(ps.boot))
}
