#' bboot_phyloseq
#'
#' Creates block bootstrap phyloseq object.
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param b Block lenght to account for dependence within-subject.
#' @param time Time variable at repeated measures.
#'
#' @return Block bootstrap phyloseq object as a list.
#' @export
bboot_phyloseq <- function(ps,b,time){

    #   add ``Index`` variable to ps sample data
    sample_data(ps)$Index <- seq(1,nsamples(ps),by=1)

    #   sample data of observed phyloseq
    samdf <- data.frame(sample_data(ps))

    #   otu table of observed phyloseq
    otu.tab <- data.frame(otu_table(ps))

    #   tax table of observed phyloseq
    tax.tab <- tax_table(ps)

    #   split observed data by subjects ---
    #   get a dataframe for each subject to do block bootstrap sampling
    samdf.split.by.subjects <- split(samdf,samdf$SubjectID)

    #   Number of repeated biological samples from a subject
    num.of.rep.obs <- lapply(samdf.split.by.subjects,function(x){dim(x)[1]})

    #   number of ovelapping blocks within a subject
    num.of.blks <- lapply(num.of.rep.obs,function(x){x-b+1})
    k <- max(unlist(num.of.blks))
    blks.first.index <- sample(1:k,k,replace = T)

    #   Sampling blocks within each subject
    sampling.blks.within.subject.indices <- lapply(samdf.split.by.subjects,FUN=function(q){bboot_indices(q,b,time,k,blks_first_index =blks.first.index)})

    #   Indices are first element of the list
    sampling.blks.within.subject.indices <- lapply(sampling.blks.within.subject.indices,"[[",1)

    # samples' indices in the block bootstrap sample
    boot.sample.indices <- unlist(sampling.blks.within.subject.indices)

    #   make a bootstrap sample of phyoloseq
    blk.boot.otu.tab <- otu.tab[,boot.sample.indices]
    blk.boot.samdf <- samdf[boot.sample.indices,]
    colnames(blk.boot.otu.tab) <- rownames(blk.boot.samdf)
    rownames(blk.boot.otu.tab) <- taxa_names(ps)
    #sum(rownames(blk.boot.otu.tab)==rownames(tax.tab))
    ps.boot <- merge_phyloseq(X=otu_table(blk.boot.otu.tab,taxa_are_rows = T),sample_data(blk.boot.samdf),tax_table(tax.tab))
    return(list(ps.boot))
}
