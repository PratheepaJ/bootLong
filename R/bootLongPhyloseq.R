#' bootLongPhyloseq
#'
#' Creates block bootstrap realizations as \code{phyloseq} objects.
#'
#' @param ps a \code{phyloseq} object.
#' @param b a numeric. Block size.
#' @param time a character. Name of the time variable.
#'
#' @return a \code{phyloseq} object. This will be a block bootstrap realization of input \code{ps}.
#' @export
bootLongPhyloseq <- function(ps,b,time,SubjectID_n="SubjectID"){

    if(dim(otu_table(ps))[1]==nsamples(ps)){
        otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
    }

   
    #   sample data
    samdf <- data.frame(sample_data(ps))
    
    #   add ``Index``
    samdf$Index <- seq(1,nsamples(ps),by=1)
    
    if(!is.factor(samdf[,SubjectID_n])){samdf[,SubjectID_n] <- as.factor(samdf[,SubjectID_n])}

    #   otu table
    otu.tab <- data.frame(otu_table(ps))

    #   split data frame by subjects ---
    g <- samdf[,SubjectID_n]
    samdf.split.by.subjects <- split(samdf,g)

    #   number of repeated biological samples per subject
    num.of.rep.obs <- lapply(samdf.split.by.subjects,function(x){dim(x)[1]})

    #   number of ovelapping blocks per subject
    num.of.blks <- lapply(num.of.rep.obs,function(x){x-b+1})

    #   number of blocks per Subject if it is a balanced-design
    L <- max(unlist(num.of.blks))
    blks.first.index <- sample(1:L,L,replace = T)

    #   block bootstrap samples' indices
    sampling.blks.within.subject.indices <- lapply(samdf.split.by.subjects,FUN=function(q){bootLongIndices(x=q,b=b,time=time,L=L,blks_first_index =blks.first.index)})

    #   indices are first element of the list
    sampling.blks.within.subject.indices <- lapply(sampling.blks.within.subject.indices,"[[",1)

    #   all sample indices
    boot.sample.indices <- as.numeric(unlist(sampling.blks.within.subject.indices))

    #   block bootstrap realization
    blk.boot.otu.tab <- otu.tab[,boot.sample.indices]
    blk.boot.samdf <- samdf[boot.sample.indices,]
    #   split by subject
    g2 <- blk.boot.samdf[,SubjectID_n]
    boot.samdf <- split(blk.boot.samdf,g2)
    
    #   create new order of Time
    boot.samdf <- lapply(boot.samdf,function(x){
        tim <- seq(1,dim(x)[1])
        x[,time] <- tim
        return(x)
    })
    
    blk.boot.samdf <- do.call("rbind",boot.samdf)
    blk.boot.samdf <- dplyr::select(blk.boot.samdf,-Index)
    colnames(blk.boot.otu.tab) <- rownames(blk.boot.samdf)
    rownames(blk.boot.otu.tab) <- taxa_names(ps)
    ps.boot <- phyloseq::merge_phyloseq(X=otu_table(blk.boot.otu.tab,taxa_are_rows = T),sample_data(blk.boot.samdf),tax_table(ps))
    return(list(ps.boot))
}
