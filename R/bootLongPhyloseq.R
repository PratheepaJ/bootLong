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

        samdf <- data.frame(sample_data(ps))

        samdf$Index <- seq(1,nsamples(ps),by=1)

        if(!is.factor(samdf[,SubjectID_n])){samdf[,SubjectID_n] <- as.factor(samdf[,SubjectID_n])}

        otu.tab <- data.frame(otu_table(ps))

        g <- samdf[,SubjectID_n]
        samdf.split.by.subjects <- split(samdf,g)

        num.of.rep.obs <- lapply(samdf.split.by.subjects,function(x){dim(x)[1]})

        num.of.blks <- lapply(num.of.rep.obs,function(x){x-b+1})

        L <- max(unlist(num.of.blks))
        blks.first.index <- sample(1:L,L,replace = T)

        sampling.blks.within.subject.indices <- lapply(samdf.split.by.subjects,FUN=function(q){bootLongIndices(x=q,b=b,time=time,L=L,blks_first_index =blks.first.index)})

        sampling.blks.within.subject.indices <- lapply(sampling.blks.within.subject.indices,"[[",1)

        boot.sample.indices <- as.numeric(unlist(sampling.blks.within.subject.indices))

        blk.boot.otu.tab <- otu.tab[,boot.sample.indices]
        blk.boot.samdf <- samdf[boot.sample.indices,]

        g2 <- blk.boot.samdf[,SubjectID_n]
        boot.samdf <- split(blk.boot.samdf,g2)

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
