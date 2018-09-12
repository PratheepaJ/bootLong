#' bootLongPhyloseq
#'
#' Creates block bootstrap realization of a \code{phyloseq}.
#'
#' @inheritParams computeStat
#'
#' @return a \code{phyloseq} object. A block bootstrap realization of \code{ps}.
#' @export
bootLongPhyloseq = function(ps,
                            time_var,
                            subjectID_var,
                            b){

        if(dim(otu_table(ps))[1]==nsamples(ps)){
            otu_table(ps) = t(otu_table(ps,taxa_are_rows = T))
        }

        samdf = data.frame(sample_data(ps))

        samdf$Index = seq(1,nsamples(ps),by=1)

        if(!is.factor(samdf[,subjectID_var])){
            samdf[,subjectID_var] = as.factor(samdf[,subjectID_var])
            }

        otu.tab = otu_table(ps) %>% data.frame

        g = samdf[,subjectID_var]
        samdf.split.by.subjects = split(samdf,g)

        num.of.rep.obs = lapply(samdf.split.by.subjects,function(x){dim(x)[1]})

        num.of.blks = lapply(num.of.rep.obs,function(x){x-b+1})

        L = max(unlist(num.of.blks))
        blks.first.index = sample(1:L,L,replace = T)

        sampling.blks.within.subject.indices = lapply(samdf.split.by.subjects,FUN=function(q){
            bootLongIndices(x=q,
                            time_var=time_var,
                            b=b,
                            L=L,
                            blks_first_index =blks.first.index)
            })

        sampling.blks.within.subject.indices = lapply(sampling.blks.within.subject.indices,"[[",1)

        boot.sample.indices = as.numeric(unlist(sampling.blks.within.subject.indices))

        blk.boot.otu.tab = otu.tab[,boot.sample.indices]
        blk.boot.samdf = samdf[boot.sample.indices,]

        g2 = blk.boot.samdf[,subjectID_var]
        boot.samdf = split(blk.boot.samdf,g2)

        boot.samdf = lapply(boot.samdf,function(x){
            tim = seq(1,dim(x)[1])
            x[,time_var] = tim
            return(x)
        })

        blk.boot.samdf = do.call("rbind",boot.samdf)
        blk.boot.samdf = dplyr::select(blk.boot.samdf,-Index)
        colnames(blk.boot.otu.tab) = rownames(blk.boot.samdf)
        rownames(blk.boot.otu.tab) = taxa_names(ps)
        ps.boot = phyloseq::merge_phyloseq(X=otu_table(blk.boot.otu.tab,taxa_are_rows = T),sample_data(blk.boot.samdf),tax_table(ps))
        return(list(ps.boot))
}
