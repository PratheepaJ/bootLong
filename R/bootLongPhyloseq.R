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

        # if(dim(otu_table(ps))[1]==nsamples(ps)){
        #     otu_table(ps) = t(otu_table(ps,taxa_are_rows = T))
        # }


        sam_ps = sample_data(ps) %>% data.frame
        sam_ps %<>% mutate(Index = seq(1,nsamples(ps),by=1))
        #sam_ps$Index = seq(1,nsamples(ps),by=1)

        if(!is.factor(sam_ps[,subjectID_var])){
            sam_ps[,subjectID_var] = as.factor(sam_ps[,subjectID_var])
            }

        if(!is.numeric(sam_ps[,time_var])){
            sam_ps[,time_var] = as.numeric(sam_ps[,time_var])
        }

        ot = otu_table(ps) %>% data.frame

        g = sam_ps[,subjectID_var]
        sam_ps.split.by.subjects = split(sam_ps, g)

        num.of.rep.obs = lapply(sam_ps.split.by.subjects, function(x){
            dim(x)[1]
            })

        num.of.blks = lapply(num.of.rep.obs, function(x){
            x-b+1
            })

        L = num.of.blks %>% unlist %>% max
        blks.first.index = sample(1:L, L, replace = TRUE)

        sampling.blks.within.subject.indices = lapply(sam_ps.split.by.subjects,FUN=function(q){
            bootLongIndices(x = q,
                            time_var = time_var,
                            b = b,
                            L = L,
                            blks_first_index = blks.first.index)
            })

        sampling.blks.within.subject.indices = lapply(sampling.blks.within.subject.indices,"[[",1)

        boot.sample.indices = as.numeric(unlist(sampling.blks.within.subject.indices))

        blk.boot.ot = ot[,boot.sample.indices]
        blk.boot.sam_ps = sam_ps[boot.sample.indices,]

        g2 = blk.boot.sam_ps[,subjectID_var]
        boot.sam_ps = split(blk.boot.sam_ps,g2)

        boot.sam_ps = lapply(boot.sam_ps,function(x){
            tim = seq(1,dim(x)[1])
            x[,time_var] = tim
            return(x)
        })

        blk.boot.sam_ps = do.call("rbind",boot.sam_ps)
        blk.boot.sam_ps = dplyr::select(blk.boot.sam_ps,-Index)
        colnames(blk.boot.ot) = rownames(blk.boot.sam_ps)
        rownames(blk.boot.ot) = taxa_names(ps)
        ps.boot = phyloseq::merge_phyloseq(X=otu_table(blk.boot.ot,taxa_are_rows = T),sample_data(blk.boot.sam_ps),tax_table(ps))
        return(list(ps.boot))
}
