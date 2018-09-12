#' bootLongIndices
#'
#' This function is called within \code{bootLongPhyloseq} for
#' each subject's data frame.
#'
#'
#' @param x a data frame. It must have \code{Index} column that will be added by \code{bootLongPhyloseq}.
#' @param L a numeric. Number of overlapping blocks per subject that will be added by \code{bootLongPhyloseq}.
#' @param blks_first_index a vector of positive integer. This will be added by \code{bootLongPhyloseq} (same for all subjects.)
#' @inheritParams computeStat
#'
#' @return A list of indices to include in the block bootstrap realization.
#' @export
bootLongIndices = function(x,
                            time_var,
                            b,
                            L,
                            blks_first_index){

        if(!is.numeric(x[,time_var])){x[,time_var] = as.numeric(x[,time_var])}

        if(!is.unsorted(x[,time_var])){
            x = x
        }else{
            x = arrange_(x,time_var)
        }


        num.of.rep.obs.x = dim(x)[1]

        num.of.blks = num.of.rep.obs.x-b+1

        if(L > num.of.blks){
            if((L-num.of.blks)%%num.of.rep.obs.x==0){
                howrep = rep(1:num.of.rep.obs.x, times=(L-num.of.blks)/num.of.rep.obs.x)
            }else{
                    howrep = c(rep(1:num.of.rep.obs.x, times=(L-num.of.blks)/num.of.rep.obs.x),1:((L-num.of.blks)%%num.of.rep.obs.x))
            }
            expand.x = bind_rows(x,x[howrep,])
            x = expand.x
        }

        subject.sample.indices = numeric(0)

        if(L==0){
            subject.sample.indices = x$Index
        }else{
            subject.sample.indices = x$Index[unlist(lapply(blks_first_index,FUN=function(y){y:(y+b-1)}))]
            subject.sample.indices = subject.sample.indices[1:num.of.rep.obs.x]
        }

        return(list(subject.sample.indices))
}
