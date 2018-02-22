#' bootLongIndices
#'
#' This function is called within \code{bootLongPhyloseq} for each subject's data frame.
#'
#'
#' @param x a data frame. It must have \code{Index} column that will be added by \code{bootLongPhyloseq}.
#' @param b a numeric. Block size.
#' @param time a character. Name of the time variable.
#' @param L a numeric. Number of overlapping blocks per subject that will be added by \code{bootLongPhyloseq}.
#' @param blks_first_index a vector of positive integer. This will be added by \code{bootLongPhyloseq} (same for all subjects.)
#'
#' @return A list of indices to include in the block bootstrap realization.
#' @export
bootLongIndices <- function(x,b,time,L,blks_first_index){

    #   if repeated samples are not ordered by 'time'
    if(!is.unsorted(x[,time])){# argument for `is.unsorted()` must be 1-dim
        x <- x
    }else{
        x <- arrange_(x,time)
    }

    #   number of repeated biological samples per subject
    num.of.rep.obs.x <- dim(x)[1]

    #   number of ovelapping blocks per subject
    num.of.blks <- num.of.rep.obs.x-b+1

    if(L>num.of.blks){
        # how many rows should be added or how many times repeated
        if((L-num.of.blks)%%num.of.rep.obs.x==0){
            howrep <- rep(1:num.of.rep.obs.x,times=(L-num.of.blks)/num.of.rep.obs.x)
        }else{
                howrep <- c(rep(1:num.of.rep.obs.x,times=(L-num.of.blks)/num.of.rep.obs.x),1:((L-num.of.blks)%%num.of.rep.obs.x))
        }

        expand.x <- bind_rows(x,x[howrep,])
        x <- expand.x
    }


    #  indices of biological samples selected in the bootstrap sample
    subject.sample.indices <- numeric(0)

    #   if there is less number of repeated observations
    #   than the block size (can happen for an unbalanced design),
    #   choose all the observations to the bootstrap sample
    if(L==0){#no need because we make up the series to account for the larger number of blocks
        subject.sample.indices <- x$Index

    }else{
        #   consecutive indices
        subject.sample.indices <- x$Index[unlist(lapply(blks_first_index,FUN=function(y){y:(y+b-1)}))]
        #   number of repeated observations in the bootstrap samples
        #   should be same as in the observed number of repeated observations
        subject.sample.indices <- subject.sample.indices[1:num.of.rep.obs.x]
    }
    return(list(subject.sample.indices))
}
