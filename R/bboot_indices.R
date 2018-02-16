#' bboot_indices
#'
#' This function is called within \code{bboot_phyloseq} for each subject's data frame.
#'
#'
#' @param x a data frame. It must have \code{Index} column that will be added by \code{bboot_phyloseq}.
#' @param b a numeric. Block size.
#' @param time a character. Name of the time variable.
#' @param k a numeric. Number of overlapping blocks per subject that will be added by \code{bboot_phyloseq}.
#' @param blks_first_index a vector of positive integer. This will be added by \code{bboot_phyloseq} (same for all subjects.)
#'
#' @return A list of indices to include in the block bootstrap realization.
#' @export
bboot_indices <- function(x,b,time,k,blks_first_index){

    #   if samples are not ordered by 'time'
    if(!is.unsorted(x[,time])){# argument for `is.unsorted()` must be 1-dim
        x <- x
    }else{
        x <- arrange_(x,time)
    }

    #   Number of repeated biological samples from a subject
    num.of.rep.obs.x <- dim(x)[1]

    #   number of ovelapping blocks within a subject
    num.of.blks <- num.of.rep.obs.x-b+1
    # k <- num.of.blks

    if(k>num.of.blks){
        # how many rows should be added or how many times repeated
        # rows.x <- num.of.rep.obs.x
        if((k-num.of.blks)%%num.of.rep.obs.x==0){
            howrep <- rep(1:num.of.rep.obs.x,times=(k-num.of.blks)/num.of.rep.obs.x)
        }else{
                howrep <- c(rep(1:num.of.rep.obs.x,times=(k-num.of.blks)/num.of.rep.obs.x),1:((k-num.of.blks)%%num.of.rep.obs.x))
        }

        expand.x <- bind_rows(x,x[howrep,])
        x <- expand.x
        }


    #  indices of biological samples selected in the bootstrap sample
    subject.sample.indices <- numeric(0)

    #   if there is less number of repeated observations
    #   than the block length (can happen for unbalanced design),
    #   choose all the observations to the bootstrap sample
    if(k<=0){
        subject.sample.indices <- x$Index

    }else{
        # randomly select first index of each block
        #blks.first.index <- sample(1:k,k,replace = T)
        #   consecutive indices
        subject.sample.indices <- x$Index[unlist(lapply(blks_first_index,FUN=function(y){y:(y+b-1)}))]
        #   number of repeated observations in the bootstrap samples
        #   should be same as in the observed number of repeated observations
        subject.sample.indices <- subject.sample.indices[1:num.of.rep.obs.x]
    }
    return(list(subject.sample.indices))
}
