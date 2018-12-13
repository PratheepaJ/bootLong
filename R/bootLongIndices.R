#' bootLongIndices
#'
#' This function is called within \code{bootLongPhyloseq} for
#' each subject's data frame.
#'
#'
#' @param x A dataframe. It must have \code{Index} column that will be added by \code{bootLongPhyloseq}.
#' @param L A numeric. Number of overlapping blocks per subject that will be added by \code{bootLongPhyloseq}.
#' @param blks.first.index A vector of positive integer. This will be added by \code{bootLongPhyloseq} (same for all subjects.)
#' @inheritParams computeStat
#'
#' @return A list of indices to include in the block bootstrap realization.
#' @export
bootLongIndices <- function(x, time_var, b, L, blks.first.index) {

    blks.first.index <- as.list(blks.first.index)

    if(!is.numeric(x[, time_var])){
        x[, time_var] <- as.numeric(x[, time_var])
        }

    if (!is.unsorted(x[, time_var])) {
        x <- x
    } else {
        x <- dplyr::arrange_(x, time_var)
    }

    num.of.rep.obs.x <- dim(x)[1]

    num.of.blks <- num.of.rep.obs.x - b + 1

    if ((L > num.of.blks)) {#unbalanced design of repeated observations
        if ((L - num.of.blks)%%num.of.rep.obs.x == 0) {
            howrep <- rep(1:num.of.rep.obs.x, times = (L - num.of.blks)/num.of.rep.obs.x)
        } else {
            howrep <- c(rep(1:num.of.rep.obs.x, times = (L - num.of.blks)/num.of.rep.obs.x),
                1:((L - num.of.blks)%%num.of.rep.obs.x))
        }

        expand_x <- dplyr::bind_rows(x, x[howrep, ])
        x <- expand_x
    }

    subject.sample.indices <- numeric(0)

    if (L==0) {
        subject.sample.indices <- x$Index
    } else {
        subject.sample.indices <- x$Index[unlist(lapply(blks.first.index, FUN = function(y) {
            y:(y + b - 1)
        }))]
        subject.sample.indices <- subject.sample.indices[1:num.of.rep.obs.x]
    }

    return(subject.sample.indices)
}
