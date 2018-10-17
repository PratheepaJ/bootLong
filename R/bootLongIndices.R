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

    x[,time_var] <- as.numeric(x[,time_var])

    if (!is.unsorted(x[, time_var])) {
        x <- x
    } else {
        x <- arrange_(x, time_var)
    }

    num_of_rep_obs_x <- dim(x)[1]

    num_of_blks <- num_of_rep_obs_x - b + 1

    if (num_of_blks > 0) {
        if ((L - num_of_blks)%%num_of_rep_obs_x == 0) {
            howrep <- rep(1:num_of_rep_obs_x, times = (L - num_of_blks)/num_of_rep_obs_x)
        } else {
            howrep <- c(rep(1:num_of_rep_obs_x, times = (L - num_of_blks)/num_of_rep_obs_x),
                1:((L - num_of_blks)%%num_of_rep_obs_x))
        }

        expand_x <- bind_rows(x, x[howrep, ])
        x <- expand_x
    }

    subject.sample.indices <- numeric(0)

    if (num_of_blks <= 0) {
        subject.sample.indices <- x$Index
    } else {
        subject.sample.indices <- x$Index[unlist(lapply(blks.first.index, FUN = function(y) {
            y:(y + b - 1)
        }))]
        subject.sample.indices <- subject.sample.indices[1:num_of_rep_obs_x]
    }

    return(subject.sample.indices)
}
