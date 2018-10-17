#' bootLongPsi
#'
#' Compute \eqn{\psi} = two-sided probability with a given block size.
#'
#' @param R A numeric. The number of block bootstrap realizations.
#' @param RR A numeric. The number of double block bootstrap realizations.
#' @param T.obs.full A numeric vector. The numeric vector of observed statistic.
#' @param ncores A numeric. The number of cores to use in \code{\link[parallel]{mclapply}}.
#' @inheritParams computeStat
#' @inheritParams bootLongPhyloseq
#'
#' @return a list with first element ``K.val`` - two-sided significance probability.
#'         second element ``observed statistic``
#'
#' @export
bootLongPsi <- function(ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b, R, RR, T.obs.full = NULL, ncores){

    res.obs <- computeStat(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

    boot.results <- mclapply(seq_len(R), FUN = function(i){
        ps.boot <- bootLongPhyloseq(ps = ps, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)
        ps.boot <- ps.boot[[1]]

        df.boot <- computeStat(ps = ps.boot, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

            boot.results.bb <- lapply(seq_len(RR), FUN = function(j){
                ps.boot.bb <- bootLongPhyloseq(ps = ps.boot, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)
                ps.boot.bb <- ps.boot.bb[[1]]
                df.boot.bb <- computeStat(ps = ps.boot.bb, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

                return(df.boot.bb)
            })

            rt <- list(df.boot, boot.results.bb)
            return(rt)
    }, mc.cores = ncores)

    boot.results.all <- lapply(boot.results, "[[", 1)

    boot.results.bb <- lapply(boot.results, "[[", 2)

    rm(boot.results)

    stat.star <- lapply(boot.results.all, FUN = function(x){
        x[, 2]
    }) %>% data.frame

    colnames(stat.star) <- paste0("stat.star", seq(1, dim(stat.star)[2]))

    sd.stat <- apply(stat.star, 1, FUN = sd, na.rm = FALSE)

    beta.obs <- res.obs[, 2]

    shrink.beta.obs <- suppressMessages(ash(beta.obs, sebetahat = sd.stat,
        mixcompdist = "normal"))
    shrink.beta.est <- shrink.beta.obs$result$PosteriorMean

    rm(boot.results.all)

    stat.obs <- data.frame(stat.obs = shrink.beta.est)

    T.obs <- data.frame(T.obs = stat.obs/sd.stat)

    if (is.null(T.obs.full)) {
        T.obs <- T.obs
    } else {
        T.obs <- as.data.frame(T.obs.full)
    }

    colnames(T.obs) <- "T.obs"

    stat.star.star <- lapply(boot.results.bb, function(x) {
        lapply(x, FUN = function(y) {
            y[, 2]
        }) %>% data.frame
    })

    sd.stat.star <- lapply(stat.star.star, function(x) {
        data.frame(sd.stat.star = apply(x, 1, FUN = sd, na.rm = TRUE))
    }) %>% data.frame

    shrink.beta.boot <- lapply(seq_len(R), function(i) {
        suppressMessages(ash(stat.star[, i], sebetahat = sd.stat.star[, i], mixcompdist = "normal"))$result
    })

    shrink.beta.boot.est <- lapply(shrink.beta.boot, function(ii) {
        ii$PosteriorMean
    })

    stat.star <- shrink.beta.boot.est %>% data.frame

    T.num.star <- data.frame(apply(stat.star, 2, function(x) {
        x - stat.obs
    }))

    T.star <- T.num.star/sd.stat.star

    colnames(T.star) <- paste0("T.star",seq(1, dim(T.star)[2]))

    df.stat <- bind_cols(T.star, T.obs = T.obs[, 1])

    K.val <- apply(df.stat, 1, function(x) {
        sum(abs(x[1:R]) >= abs(x[(R + 1)]))/R
    })

    rt <- list(K.val, T.obs)
    names(rt) <- c("K.val", "T.obs")

    rm(list = c("stat.star", "stat.obs", "sd.stat", "stat.star.star", "sd.stat.star", "T.num.star", "df.stat"))

    return(rt)
}
