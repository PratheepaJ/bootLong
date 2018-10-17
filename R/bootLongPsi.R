#' bootLongPsi
#'
#' Compute \eqn{\psi} = two-sided probability with a given block size.
#'
#' @inheritParams computeStat
#' @inheritParams bootLongSubsampling
#' @param T.obs.full A numeric vector. The numeric vector of observed statistic.
#'
#' @return a list with first element ``K.val`` - two-sided significance probability.
#'         second element ``observed statistic``
#'
#' @export
bootLongPsi <- function(ps, main_factor, time_var, subjectID_var, sampleID_var, b, R, RR, T.obs.full = NULL, ncores) {

    res.obs <- computeStat(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

    boot.results <- list()

    boot.results <- mclapply(seq_len(R), FUN = function(i) {
        ps.boot <- bootLongPhyloseq(ps = ps, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)

        ps.boot <- ps.boot[[1]]

        df.boot <- computeStat(ps = ps.boot, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

        boot.results.bb <- lapply(seq_len(RR), FUN = function(j) {
            ps.boot.bb <- bootLongPhyloseq(ps.boot, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)

            ps.boot.bb <- ps.boot.bb[[1]]
            df.boot.bb <- computeStat(ps = ps.boot.bb, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)
            # rm(ps.boot.bb)
            return(df.boot.bb)
        })

        #rm(ps.boot)
        return(list(df.boot, boot.results.bb))
    }, mc.cores = ncores)

    boot.results.all <- lapply(boot.results, "[[", 1)

    boot.results.bb <- lapply(boot.results, "[[", 2)

    rm(boot.results)

    stat.star <- do.call("cbind", lapply(boot.results.all, FUN = function(x) {
        x[, 2]
    }))

    stat.star <- as.data.frame(stat.star)

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

    stat.star.star <- lapply(boot.results.bb, function(x) {
        do.call("cbind", lapply(x, FUN = function(x) {
            x[, 2]
        }))

    })

    sd.stat.star <- do.call("cbind", lapply(stat.star.star, function(x) {
        data.frame(sd.stat.star = apply(x, 1, FUN = sd, na.rm = TRUE))
    }))

    shrink.beta.boot <- lapply(seq_len(R), function(i) {
        suppressMessages(ash(stat.star[, i], sebetahat = sd.stat.star[, i],
            mixcompdist = "normal"))$result
    })

    shrink.beta.boot.est <- lapply(shrink.beta.boot, function(ii) {
        ii$PosteriorMean
    })

    stat.star <- do.call("cbind", shrink.beta.boot.est)

    T.num.star <- data.frame(apply(stat.star, 2, function(x) {
        x - stat.obs
    }))

    T.star <- T.num.star/sd.stat.star

    df.stat <- bind_cols(T.star, T.obs = T.obs[, 1])

    K.val <- apply(df.stat, 1, function(x) {
        sum(abs(x[1:R]) >= abs(x[(R + 1)]))/R
    })

    rt <- list(K.val, T.obs)

    rm(list = c("stat.star", "stat.obs", "sd.stat", "stat.star.star", "sd.stat.star",
        "T.num.star", "df.stat"))
    gc(reset = TRUE)

    return(rt)

}

