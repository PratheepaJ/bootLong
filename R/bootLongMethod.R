#' bootLongMethod
#'
#' @param FDR A numeric. False discovery rate.
#' @inheritParams bootLongPhyloseq
#' @inheritParams bootLongSubsampling
#'
#' @return list of dataframe with ASV, observed stat, pvalues, adjusted pvalues, lcl, ucl, observed pivotal quantity; stat.obs; stat.star; stat.star.star; T.star_obs.
#' @export
bootLongMethod <- function(ps, main_factor, time_var, subjectID_var, sampleID_var, b, R, RR, FDR = 0.05) {

    doParallel::registerDoParallel(parallel::detectCores())
    BiocParallel::register(BiocParallel::DoparParam())

    if (dim(otu_table(ps))[1] == nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps, taxa_are_rows = T))
    }

    res.obs <- computeStat(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)

    stat.name <- colnames(res.obs)[2]

    boot.results <- list()

    boot.results <- bplapply(seq_len(R), FUN = function(i) {

        ps.boot <- bootLongPhyloseq(ps, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)

        ps.boot <- ps.boot[[1]]

        # df.boot <- computeStat(ps = ps.boot, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)
        df.boot <- ps.boot

            boot.results.bb <- lapply(seq_len(RR), FUN = function(j) {
                ps.boot.bb <- bootLongPhyloseq(ps.boot, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b)
                ps.boot.bb <- ps.boot.bb[[1]]
                # df.boot.bb <- computeStat(ps = ps.boot.bb, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, b = b)
                df.boot.bb <- ps.boot.bb
                # rm(ps.boot.bb)
                return(df.boot.bb)

            })

        # rm(ps.boot)

        return(list(df.boot, boot.results.bb))


    })

    # boot.results.all <- lapply(boot.results, "[[", 1)
    #
    # boot.results.bb <- lapply(boot.results, "[[", 2)
    #
    # rm(boot.results)
    #
    # stat.star <- do.call("cbind", lapply(boot.results.all, FUN = function(x) {
    #     x[, 2]
    # }))
    # stat.star <- as.data.frame(stat.star)
    #
    # sd.stat <- apply(stat.star, 1, FUN = sd, na.rm = FALSE)
    #
    # beta.obs <- res.obs[, 2]
    # shrink.beta.obs <- suppressMessages(ash(beta.obs, sebetahat = sd.stat, mixcompdist = "normal"))
    # shrink.beta.est <- shrink.beta.obs$result$PosteriorMean
    #
    # rm(boot.results.all)
    #
    # stat.obs <- data.frame(stat.obs = shrink.beta.est)
    #
    # T.obs <- data.frame(T.obs = stat.obs/sd.stat)
    #
    # stat.star.star <- lapply(boot.results.bb, function(x) {
    #     do.call("cbind", lapply(x, FUN = function(x) {
    #         x[, 2]
    #     }))
    #
    # })
    #
    # sd.stat.star <- do.call("cbind", lapply(stat.star.star, function(x) {
    #     data.frame(sd.stat.star = apply(x, 1, FUN = sd, na.rm = TRUE))
    # }))
    #
    # shrink.beta.boot <- lapply(seq_len(R), function(i) {
    #     suppressMessages(ash(stat.star[, i], sebetahat = sd.stat.star[, i], mixcompdist = "normal"))$result
    # })
    #
    # shrink.beta.boot.est <- lapply(shrink.beta.boot, function(ii) {
    #     ii$PosteriorMean
    # })
    #
    # stat.star <- do.call("cbind", shrink.beta.boot.est)
    #
    # T.num.star <- data.frame(apply(stat.star, 2, function(x) {
    #     x - stat.obs
    # }))
    #
    # T.star <- T.num.star/sd.stat.star
    #
    # T.star_obs <- dplyr::bind_cols(T.star, T.obs = T.obs[, 1])
    #
    # pvalue <- apply(T.star_obs, 1, function(x) {
    #     sum(abs(x[1:R]) >= abs(x[(R + 1)]))/R
    # })
    #
    # pvalue.adj <- data.frame(pvalue.adj = p.adjust(pvalue, method = "BH"))
    #
    # txnames <- dplyr::select(res.obs, ASV)
    # out <- data.frame(Taxa = txnames, stat = stat.obs[, 1], pvalue = pvalue, pvalue.adj = pvalue.adj)
    #
    # lcl <- apply(stat.star, 1, FUN = function(x) {
    #     quantile(x, probs = FDR/2, na.rm = TRUE)
    # })
    # ucl <- apply(stat.star, 1, FUN = function(x) {
    #     quantile(x, probs = (1 - FDR/2), na.rm = TRUE)
    # })
    #
    # out <- dplyr::bind_cols(out, lcl = lcl, ucl = ucl, T.obs = T.obs[, 1])
    #
    # rt <- list(out, stat.obs, stat.star, stat.star.star, T.star_obs)
    #
    # names(rt) <- c("summary", "beta.hat", "beta.hat.star", "beta.hat.star.star", "T.obs")
    rt <- boot.results
    return(rt)
}
