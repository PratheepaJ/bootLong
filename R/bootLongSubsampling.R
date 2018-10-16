#' bootLongSubsampling
#'
#' Finding optimal block size
#'
#' @param lI A numeric. The initial block size.
#' @param R A numeric. The number of block bootstrap realizations.
#' @param RR A numeric. The number of double block bootstrap realizations.
#' @param omega A numeric. The proportion between 0 and 1 of repeated observations for subsampling.
#' @param lC1 A numeric. (optional) any block size to start with. defualt is 1.
#' @param lC2 A numeric. (optional) any block size to end with. defualt is (lI-1).
#' @inheritParams bootLongMSEPsi
#' @inheritParams computeStat
#'
#' @return A list of ``MSE`` in calculating K with block sizes lC=1:(lI-1), ``Khat`` values with block sizes lC=1:(lI-1) and ``Khat.obs`` with lI
#' @export
#' @importFrom parallel mclapply

bootLongSubsampling <- function(ps, main_factor, time_var, subjectID_var, sampleID_var,
    lI, R, RR, omega = 0.6, lC1 = 1, lC2 = NULL, ncores) {

    if (dim(otu_table(ps))[1] == nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps, taxa_are_rows = T))
    }

    sam_ps <- sample_data(ps) %>% data.frame

    if (!is.numeric(sam_ps[, time_var])) {
        sam_ps[, time_var] <- as.numeric(sam_ps[, time_var])
    }

    g <- sam_ps[, subjectID_var]
    sam.ps.by.sub <- split(sam_ps, g)

    sam.ps.by.sub.mod <- lapply(sam.ps.by.sub, function(x){
        if (!is.unsorted(x[, time_var])) {
            x <- x
        } else {
            x <- arrange_(x, time_var)
        }
    })

    sam.ps.by.sub.mod <- do.call("rbind", sam.ps.by.sub.mod)
    rownames(sam.ps.by.sub.mod) <- sam.ps.by.sub.mod[,sampleID_var]

    sample_data(ps) <- sam.ps.by.sub.mod

    if (!is.factor(sam_ps[, subjectID_var])) {
        sam_ps[, subjectID_var] <- as.factor(sam_ps[, subjectID_var])
    }

    if (!is.factor(sam_ps[, sampleID_var])) {
        sam_ps[, sampleID_var] <- as.factor(sam_ps[, sampleID_var])
    }

    sample_data(ps) <- sam_ps

    qj <- table(sam_ps[, subjectID_var])

    Wj <- trunc(table(sam_ps[, subjectID_var]) * omega)

    if (lI > max(Wj)) {
        stop(paste("choose lI at most", max(Wj)))
    }

    if (is.null(lC2)) {
        lC2 <- (lI - 1)
    }

    lC <- seq(lC1, lC2, by = 1)

    psi.hat.lI <- bootLongPsi(ps = ps, main_factor = main_factor, time_var = time_var,
        subjectID_var = subjectID_var, b = lI, R = R, RR = RR, T.obs.full = NULL, ncores = ncores)


    Khat.obs <- psi.hat.lI[[1]]
    T.obs <- psi.hat.lI[[2]]

    blk_size_choice <- as.list(c(1:length(lC)))

    mseKhatKobs <- lapply(blk_size_choice, function(y) {
        bootLongMSEPsi(ps = ps, main_factor = main_factor, time_var = time_var,
            subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = y,
            R = R, RR = RR, qj = qj, Wj = Wj, Khat.obs = Khat.obs, T.obs.full = T.obs,
            ncores = ncores)
    })

    return(mseKhatKobs)
}
