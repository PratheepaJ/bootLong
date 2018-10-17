#' bootLongSubsampling
#'
#' Finding optimal block size.
#'
#' @param lI A numeric. The initial block size.
#' @param omega  A numeric. The proportion between 0 and 1 of repeated observations for subsampling.
#' @param lC1 A numeric. (optional) any block size to start with. Defualt is 1.
#' @param lC2 A numeric. (optional) any block size to end with. defualt is (lI-1).
#' @inheritParams bootLongMSEPsi
#'
#' @return A list of ``MSE`` in calculating K with block sizes lC=1:(lI-1), ``Khat`` values with block sizes lC=1:(lI-1) and ``Khat.obs`` with lI
#' @export
#' @importFrom parallel mclapply
bootLongSubsampling <- function(ps, main_factor, time_var, subjectID_var, sampleID_var, lI, R, RR, omega = .6, lC1 = 1, lC2 = NULL, ncroes){

    if (dim(otu_table(ps))[1] == nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps, taxa_are_rows = T))
    }

    sam.ps <- sample_data(ps) %>% data.frame

    # if (!is.numeric(sam.ps[, time_var])) {
    #     sam.ps[, time_var] <- as.numeric(sam.ps[, time_var])
    # }

    sam.ps[, subjectID_var] <- factor(sam.ps[, subjectID_var], levels = unique(sam.ps[, subjectID_var] ))

    # if (!is.factor(sam.ps[, sampleID_var])) {
    #     sam.ps[, sampleID_var] <- as.factor(sam.ps[, sampleID_var])
    # }
    #
    # g <- sam.ps[, subjectID_var]
    # sam.ps.by.sub <- split(sam.ps, g)
    #
    # sam.ps.by.sub.mod <- lapply(sam.ps.by.sub, function(x){
    #     if (!is.unsorted(x[, time_var])) {
    #         x <- x
    #     } else {
    #         x <- arrange_(x, time_var)
    #     }
    # })
    #
    # sam.ps.by.sub.mod <- do.call("rbind", sam.ps.by.sub.mod)
    # rownames(sam.ps.by.sub.mod) <- sam.ps.by.sub.mod[,sampleID_var] %>% as.character
    #
    # ps <- merge_phyloseq(otu_table(ps, taxa_are_rows = TRUE), sample_data(sam.ps.by.sub.mod), tax_table(ps))
    qj <- table(sam.ps[, subjectID_var])

    Wj <- trunc(table(sam.ps[, subjectID_var]) * omega)

    if (lI > max(Wj)) {
        stop(paste("choose lI at most", max(Wj)))
    }

    if (is.null(lC2)) {
        lC2 <- (lI - 1)
    }

    lC <- seq(lC1, lC2, by = 1)

    psi.hat.lI <- bootLongPsi(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = lI, R = R, RR = RR, T.obs.full = NULL, ncores = ncores)

    Khat.obs <- psi.hat.lI[[1]]
    T.obs <- psi.hat.lI[[2]]

    blk.size.choice <- as.list(c(1:length(lC)))

    mse.Khat.Kobs <- lapply(blk.size.choice, function(y){
        bt <- bootLongMSEPsi(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = y, R = R, RR = RR, qj = qj, Wj = Wj, Khat.obs = Khat.obs, T.obs.full = T.obs, ncores = ncores)
        return(bt)
    })

    return(mse.Khat.Kobs)
}



