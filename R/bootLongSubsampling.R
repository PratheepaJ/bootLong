#' bootLongSubsampling
#'
#' Finding optimal block size.
#'
#' @param lI A numeric. The initial block size.
#' @param omega  A numeric. The proportion between 0 and 1 of repeated observations for subsampling.
#' @param lC1 A numeric. (optional) any block size to start with. Defualt is 1.
#' @param lC2 A numeric. (optional) any block size to end with. defualt is (lI-1).
#' @param psi.hat.lI A logical. psi is estiamted already with the initial block size.
#' @param psi.hat.lI.val A list. If psi is calculated already, give the list of values.
#' @inheritParams bootLongMSEPsi
#'
#' @return A list of ``MSE`` in calculating K with block sizes lC=1:(lI-1), ``Khat`` values with block sizes lC=1:(lI-1) and ``Khat.obs`` with lI
#' @export
#' @importFrom parallel mclapply
bootLongSubsampling <- function(ps, main_factor, time_var, subjectID_var, sampleID_var, lI, R, RR, omega = .6, lC1 = 1, lC2 = NULL, ncores, psi.hat.lI = FALSE, psi.hat.lI.val = NULL){

    # doParallel::registerDoParallel(parallel::detectCores())
    # BiocParallel::register(BiocParallel::DoparParam())

    if (dim(otu_table(ps))[1] == nsamples(ps)) {
        otu_table(ps) <- t(otu_table(ps, taxa_are_rows = T))
    }

    sam.ps <- sample_data(ps) %>% data.frame

    if(!(all(as.character(sam.ps[, sampleID_var]) == sample_names(ps)))){
        stop(paste0(sampleID_var, " must be same as sample names in the phyloseq"))
    }

    sam.ps[, subjectID_var] <- factor(sam.ps[, subjectID_var], levels = unique(sam.ps[, subjectID_var]))


    qj <- table(sam.ps[, subjectID_var])

    Wj <- trunc(table(sam.ps[, subjectID_var]) * omega)

    if (lI > max(Wj)) {
        stop(paste("choose lI at most", max(Wj)))
    }

    if (is.null(lC2)) {
        lC2 <- (lI - 1)
    }

    lC <- seq(lC1, lC2, by = 1)

    if(psi.hat.lI){# if this has been calculated already and given as psi.hat.lI.val

        Khat.obs <- psi.hat.lI.val[[1]]
        T.obs <- psi.hat.lI.val[[2]]

        blk.size.choice <- as.list(lC)

        mse.Khat.Kobs <- lapply(blk.size.choice, function(y){
            bt <- bootLongMSEPsi(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = y, R = R, RR = RR, qj = qj, Wj = Wj, Khat.obs = Khat.obs, T.obs.full = T.obs, ncores = ncores)
            return(bt)
        })

        rt <- mse.Khat.Kobs

    }else{
        psi.hat.lI <- bootLongPsi(ps = ps, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = lI, R = R, RR = RR, T.obs.full = NULL, ncores = ncores)

        rt <- psi.hat.lI

    }

    return(rt)
}



