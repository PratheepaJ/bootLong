#' bootLongSubsampling
#'
#' Finding optimal block size
#'
#' @param ps a \code{phyloseq} object.
#' @param R  numeric. Number of block bootstrap realizations.
#' @param RR numeric. Number of double block bootstrap realization.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param time Time variable at repeated observations.
#' @param subjectidvar variable name for SubjectID
#' @param lI initial block size
#' @param omega percentage of repeated observations for subsampling
#'
#' @return list of ``MSE`` in calculating K with block sizes lC=1:(lI-1), ``Khat`` values with block sizes lC=1:(lI-1) and ``Khat.obs`` with lI
#' @export
#' @import "BiocParallel"
#' @import "doParallel"
#' @import "dplyr"
#' @import "parallel"
bootLongSubsampling <- function(ps,R,RR,factors,time,SubjectID_n="SubjectID",lI,omega=.6,lC1=1,lC2=NULL){
    # doParallel::registerDoParallel(parallel::detectCores())
    # BiocParallel::register(BiocParallel::DoparParam())

    #   number of repeated samples per subject
    qj <- table(sample_data(ps)[,SubjectID_n])

    #   number of subsamples per subject
    Wj <- trunc(table(sample_data(ps)[,SubjectID_n])*omega)

    #   lI must be less than the largest number of subsample
    if(lI > max(Wj)){stop(paste("choose lI at most",max(Wj)))}

    if(is.null(lC2)){lC2 <- (lI-1) }
    #   lC
    lC <- seq(lC1,lC2,by=1)

    #   compute K with initial block size and all whole data
    psi.hat.lI  <- bootLongPsi(ps=ps,b=lI,R=R,RR=RR,factors=factors,time=time,T.obs.full=NULL,SubjectID_n=SubjectID_n)

    Khat.obs <- psi.hat.lI[[1]]
    T.obs <- psi.hat.lI[[2]]

    mseKhatKobs <- lapply(seq_len(length(lC)),function(y){bootLongMSEPsi(ps=ps,qj=qj,Wj=Wj,b=y,R=R,RR=RR,factors=factors,time=time,Khat.obs=Khat.obs,T.obs.full=T.obs,SubjectID_n=SubjectID_n)})
    return(mseKhatKobs)
}
