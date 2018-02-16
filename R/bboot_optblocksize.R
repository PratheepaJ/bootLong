#' bboot_optblocksize
#'
#' Finding optimal block size
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param R  Number of block bootstrap realization.
#' @param RR Number of double block bootstrap realization.
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
bboot_optblocksize <- function(ps,R,RR,factors,time,subjectidvar="SubjectID",lI,omega=.6){
    doParallel::registerDoParallel(parallel::detectCores())
    BiocParallel::register(BiocParallel::DoparParam())

    names(sample_data(ps))[names(sample_data(ps))==subjectidvar] <- "SubjectID"
    #names(sample_data(ps))[names(sample_data(ps))==time] <- "Time"

    #   number of repeated samples per subject
    qj <- table(sample_data(ps)$SubjectID)

    #   number of subsamples per subject
    Wj <- trunc(table(sample_data(ps)$SubjectID)*omega)

    #   lI must be less than the largest subsample
    if(lI>=max(Wj)){stop("choose smaller initial block size")}

    #   lC
    lC <- seq(1,(lI-1),by=1)

    #   compute K with initial block size and all whole data
    computeK.res  <- ComputeK(ps=ps,b=lI,R=R,RR=RR,factors=factors,time=time)

    Khat.obs <- computeK.res[[1]]
    T.obs <- computeK.res[[2]]

    mseKhatKobs <- bplapply(seq_len(length(lC)),function(y){ComputeMSE(ps=ps,qj=qj,Wj=Wj,b=lC[y],R=R,RR=RR,factors=factors,time=time,Khat.obs=Khat.obs,T.obs.full=T.obs)})
    return(mseKhatKobs)
}
