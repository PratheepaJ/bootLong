#' bootLongSubsampling
#'
#' Finding optimal block size
#'
#' @param R A numeric. The number of block bootstrap realizations.
#' @param RR A numeric. The number of double block bootstrap realizations.
#' @param lI A numeric. The initial block size.
#' @param omega A numeric. The proportion between 0 and 1 of repeated observations for subsampling.
#' @param lC1 A numeric. (optional) any block size to start with. defualt is 1.
#' @param lC2 A numeric. (optional) any block size to end with. defualt is (lI-1).
#' @inheritParams computeStat
#'
#' @return A list of ``MSE`` in calculating K with block sizes lC=1:(lI-1), ``Khat`` values with block sizes lC=1:(lI-1) and ``Khat.obs`` with lI
#' @export
#' @import "BiocParallel"
#' @import "doParallel"
#' @import "parallel"

bootLongSubsampling = function(ps,
                               main_factor,
                               time_var,
                               subjectID_var,
                               R,
                               RR,
                               lI,
                               omega=.6,
                               lC1=1, lC2=NULL){

        qj = table(sample_data(ps)[,subjectID_var])

        Wj = trunc(table(sample_data(ps)[,subjectID_var])*omega)

        if(lI > max(Wj)){
            stop(paste("choose lI at most", max(Wj)))
            }

        if(is.null(lC2)){
            lC2 = (lI-1)
            }

        lC = seq(lC1, lC2, by=1)

        psi.hat.lI  = bootLongPsi(ps=ps,b=lI,R=R,RR=RR,main_factor=main_factor,time_var=time_var,T.obs.full=NULL,subjectID_var=subjectID_var)

        Khat.obs = psi.hat.lI[[1]]
        T.obs = psi.hat.lI[[2]]

        mseKhatKobs = lapply(seq_len(length(lC)),function(y){bootLongMSEPsi(ps=ps,qj=qj,Wj=Wj,b=y,R=R,RR=RR,main_factor=main_factor,time_var=time_var,Khat.obs=Khat.obs,T.obs.full=T.obs,subjectID_var = subjectID_var)})
        return(mseKhatKobs)
}
