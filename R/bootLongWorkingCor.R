#' compute working correlation using block bootstrap
#'
#' @param x A numeric vector.
#' @inheritParams computeStat
#'
#' @return A numeric matrix. The working autocorrelation for GEE using block bootstrap method.
#' @export
bootLongWorkingCor <- function(x, b){

    # corrK <- numeric(0)
    #
    # y.list0 <- lapply(seq_len((length(x))), function(i){
    #
    #     Yi0 <- c(x[i], x[i]^2)
    #
    #     return(Yi0)
    # })
    #
    # y.star.mat0 <- do.call(rbind, y.list0)
    #
    # y.star.mat.colSum0 <- colSums(y.star.mat0)
    #
    # covK0 <- y.star.mat.colSum0[2] - (length(x))*((y.star.mat.colSum0[1])/length(x))^2
    #
    # covK <- numeric(0)
    #
    # for(k in 1:(b-1)){
    #     y.list <- lapply(seq_len((length(x)-k)), function(i){
    #         Yi <- x[i]
    #
    #         for(j in 2:(k+2)){
    #             Yi[j] <- x[i]*x[i+(j-2)]
    #         }
    #         return(Yi)
    #     })
    #
    #
    #     y.star.mat <- do.call(rbind, y.list)
    #
    #     y.star.mat.colSum <- colSums(y.star.mat)
    #
    #
    #     covK[k] <- y.star.mat.colSum[(k+2)] - (length(x)-k)*(y.star.mat.colSum0[1]/length(x))^2
    # }
    #
    #
    #
    # for(k in b:(length(x)-1)){
    #     covK[k] <- 0
    # }
    #
    # covK <- c(covK0, covK)
    #
    # corK <- covK/covK0
    #
    #
    # acf.res <- corK
    # workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
    # for (rw in 1:length(acf.res)) {
    #     for (nc in 1:length(acf.res)) {
    #         workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
    #     }
    # }

        ###method 2
        # corrK <- numeric(0)
        #
        # k <- (b-1)
        # y.list <- lapply(seq_len((length(x)-k)), function(i){
        #     Yi <- x[i]
        #
        #     for(j in 2:(k+2)){
        #         Yi[j] <- x[i]*x[i+(j-2)]
        #     }
        #     return(Yi)
        # })
        #
        #
        # y.star.mat <- do.call(rbind, y.list)
        #
        # y.star.mat.colmean <- colMeans(y.star.mat)
        #
        # covK0 <- y.star.mat.colmean[2] - (y.star.mat.colmean[1])^2
        #
        # covK <- numeric(0)
        # for(kk in 1:k){
        #     covK[kk] <- y.star.mat.colmean[(kk+2)] - (y.star.mat.colmean[1])^2
        # }
        #
        #
        # for(k in b:(length(x)-1)){
        #     covK[k] <- 0
        # }
        #
        # covK <- c(covK0, covK)
        #
        # corK <- covK/covK0
        #
        #
        # acf.res <- corK
        # workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
        # for (rw in 1:length(acf.res)) {
        #     for (nc in 1:length(acf.res)) {
        #         workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
        #     }
        # }

        #method 1
        acf.res <- as.numeric(acf(x, plot = F,lag.max = length(x))$acf)
        #acf.res[(b+1):length(acf.res)] <- 0
        workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
        for (rw in 1:length(acf.res)) {
            for (nc in 1:length(acf.res)) {
                workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
            }
        }

    return(workCorr)

}
