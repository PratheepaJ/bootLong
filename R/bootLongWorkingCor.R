#' compute working correlation using block bootstrap
#'
#' @param x A numeric vector.
#' @inheritParams computeStat
#'
#' @return A numeric matrix. The working autocorrelation for GEE using block bootstrap method.
#' @export
bootLongWorkingCor <- function(x, b){

    ## naive method
    q <- length(x)
    k <- (b-1)

    autoCorr <- function(gk, first.ind){

        X.boot <- lapply(first.ind, function(y){
                    y:(y + b - 1)
                })
        X.boot <- x[unlist(X.boot)[1:q]]

        if(gk == 0){
            rt <- 1
        }else{
            rt <- (sum((X.boot[-c((q-gk+1):q)]-mean(X.boot))*(X.boot[-c(1:gk)]-mean(X.boot)))/(q))/(sum((X.boot-mean(X.boot))^2)/q)
        }

        return(rt)

    }

    workCorr <- matrix(-1, nrow = length(acf.res), ncol = length(acf.res))
    while (det(workCorr) <= 0) {
        L <- q - b + 1
        L0 <- ceiling(q/b)
        first.ind <- sample(1:L, L0, replace = TRUE)
        corK <- lapply(1:k, FUN = autoCorr, first.ind = first.ind)
        corK <- unlist(corK)

        for(k in b:(length(x)-1)){
            corK[k] <- 0
        }

        corK <- c(1, corK)
        acf.res <- corK
        workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
        for (rw in 1:length(acf.res)) {
            for (nc in 1:length(acf.res)) {
                workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
            }
        }


    }

    return(workCorr)

}
