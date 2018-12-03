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

    if(q <= b){
        b <- q
    }
    maxk <- (b-1)

    autoCorr <- function(gk, first.ind, b, q, x){

        X.boot <- lapply(first.ind, function(y){
                    y:(y + b - 1)
                })

        X.boot <- x[unlist(X.boot)[1:q]]

        if(all(X.boot == 0)){
            if(gk == 0){
                rt <- 1
            }else{
                rt <- 0
            }
        }else{
            if(gk == 0){
                rt <- 1
            }else{
                rt <- (sum((X.boot[-c((q-gk+1):q)]-mean(X.boot))*(X.boot[-c(1:gk)]-mean(X.boot)))/(q))/(sum((X.boot-mean(X.boot))^2)/q)
            }
        }

        return(rt)

    }

    if(all(x == 0) | maxk == 0){
        workCorr <- diag(q)
    }else{
        workCorr <- matrix(-1, nrow = q, ncol = q)
        L <- q - b + 1

        L0 <- ceiling(q/b)

        ite = 0

        while (det(workCorr) <= 0 | ite < 51) {
            first.ind <- sample(1:L, L0, replace = TRUE)
            corK <- lapply(1:maxk, FUN = autoCorr, first.ind = first.ind, b = b, q = q, x = x)
            corK <- unlist(corK)

            for(k in b:(q-1)){
                corK[k] <- 0
            }

            corK <- c(1, corK)
            acf.res <- corK

            for (rw in 1:length(acf.res)) {
                for (nc in 1:length(acf.res)) {
                    workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
                }
            }

            ite <- ite + 1
        }

        if(ite >= 50){
            workCorr <- diag(q)
        }
    }


    return(workCorr)

}
