#' Sample lag autocorrelation
#'
#' @param x A numeric vector.
#' @param b A numeric scalar. The block size.
#'
#' @return A matrix of autocorrelation.
#' @export
sample.lag.correlation <- function(x, b=3){
    aCovk <- numeric(0)
    autCov <- matrix(nrow = length(x), ncol = length(x))

    xbar.n.0 <- (length(x))^(-1)*sum(x[1:(length(x))])
    prod0 = 0
    for(j in 1:(length(x))){
        prod0 = prod0 + (x[(j)] - xbar.n.0)*(x[j]-xbar.n.0)
    }

    aCov0 <- (length(x))^(-1)*prod0

    for(k in 1:(b+1)){
        xbar.n.k <- (length(x)-k)^(-1)*sum(x[1:(length(x)-k)])
        prod = 0
        for(j in 1:(length(x)-k)){
            prod = prod + (x[(j+k)] - xbar.n.k)*(x[j]-xbar.n.k)
        }

        aCovk[k] <- (length(x)-k)^(-1)*prod
    }

    for(k in (b+2):(length(x)-1)){
        aCovk[k] <- 0
    }

    aCovk <- c(aCov0, aCovk)
    aCork <- aCovk/aCov0
    acf.res <- aCork
    workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
    for (rw in 1:length(acf.res)) {
        for (nc in 1:length(acf.res)) {
            workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
        }
    }

    return(workCorr)
}
