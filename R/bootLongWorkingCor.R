x <- as.numeric(arima.sim(n=10, list(ar=c(.5))))
b <- 3
bootLongWorkingCor <- function(x, b){

    doit <- TRUE

    while(doit){
        # x <- as.numeric(arima.sim(n=10, list(ar=c(.5))))
        #
        # b <- 2

        corrK <- numeric(0)

        k <- (b-1)
        y.list <- lapply(seq_len((length(x)-k)), function(i){
            Yi <- x[i]

            for(j in 2:(k+2)){
                Yi[j] <- x[i]*x[i+(j-2)]
            }
            return(Yi)
        })

        y.blks <- length(y.list) - b + 1

        sample.blks.first.ind <- sample(1:y.blks, y.blks, replace = TRUE)
        sample.blks.first.ind.list <- as.list(sample.blks.first.ind)
        y.sample.list <- lapply(sample.blks.first.ind.list, function(yst){
            yst:(yst+b-1)
        })
        y.sample <- y.sample.list %>% unlist
        y.sample <- y.sample[1:length(y.list)]

        y.star.list <- y.list[y.sample]

        y.star.mat <- do.call(rbind, y.star.list)

        y.star.mat.colmean <- colMeans(y.star.mat)

        covK0 <- y.star.mat.colmean[2] - (y.star.mat.colmean[1])^2

        covK <- numeric(0)
        for(kk in 1:k){
            covK[kk] <- y.star.mat.colmean[(kk+2)] - (y.star.mat.colmean[1])^2
        }


        for(k in b:(length(x)-1)){
            covK[k] <- 0
        }

        covK <- c(covK0, covK)

        corK <- covK/covK0


        acf.res <- corK
        workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
        for (rw in 1:length(acf.res)) {
            for (nc in 1:length(acf.res)) {
                workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
            }
        }

        val <- det(workCorr)
        all.less.one <- all(abs(workCorr) <= 1)

        doit <- (val < 0) | (!all.less.one)

    }


    return(workCorr)

}


# covK0 <- y.star.mat.colmean[2] - (y.star.mat.colmean[1])^2
# covK1 <- y.star.mat.colmean[3] - (y.star.mat.colmean[1])^2
# covK2 <- y.star.mat.colmean[4] - (y.star.mat.colmean[1])^2


# for(k in 1:(b-1)){
#
#     y.list <- lapply(seq_len((length(df$x)-k)), function(i){
#         c(df$x[i], df$x[i]^2, df$x[i]*df$x[(i+k)])
#     })
#
#     y.blks <- length(y.list) - b + 1
#
#     sample.blks.first.ind <- sample(1:y.blks, y.blks, replace = TRUE)
#     sample.blks.first.ind.list <- as.list(sample.blks.first.ind)
#     y.sample.list <- lapply(sample.blks.first.ind.list, function(yst){
#         yst:(yst+b-1)
#     })
#     y.sample <- y.sample.list %>% unlist
#     y.sample <- y.sample[1:length(y.list)]
#
#     y.star.list <- y.list[y.sample]
#
#     y.star.mat <- do.call(rbind, y.star.list)
#
#     y.star.mat.colmean <- colMeans(y.star.mat)
#
#     val <- (y.star.mat.colmean[3]-y.star.mat.colmean[1]^2)/(y.star.mat.colmean[2]-y.star.mat.colmean[1]^2)
#
#     if((y.star.mat.colmean[2] > (y.star.mat.colmean[1])^2) & (abs(val) <= 1)){
#         corrK[k] <- (y.star.mat.colmean[3]-y.star.mat.colmean[1]^2)/(y.star.mat.colmean[2]-y.star.mat.colmean[1]^2)
#     }else {
#         corrK[k] <- 0
#     }
#
#
#
# }
#
# for(k in b:(length(x)-1)){
#     corrK[k] <- 0
# }
#
# corrK <- c(1, corrK)
#
# acf.res <- corrK
# workCorr <- matrix(nrow = length(acf.res), ncol = length(acf.res))
# for (rw in 1:length(acf.res)) {
#     for (nc in 1:length(acf.res)) {
#         workCorr[rw, nc] <- acf.res[(abs(rw - nc) + 1)]
#     }
# }
#
#
#
