x <- as.numeric(arima.sim(n=10, list(ar=c(.5))))
df <- data.frame( x = x)
b <- 2

k <- 0

# blks <- length(x)-b+1
# sample.blk.first.ind0 <- sample(1:blks, blks, replace = TRUE)
# sample.blk.first.ind.list0 <- as.list(sample.blk.first.ind0)
# sample.blks.all <- lapply(sample.blk.first.ind.list0, function(xst){
#     xst:(xst+b-1)
# })
#
# sample.blks.all <- sample.blks.all %>% unlist
#
# sample.blks.all <- sample.blks.all[1:length(x)]
#
# x.star <- x[sample.blks.all]
#
# covK0 <- var(x.star)*(length(x)-1)/length(x)
#
# covK0

y.list0 <- lapply(seq_len((length(df$x)-k)), function(i){
    c(df$x[i], df$x[i]*df$x[(i+k)])
})


y.blks0 <- length(y.list0) - b + 1
sample.blks.first.ind0 <- sample(1:y.blks0, y.blks0, replace = TRUE)
sample.blks.first.ind.list0 <- as.list(sample.blks.first.ind0)
y.sample.list0 <- lapply(sample.blks.first.ind.list0, function(yst){
    yst:(yst+b-1)
})
y.sample0 <- y.sample.list0 %>% unlist
y.sample0 <- y.sample0[1:length(y.list0)]

y.star.list0 <- y.list0[y.sample0]

y.star.mat0 <- do.call(rbind, y.star.list0)

covK0 <- (length(x)-k)^(-1)*sum(y.star.mat0[,2]) - ((length(x)-k)^(-1)*sum(y.star.mat0[,1]))^2

covK0



covK <- numeric(0)

for(k in 1:(b-1)){
    y.list <- lapply(seq_len((length(df$x)-k)), function(i){
        c(df$x[i],df$x[i]*df$x[(i+k)])
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

    covK[k] <- (length(x)-k)^(-1)*sum(y.star.mat[,2]) - ((length(x)-k)^(-1)*sum(y.star.mat[,1]))^2

}

for(k in b:(length(x)-1)){
    covK[k] <- 0
}

covK <- c(covK0, covK)
corK <- covK/covK0

corK
