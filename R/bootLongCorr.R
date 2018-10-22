x <- seq(1,10,by=1)
df <- data.frame(x = x, Index = 1:length(x))

Ly0 <- length(x)
first.index.Ly0 <- sample(1:Ly0, Ly0, replace = TRUE)
first.index.Ly0.lst <- as.list(first.index.Ly0)
x.index.0 <- df$Index[lapply(first.index.Ly0.lst, function(y){y:(y+0)}) %>% unlist]
x.star.0 <- df$x[x.index.0[1:length(x)]]

x.star <- list()
for(k in 1:(b-1)){
    Ly <- length(x) - (k+1) + 1
    first.index.Ly <- sample(1:Ly, Ly, replace = TRUE)
    first.index.Ly.lst <- as.list(first.index.Ly)
    x.index <- df$Index[lapply(first.index.Ly.lst, function(y){y:(y+k)}) %>% unlist]
    x.star[[k]] <- df$x[x.index[1:length(x)]]
}
######

aCovk <- numeric(0)
xbar.n.0 <- (length(x.star.0))^(-1)*sum(x.star.0[1:(length(x.star.0))])
prod0 <- 0
for(j in 1:(length(x.star.0))){
    prod0 <- prod0 + (x.star.0[(j)] - xbar.n.0)*(x.star.0[j]-xbar.n.0)
}
aCov0 <- (length(x.star.0))^(-1)*prod0


aCovk <- list()

aCovk <- lapply(seq_len((b-1)), function(k){
    xbar.n.k <- (length(x.star[[k]])-k)^(-1)*sum(x.star[[k]][1:(length(x.star[[k]])-k)])
    prod = 0
    for(j in 1:(length(x.star[[k]])-k)){
        prod <- prod + (x.star[[k]][(j+k)] - xbar.n.k)*(x.star[[k]][j]-xbar.n.k)
    }

    aCovk.val <- (length(x.star[[k]])-k)^(-1)*prod

    return(aCovk.val)
})

aCovk <- aCovk %>% unlist

for(k in b:(length(x)-1)){
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

det(workCorr)
