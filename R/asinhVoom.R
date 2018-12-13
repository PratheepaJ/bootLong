#' Modifications to \code{\link[limma]{voom}}
#'
#'Computes the observation level weight based on the mean-variance relationship.
#'
#'
#' @param counts A numeric matrix of otu_table.
#' @param design A design matrix with rows corresponding to samples and columns to coefficients to be estimated. Defaults to the unit vector meaning that samples are treated as replicates.
#' @param sj A numeric vector containing library size normalization factor based on the median-ratio method.
#' @param span A numeric. Width of the lowess smoothing window as a proportion.
#' @param plot A logical, should a plot of the mean-variance trend be displayed?
#'
#' @return A matrix of variance stablized counts.
#' @export
#' @importFrom limma lmFit
#' @importFrom stats lowess approxfun
#' @importFrom graphics title lines
asinhVoom <- function(counts, design = NULL, sj, span = 0.5,
    plot = FALSE){

    inv_asinh <- function(x) {
        y <- 0.5*exp(-x)*(exp(2*x)-1)
        return(y)
    }

    out <- list()

    counts <- as.matrix(counts)

    n <- nrow(counts)

    if (n < 2L){
        stop("Need at least two genes to fit a mean-variance trend")
    }


    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

    y <- t(asinh(t(counts)/sj))

    fit <- lmFit(y, design)

    if (is.null(fit$Amean)){
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    }


    sx <- fit$Amean
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(counts) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }

    l <- lowess(sx, sy, f = span)

    if (plot) {
        plot(sx, sy, xlab = "asinh(count)", ylab = "Sqrt( standard deviation )",
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }

    f <- approxfun(l, rule = 2)

    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, j,
            drop = FALSE])
    } else {
        fitted.values <- fit$coef %*% t(fit$design)
    }

    fitted.adj.library <- inv_asinh(fitted.values)
    fitted.count <- t(t(fitted.adj.library)*(sj))
    fitted.logcount <- asinh(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)

    out$E <- y
    out$weights <- w
    out$design <- design

    if (is.null(out$targets)){
        out$targets <- data.frame(sj = sj)
    }else{
        out$targets$sj <- sj
    }


    methods::new("EList", out)
}
