#' arcsinhTransform
#'
#' This is a modified \code{voom} utility function used by the higher level function \code{computeStat()}
#' to compute arcsinh transformation and weights for the \code{glm.nb} and \code{geem}.
#'
#' @param counts \code{DGEList}.
#' @param design a model matrix created by \code{computeStat()}.
#' @param lib.size numeric vector containing library size normalization factor computed by \code{computeStat()}.
#' @param span width of the lowess smoothing window as a proportion.

#' @return An EList object. E is numeric matrix of transformed counts on the arcsine scale
#' @export
arcsinhTransform <- function (counts, design = NULL, lib.size = NULL,span = 0.5,...){
        #   arcsine transformation
        arcs <- function(x) {
        transformed <- log(x + sqrt(x ^ 2 + 1))
        return(transformed)
        }

        #   inverse of arcsine transformation
        inv_asinh <- function(x) {
        y <- 0.5*exp(-x)*(exp(2*x)-1)
        return(y)
        }

        out <- list()

        if (is(counts, "DGEList")){
            out$genes <- counts$genes
            out$targets <- counts$samples
            if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
                0)
                design <- model.matrix(~group, data = counts$samples)
            if (is.null(lib.size))
                lib.size <- with(counts$samples, lib.size * norm.factors)
            counts <- counts$counts
        }else{
        isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
        if (isExpressionSet) {
            if (length(Biobase::fData(counts)))
                out$genes <- Biobase::fData(counts)
            if (length(Biobase::pData(counts)))
                out$targets <- Biobase::pData(counts)
            counts <- Biobase::exprs(counts)
        }
        else {
            counts <- as.matrix(counts)
        }
        }

        n <- nrow(counts)

        if (n < 2L){stop("Need at least two genes to fit a mean-variance trend")}

        if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
        }

        if (is.null(lib.size)){lib.size <- colSums(counts)}

        y <- t(arcs(t(counts)/lib.size))

        fit <- lmFit(y, design, ...)

        if (is.null(fit$Amean)){fit$Amean <- rowMeans(y, na.rm = TRUE)}

        #   fitting mean variance dependence
        #   fitted mean for each row (taxa)
        sx <- fit$Amean
        #   fitted standard deviation
        sy <- sqrt(fit$sigma)

        allzero <- rowSums(counts) == 0

        if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
        }

        l <- lowess(sx, sy, f = span)

        #   linearly interpolate given the smooth line from the lowess
        f <- approxfun(l, rule = 2)

        if(fit$rank < ncol(design)){
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
                                                                      j, drop = FALSE])
        }else{
        fitted.values <- fit$coef %*% t(fit$design)
        }

        fitted.cpm <- inv_asinh(fitted.values)

        fitted.count <- t(t(fitted.cpm)*(lib.size))

        fitted.arcs <- arcs(fitted.count)

        w <- 1/f(fitted.arcs)^4

        dim(w) <- dim(fitted.arcs)

        out$E <- y

        out$weights <- w

        out$design <- design

        if(is.null(out$targets)){
        out$targets <- data.frame(lib.size = lib.size)
        }else{
        out$targets$lib.size <- lib.size
        }

        new("EList", out)
}
