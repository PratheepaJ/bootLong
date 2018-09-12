asinh_voom = function (counts,
                       design = NULL,
                       sj,
                       normalize.method = "none",
                       span = 0.5,
                       plot = FALSE,
                       save.plot = FALSE, ...){

        inv_asinh <- function(x) {
            y <- 0.5*exp(-x)*(exp(2*x)-1)
            return(y)
        }

        out <- list()

        # if (is(counts, "DGEList")) {
        #     out$genes <- counts$genes
        #     out$targets <- counts$samples
        #     if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
        #         0)
        #         design <- model.matrix(~group, data = counts$samples)
        #     if (is.null(sj))
        #         sj <- counts$samples$sj * counts$samples$norm.factors
        #     counts <- counts$counts
        # }
        # else {
        #     isExpressionSet <- suppressPackageStartupMessages(is(counts,
        #                                                          "ExpressionSet"))
        #     if (isExpressionSet) {
        #         if (length(Biobase::fData(counts)))
        #             out$genes <- Biobase::fData(counts)
        #         if (length(Biobase::pData(counts)))
        #             out$targets <- Biobase::pData(counts)
        #         counts <- Biobase::exprs(counts)
        #     }
        #     else {
        #         counts <- as.matrix(counts)
        #     }
        # }
        counts <- as.matrix(counts)

        n <- nrow(counts)
        if (n < 2L)
            stop("Need at least two genes to fit a mean-variance trend")
        if (is.null(design)) {
            design <- matrix(1, ncol(counts), 1)
            rownames(design) <- colnames(counts)
            colnames(design) <- "GrandMean"
        }
        # if (is.null(sj))
        #     sj <- colSums(counts)
        y <- t(asinh(t(counts)*sj))
        y <- normalizeBetweenArrays(y, method = normalize.method)
        fit <- lmFit(y, design, ...)
        if (is.null(fit$Amean))
            fit$Amean <- rowMeans(y, na.rm = TRUE)
        sx <- fit$Amean - mean(asinh(sj))
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
            fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
                                                                          j, drop = FALSE])
        }
        else {
            fitted.values <- fit$coef %*% t(fit$design)
        }
        fitted.adj.library <- inv_asinh(fitted.values)
        fitted.count <- t(t(fitted.adj.library) / (sj))
        fitted.logcount <- asinh(fitted.count)
        w <- 1/f(fitted.logcount)^4
        dim(w) <- dim(fitted.logcount)
        out$E <- y
        out$weights <- w
        out$design <- design
        if (is.null(out$targets))
            out$targets <- data.frame(sj = sj)
        else out$targets$sj <- sj
        if (save.plot) {
            out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )",
                                ylab = "Sqrt( standard deviation )")
            out$voom.line <- l
        }

    new("EList", out)
}
