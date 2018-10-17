#' bootLongMSEPsi
#'
#' Compute MSE in computing \eqn{\psi} = two-sided probability with different block sizes.
#'
#' @param qj A numeric vector. The number of repeated observations for j-th subject.
#' @param Wj A numeric vector. The number of repeated observations in the subset for j-th subject
#' @param Khat.obs A numeric vector. The second element of the output of \code{bootLongPsi} evaluated using the initial block size and full data.
#' @inheritParams bootLongPsi
#'
#' @return A list of MSE computed with given block size b, Khat with all subsamples, Khat with initial block size.
#'
#'
#' @export
#' @importFrom parallel mclapply
bootLongMSEPsi <- function(ps, main_factor, time_var, subjectID_var, sampleID_var, b, R, RR, qj, Wj, Khat.obs = NULL, T.obs.full = NULL, ncores){
    if (is.null(Khat.obs)) {
        stop("User needs to run bootLongPsi() function with an initial block length ")
    }
    if (is.null(T.obs.full)) {
        stop("User needs to provide observed test statistic")
    }

    sam.ps <- sample_data(ps) %>% data.frame

    sam.ps[, time_var] <- as.numeric(sam.ps[, time_var])

    sam.ps[, subjectID_var] <- factor(sam.ps[, subjectID_var], levels = unique(sam.ps[, subjectID_var]))

    g <- sam.ps[, subjectID_var]

    sam.ps.lst <- split(sam.ps, g)

    num.sub.sam <- max(qj) - max(Wj) + 1

    sam.ps.q.W <- mapply(sam.ps.lst, as.list(qj), as.list(Wj), FUN = list, SIMPLIFY = FALSE)

    sam.ps.q.W.or <- lapply(sam.ps.q.W, function(x) {
        x[[1]] <- dplyr::arrange_(x[[1]], time_var)
        return(x)
    })

    if (num.sub.sam < 5) {
        stop("decrease omega")
    }

    ps.sub <- list()
    for (i in 1:num.sub.sam) {
        sub.sam.i <- lapply(sam.ps.q.W.or, function(x) {
            xd <- x[[1]]
            W <- x[[3]]
            ss <- dplyr::slice(x[[1]], i:(W + i - 1)) %>% data.frame
            return(ss)
        })
        sub.sam.i <- do.call("rbind", sub.sam.i) %>% data.frame
        rownames(sub.sam.i) <- sub.sam.i[, sampleID_var] %>% as.character
        subsam.id <- sub.sam.i[, sampleID_var] %>% as.character
        #subsam.id <- as.character(subsam.id)
        ps.sub[[i]] <- prune_samples(subsam.id, ps)
    }

    Khat <- lapply(ps.sub, function(y){
        k.hat <- bootLongPsi(ps = y, main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b, R = R, RR = RR, T.obs.full = T.obs.full, ncores = ncores)
        k.hat.v <- k.hat[[1]]
        return(k.hat.v)
    })

    # Khat <- list()
    # for(i in 1:length(ps.sub)){
    #     k.hat <- bootLongPsi(ps = ps.sub[[i]], main_factor = main_factor, time_var = time_var, subjectID_var = subjectID_var, sampleID_var = sampleID_var, b = b, R = R, RR = RR, T.obs.full = T.obs.full, ncores = ncores)
    #     Khat[[i]] <- k.hat[[1]]
    # }



    rm(ps)
    rm(sam.ps)
    rm(sam.ps.q.W)
    rm(sam.ps.q.W.or)
    rm(ps.sub)

    Khat.squared.diff <- lapply(Khat, FUN = function(w) {
        (w - Khat.obs)^2
    })

    Khat.squared.diff.df <- do.call("cbind", Khat.squared.diff)

    MSE_i <- apply(Khat.squared.diff.df, 1, FUN = function(x) {
        mean(x)
    })

    rm(Khat.squared.diff.df)

    rt <- list(MSE_i = MSE_i, Khat = Khat, Khat.obs = Khat.obs)

    gc(reset = TRUE)

    return(rt)
}

