#' bootLongPhyloseq
#'
#' Creates block bootstrap realization of a \code{phyloseq}.
#'
#' @param sampleID_var Character string. The name of the sample ID variable.
#' @inheritParams computeStat
#' @return a \code{phyloseq} object. A block bootstrap realization of \code{ps}.
#' @export
bootLongPhyloseq <- function(ps, time_var, subjectID_var, sampleID_var, b) {
    if(dim(otu_table(ps))[1]==nsamples(ps)){
        otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
    }

    sam.pss <- sample_data(ps) %>% data.frame
    sam.pss$Index <- seq(1, nsamples(ps), by = 1)

    sam.pss[, subjectID_var] <- factor(sam.pss[, subjectID_var], levels = unique(sam.pss[, subjectID_var] ))

    sam.pss[, time_var] <- as.numeric(sam.pss[, time_var])

    ot <- otu_table(ps) %>% data.frame

    g <- sam.pss[, subjectID_var]
    sam.ps.split.by.subjects <- split(sam.pss, g)

    num.of.rep.obs <- lapply(sam.ps.split.by.subjects, function(x) {
        dim(x)[1]
    })

    num.of.blks <- lapply(num.of.rep.obs, function(x) {
        x - b + 1
    })

    num.of.rep.obs.max <- num.of.rep.obs %>% unlist %>% max
    L <- num.of.blks %>% unlist %>% max
    L0 <- ceiling(num.of.rep.obs.max/b)
    blks.first.index <- sample(1:L, L0, replace = TRUE)

    sampling.blks.within.subject.indices <- lapply(sam.ps.split.by.subjects,
        FUN = function(q) {
            rt <- bootLongIndices(x = q, time_var = time_var, b = b, L = L, blks.first.index = blks.first.index)
            return(rt)
        })

    boot.sample.indices <- as.numeric(unlist(sampling.blks.within.subject.indices))

    # sum(rownames(sam.pss)==sam.pss[, sampleID_var] %>% as.character())
    rownames(sam.pss) <- sam.pss[, sampleID_var] %>% as.character()
    blk.boot.ot <- ot[, boot.sample.indices]
    blk.boot.sam.ps <- sam.pss[boot.sample.indices, ]
    blk.boot.sam.ps[, sampleID_var] <- rownames(blk.boot.sam.ps) %>% as.factor

    blk.boot.sam.ps[, subjectID_var] <- factor(blk.boot.sam.ps[, subjectID_var], levels = unique(blk.boot.sam.ps[, subjectID_var]))

    g2 <- blk.boot.sam.ps[, subjectID_var]
    boot.sam.ps <- split(blk.boot.sam.ps, g2)

    boot.sam.ps <- lapply(boot.sam.ps, function(x) {
        tim <- seq(1, dim(x)[1])
        x[, time_var] <- tim
        return(x)
    })

    blk.boot.sam.ps <- do.call("rbind", boot.sam.ps)
    blk.boot.sam.ps <- dplyr::select(blk.boot.sam.ps, -Index)
    colnames(blk.boot.ot) <- rownames(blk.boot.sam.ps)
    rownames(blk.boot.ot) <- taxa_names(ps)
    # rownames(blk.boot.sam.ps) <- blk.boot.sam.ps[, sampleID_var] %>% as.character
    ps.boot <- phyloseq::merge_phyloseq(X = otu_table(blk.boot.ot, taxa_are_rows = T), sample_data(blk.boot.sam.ps), tax_table(ps))

    return(list(ps.boot))
}
