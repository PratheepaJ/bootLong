compute_stat_expand <- function(ps,factors){
    ot <- as.matrix(otu_table(ps))
    anno <- data.frame(tax_table(ps))
    samdf <- data.frame(sample_data(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samdf)

    #   to estimate the size factors
    pse <- ps
    otu_table(pse) <- otu_table(ps)+1
    pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
    rm(pse)
    pDE <- DESeq2::estimateSizeFactors(pDE)
    sj <- sizeFactors(pDE)
    #   NOTE: if we want to use gene specific size factor, then
    #sij <- normalizationFactors(pDE) # matrix

    v <- voom_arcsine(counts=dgeList, design=mm, lib.size=sj,plot = F)

    #   estimate coefficients
    #   we have voom transformed data
    #   Use limma::lmFit to fit linear model using computed weights
    #   fit linear model for each taxa
    fit <- limma::lmFit(v, design = mm)

    fit.ebayes <- limma::eBayes(fit)

    res <- data.frame(fit.ebayes$lods)

    ASV <- taxa_names(ps)
    res <- bind_cols(res,pvalue=fit.ebayes$p.value[,2],ASV=ASV)
    res <- res[,c(2,3,4)]
    adj.pvalue <- p.adjust(res$pvalue,method = "BH")
    res$adj.pvalue <- adj.pvalue
    names(res)[1] <- "lfc"
    rm(ps)
    return(res)
}
