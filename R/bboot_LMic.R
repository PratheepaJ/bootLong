#' bboot_LMic
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param b numeric, optimal block size to account for dependence within-subject.
#' @param R Number of block bootstrap realization.
#' @param RR Number of double block bootstrap realization.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param time Time variable at repeated observations.
#' @param FDR False discovery rate
#'
#' @return list of dataframe with ASV, observed stat, adjusted pvalues, lcl, ucl, observed pivotal quantity; stat.obs; stat.star; stat.star.star; T.star_obs
#' @export
#'
bboot_LMic <- function(ps,b,R,RR,factors,time,FDR=.1){
    #   otu table of observed phyloseq: rows taxa; columns samples
    if(dim(otu_table(ps))[1]==nsamples(ps)){
        otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
    }

    res.obs <- compute_stat(ps,factors)
    stat.name <- colnames(res.obs)[1]

    boot.results <- list()

    boot.results <- lapply(seq_len(R),FUN=function(i){
        ps.boot <- bboot_phyloseq(ps,b,time)
        ps.boot <- ps.boot[[1]]

        df.boot <- compute_stat(ps.boot,factors)

        #   double MBB
        boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
            ps.boot.bb <- bboot_phyloseq(ps.boot,b,time)
            ps.boot.bb <- ps.boot.bb[[1]]
            df.boot.bb <- compute_stat(ps.boot.bb,factors)
            rm(ps.boot.bb)
            return(df.boot.bb)
        })

        rm(ps.boot)

        return(list(df.boot,boot.results.bb))

    })

    #   stat* from MBB estimates (R times)
    boot.results.all <- lapply(boot.results,"[[",1)
    #   stat** from double MBB to compute SE(stat*)
    boot.results.bb <- lapply(boot.results,"[[",2)

    rm(boot.results)

    #   stat* from MBB in a dataframe
    stat.star <- do.call("cbind",lapply(boot.results.all,FUN=function(x){x[,1]}))
    stat.star <- as.data.frame(stat.star)

    rm(boot.results.all)

    #   add observed value of stat to a dataframe
    stat.obs <- data.frame(stat.obs=res.obs[,1])

    #   compute SE(stat.obs) - from stat*
    sd.stat <- data.frame(sd.stat=apply(stat.star,1,FUN=sd,na.rm=FALSE))

    #   compute T observed
    T.obs <- data.frame(T.obs=stat.obs/sd.stat)

    #   compute   SE(stat*) from double MBB
    stat.star.star <-  lapply(boot.results.bb,function(x){
        do.call("cbind",lapply(x,FUN=function(x){x[,1]}))

    })

    sd.stat.star <- do.call("cbind",lapply(stat.star.star,function(x){data.frame(sd.stat.star=apply(x,1,FUN=sd,na.rm=TRUE))}))

    #   compute stat*-stat.obs
    T.num.star <- data.frame(apply(stat.star,2,function(x){x-stat.obs}))

    T.star <- T.num.star/sd.stat.star

    #   compute p-value
    T.star_obs <- dplyr::bind_cols(T.star ,T.obs=T.obs[,1])

    pvalue <- apply(T.star_obs,1,function(x){sum(abs(x[1:R])>=abs(x[(R+1)]))/R})

    #   adjusted p-value for multiple testing
    pvalue.adj <- data.frame(pvalue.adj=p.adjust(pvalue,method = "BH"))

    # output taxa names, stat, adj pvalues
    #txnames <- taxa_names(ps)
    txnames <- dplyr::select(res.obs,ASV)
    out <- data.frame(Taxa=txnames,stat=res.obs[,1],pvalue.adj=pvalue.adj)
    names(out)[which(names(out)=="stat")] <- stat.name
    #   compute confidence interval: not the simultaneous CI so will be wider than expected
    lcl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=FDR/2,na.rm=TRUE)})
    ucl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=(1-FDR/2),na.rm=TRUE)})

    out <- dplyr::bind_cols(out,lcl=lcl,ucl=ucl,T.obs=T.obs[,1])

    rt <- list(out,stat.obs,stat.star,stat.star.star,T.star_obs)
    return(rt)
}
