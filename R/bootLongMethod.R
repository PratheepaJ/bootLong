#' bootLongMethod
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param b numeric, optimal block size to account for dependence within-subject.
#' @param R Number of block bootstrap realization.
#' @param RR Number of double block bootstrap realization.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param time Time variable at repeated observations.
#' @param FDR False discovery rate
#'
#' @return list of dataframe with ASV, observed stat, pvalues, adjusted pvalues, lcl, ucl, observed pivotal quantity; stat.obs; stat.star; stat.star.star; T.star_obs.
#' @import "ashr"
#' @export
bootLongMethod <- function(ps,b,R,RR,factors,time,FDR=.1){
    #   otu table of observed phyloseq: rows taxa; columns samples
    if(dim(otu_table(ps))[1]==nsamples(ps)){
        otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
    }

    res.obs <- computeStat(ps,factors,time,b)

    stat.name <- colnames(res.obs)[2]

    boot.results <- list()

    boot.results <- lapply(seq_len(R),FUN=function(i){
        ps.boot <- bootLongPhyloseq(ps,b,time)
        ps.boot <- ps.boot[[1]]

        df.boot <- computeStat(ps.boot,factors,time,b)

        #   double MBB
        boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
            ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)
            ps.boot.bb <- ps.boot.bb[[1]]
            df.boot.bb <- computeStat(ps.boot.bb,factors,time,b)
            rm(ps.boot.bb)
            return(df.boot.bb)
        })

        rm(ps.boot)

        return(list(df.boot,boot.results.bb))

    })

    #   beta* from MBB estimates (R times)
    boot.results.all <- lapply(boot.results,"[[",1)
    #   beta** from double MBB to compute SE(stat*)
    boot.results.bb <- lapply(boot.results,"[[",2)

    rm(boot.results)

    #   stat* from MBB in a dataframe
    stat.star <- do.call("cbind",lapply(boot.results.all,FUN=function(x){x[,2]}))
    stat.star <- as.data.frame(stat.star)
    #   compute SE(stat.obs) - from stat*
    sd.stat <- apply(stat.star,1,FUN=sd,na.rm=FALSE)
    #   empirical Bayes shrinkage
    beta.obs <- res.obs[,2]
    shrink.beta.obs <- suppressMessages(ash(beta.obs,sebetahat = sd.stat,mixcompdist = "normal"))
    shrink.beta.est <- shrink.beta.obs$result$PosteriorMean
    shrink.beta.sd <- shrink.beta.obs$result$PosteriorSD

    rm(boot.results.all)

    #   add observed value of stat to a dataframe
    stat.obs <- data.frame(stat.obs=shrink.beta.est)

    #   compute T observed
    T.obs <- data.frame(T.obs=stat.obs/shrink.beta.sd)

    #   compute   SE(stat*) from double MBB
    stat.star.star <-  lapply(boot.results.bb,function(x){
        do.call("cbind",lapply(x,FUN=function(x){x[,2]}))

    })

    sd.stat.star <- do.call("cbind",lapply(stat.star.star,function(x){data.frame(sd.stat.star=apply(x,1,FUN=sd,na.rm=TRUE))}))
    #   empirical Bayes shrinkage
    shrink.beta.boot <- list()
    shrink.beta.boot.est <- list()
    shrink.beta.boot.sd <- list()
    for(i in 1:R){
        shrink.beta.boot[[i]] <- suppressMessages(ash(stat.star[,i],sebetahat = sd.stat.star[,i],mixcompdist = "normal"))
        shrink.beta.boot.est[[i]] <- shrink.beta.boot[[i]]$result$PosteriorMean
        shrink.beta.boot.sd[[i]] <- shrink.beta.boot[[i]]$result$PosteriorSD
    }

    stat.star <- do.call("cbind",shrink.beta.boot.est)
    sd.stat.star <- do.call("cbind",shrink.beta.boot.sd)
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
    out <- data.frame(Taxa=txnames,stat=stat.obs[,1],pvalue=pvalue,pvalue.adj=pvalue.adj)
    names(out)[which(names(out)=="stat")] <- stat.name
    #   compute confidence interval: not the simultaneous CI so will be wider than expected
    lcl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=FDR/2,na.rm=TRUE)})
    ucl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=(1-FDR/2),na.rm=TRUE)})

    out <- dplyr::bind_cols(out,lcl=lcl,ucl=ucl,T.obs=T.obs[,1])

    rt <- list(out,stat.obs,stat.star,stat.star.star,T.star_obs)
    names(rt) <- c("summary","beta.hat","beta.hat.star","beta.hat.star.star","T.obs")
    return(rt)
}

### attempt 1
# bootLongMethod <- function(ps,b,R,RR,factors,time,FDR=.1){
#     #   otu table of observed phyloseq: rows taxa; columns samples
#     if(dim(otu_table(ps))[1]==nsamples(ps)){
#         otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
#     }
#
#     res.obs <- computeStat(ps,factors)
#     stat.name <- colnames(res.obs)[1]
#
#     boot.results <- list()
#
#     boot.results <- lapply(seq_len(R),FUN=function(i){
#         ps.boot <- bootLongPhyloseq(ps,b,time)
#         ps.boot <- ps.boot[[1]]
#
#         df.boot <- computeStat(ps.boot,factors)
#
#         #   double MBB
#         boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
#             ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)
#             ps.boot.bb <- ps.boot.bb[[1]]
#             df.boot.bb <- computeStat(ps.boot.bb,factors)
#             rm(ps.boot.bb)
#             return(df.boot.bb)
#         })
#
#         rm(ps.boot)
#
#         return(list(df.boot,boot.results.bb))
#
#     })
#
#     #   stat* from MBB estimates (R times)
#     boot.results.all <- lapply(boot.results,"[[",1)
#     #   stat** from double MBB to compute SE(stat*)
#     boot.results.bb <- lapply(boot.results,"[[",2)
#
#     rm(boot.results)
#
#     #   stat* from MBB in a dataframe
#     stat.star <- do.call("cbind",lapply(boot.results.all,FUN=function(x){x[,1]}))
#     stat.star <- as.data.frame(stat.star)
#
#     rm(boot.results.all)
#
#     #   add observed value of stat to a dataframe
#     stat.obs <- data.frame(stat.obs=res.obs[,1])
#
#     #   compute SE(stat.obs) - from stat*
#     sd.stat <- data.frame(sd.stat=apply(stat.star,1,FUN=sd,na.rm=FALSE))
#
#     #   compute T observed
#     T.obs <- data.frame(T.obs=stat.obs/sd.stat)
#
#     #   compute   SE(stat*) from double MBB
#     stat.star.star <-  lapply(boot.results.bb,function(x){
#         do.call("cbind",lapply(x,FUN=function(x){x[,1]}))
#
#     })
#
#     sd.stat.star <- do.call("cbind",lapply(stat.star.star,function(x){data.frame(sd.stat.star=apply(x,1,FUN=sd,na.rm=TRUE))}))
#
#     #   compute stat*-stat.obs
#     T.num.star <- data.frame(apply(stat.star,2,function(x){x-stat.obs}))
#
#     T.star <- T.num.star/sd.stat.star
#
#     #   compute p-value
#     T.star_obs <- dplyr::bind_cols(T.star ,T.obs=T.obs[,1])
#
#     pvalue <- apply(T.star_obs,1,function(x){sum(abs(x[1:R])>=abs(x[(R+1)]))/R})
#
#     #   adjusted p-value for multiple testing
#     pvalue.adj <- data.frame(pvalue.adj=p.adjust(pvalue,method = "BH"))
#
#     # output taxa names, stat, adj pvalues
#     #txnames <- taxa_names(ps)
#     txnames <- dplyr::select(res.obs,ASV)
#     out <- data.frame(Taxa=txnames,stat=res.obs[,1],pvalue=pvalue,pvalue.adj=pvalue.adj)
#     names(out)[which(names(out)=="stat")] <- stat.name
#     #   compute confidence interval: not the simultaneous CI so will be wider than expected
#     lcl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=FDR/2,na.rm=TRUE)})
#     ucl <- apply(stat.star,1,FUN=function(x){quantile(x,probs=(1-FDR/2),na.rm=TRUE)})
#
#     out <- dplyr::bind_cols(out,lcl=lcl,ucl=ucl,T.obs=T.obs[,1])
#
#     rt <- list(out,stat.obs,stat.star,stat.star.star,T.star_obs)
#     names(rt) <- c("summary","beta.hat","beta.hat.star","beta.hat.star.star","T.obs")
#     return(rt)
# }


