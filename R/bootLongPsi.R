#' bootLongPsi
#'
#' Compute \eqn{\psi} = two-sided probability with a given block size
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param b numeric, block size to account for dependence within-subject.
#' @param R Number of block bootstrap realization.
#' @param RR Number of double block bootstrap realization.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param time Time variable at repeated observations.
#' @param T.obs.full If observed statistic is already computed
#'
#' @return a list with first element ``K.val`` - two-sided significance probability.
#'         second element ``observed statistic``
#'
#' @export
bootLongPsi <- function(ps,b,R,RR,factors,time,T.obs.full=NULL){
    #   otu table of observed phyloseq: rows taxa; columns samples
    if(dim(otu_table(ps))[1]==nsamples(ps)){
        otu_table(ps) <- t(otu_table(ps,taxa_are_rows = T))
    }

    #   compute observed statistic
    res.obs <- computeStat(ps=ps,factors=factors)

    boot.results <- list()

    # boot.results <- lapply(seq_len(R),FUN=function(i){
    #     ps.boot <- bootLongPhyloseq(ps,b,time)[[1]]
    # 
    # 
    #     df.boot <- computeStat(ps.boot,factors)
    # 
    #     #   double MBB
    #     boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
    #         ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)[[1]]
    # 
    #         df.boot.bb <- computeStat(ps.boot.bb,factors)
    #         rm(ps.boot.bb)
    #         return(df.boot.bb)
    #     })
    # 
    #     rm(ps.boot)
    # 
    #     return(list(df.boot,boot.results.bb))
    # 
    # })

    boot.results <- lapply(seq_len(R),FUN=function(i){
            ps.boot <- bootLongPhyloseq(ps,b,time)
            ps.boot <- ps.boot[[1]]
            
            df.boot <- computeStat(ps=ps.boot,factors=factors,time=time,b=b)
            
            #   double MBB
            boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
                    ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)
                    ps.boot.bb <- ps.boot.bb[[1]]
                    df.boot.bb <- computeStat(ps=ps.boot.bb,factors=factors,time=time,b=b)
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

    # #   stat* from MBB in a dataframe
    # stat.star <- do.call("cbind",lapply(boot.results.all,FUN=function(x){x[,1]}))
    # stat.star <- as.data.frame(stat.star)
    #   stat* from MBB in a dataframe
    stat.star <- do.call("cbind",lapply(boot.results.all,FUN=function(x){x[,2]}))
    stat.star <- as.data.frame(stat.star)

    #   ADDED
    #   compute SE(stat.obs) - from stat*
    sd.stat <- apply(stat.star,1,FUN=sd,na.rm=FALSE)
    #   empirical Bayes shrinkage
    beta.obs <- res.obs[,2]
    shrink.beta.obs <- suppressMessages(ash(beta.obs,sebetahat = sd.stat,mixcompdist = "normal"))
    shrink.beta.est <- shrink.beta.obs$result$PosteriorMean
    shrink.beta.sd <- shrink.beta.obs$result$PosteriorSD
    rm(boot.results.all)

    # #   add observed value of stat to a dataframe
    # stat.obs <- data.frame(stat.obs=res.obs[,1])
    #   add observed value of stat to a dataframe
    stat.obs <- data.frame(stat.obs=shrink.beta.est)

    #   compute SE(stat.obs) - from stat*
    sd.stat <- data.frame(sd.stat=apply(stat.star,1,FUN=sd,na.rm=FALSE))

    # #   compute T observed
    # T.obs <- data.frame(T.obs=stat.obs/sd.stat)
    #   compute T observed
    T.obs <- data.frame(T.obs=stat.obs/shrink.beta.sd)
    if(is.null(T.obs.full)){
        T.obs <- T.obs
    }else{
        T.obs <- as.data.frame(T.obs.full)
    }

    # #   compute   SE(stat*) from double MBB
    # stat.star.star <-  lapply(boot.results.bb,function(x){
    #     do.call("cbind",lapply(x,FUN=function(x){x[,1]}))
    # 
    # })
    # 
    #   compute   SE(stat*) from double MBB
    stat.star.star <-  lapply(boot.results.bb,function(x){
            do.call("cbind",lapply(x,FUN=function(x){x[,2]}))
            
    })

    # sd.stat.star <- do.call("cbind",lapply(stat.star.star,function(x){data.frame(sd.stat.star=apply(x,1,FUN=sd,na.rm=TRUE))}))
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
    
    
    # #   compute stat*-stat.obs
    # T.d.star <- data.frame(apply(stat.star,2,function(x){x-stat.obs}))
    # 
    # T.t.star <- T.d.star/sd.stat.star
    #   compute stat*-stat.obs
    T.num.star <- data.frame(apply(stat.star,2,function(x){x-stat.obs}))
    
    T.star <- T.num.star/sd.stat.star
    
    #df.stat <- bind_cols(T.t.star,T.obs=T.obs[,1])
    df.stat <- bind_cols(T.star,T.obs=T.obs[,1])

    K.val <- apply(df.stat,1,function(x){sum(abs(x[1:R])>=abs(x[(R+1)]))/R})

    rt <- list(K.val, T.obs)

    rm(list=c("stat.star","stat.obs","sd.stat","stat.star.star","sd.stat.star","T.d.star","df.stat"))
    gc(reset = TRUE)

    return(rt)

}

