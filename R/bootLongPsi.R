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
    res.obs <- compute_stat(ps=ps,factors=factors)

    boot.results <- list()

    # ps.boot <- lapply(seq_len(R),FUN=function(i){
    #     bootLongPhyloseq(ps=ps,b=b,time=time)[[1]]
    # })
    #
    # df.boot <- lapply(ps.boot,FUN=function(pb){
    #     compute_stat(pb,factors)
    # })
    #
    # #ps.boot <- mapply(append, ps.boot, RR, SIMPLIFY = FALSE)
    #
    # boot.results.bb <- lapply(ps.boot,FUN=function(k){
    #     ps.boot.bb <- bootLongPhyloseq(k,b,time)
    #     return(ps.boot.bb)
    #
    # })

    # boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
    #     ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)
    #     ps.boot.bb <- ps.boot.bb[[1]]
    #     df.boot.bb <- compute_stat(ps.boot.bb,factors)
    #     rm(ps.boot.bb)
    #     return(df.boot.bb)
    # })


    boot.results <- lapply(seq_len(R),FUN=function(i){
        ps.boot <- bootLongPhyloseq(ps,b,time)[[1]]


        df.boot <- compute_stat(ps.boot,factors)

        #   double MBB
        boot.results.bb <- lapply(seq_len(RR),FUN=function(j){
            ps.boot.bb <- bootLongPhyloseq(ps.boot,b,time)[[1]]

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
    if(is.null(T.obs.full)){
        T.obs <- T.obs
    }else{
        T.obs <- as.data.frame(T.obs.full)
    }

    #   compute   SE(stat*) from double MBB
    stat.star.star <-  lapply(boot.results.bb,function(x){
        do.call("cbind",lapply(x,FUN=function(x){x[,1]}))

    })

    sd.stat.star <- do.call("cbind",lapply(stat.star.star,function(x){data.frame(sd.stat.star=apply(x,1,FUN=sd,na.rm=TRUE))}))

    #   compute stat*-stat.obs
    T.d.star <- data.frame(apply(stat.star,2,function(x){x-stat.obs}))

    T.t.star <- T.d.star/sd.stat.star

    df.stat <- bind_cols(T.t.star,T.obs=T.obs[,1])

    K.val <- apply(df.stat,1,function(x){sum(abs(x[1:R])>=abs(x[(R+1)]))/R})

    rt <- list(K.val, T.obs)

    rm(list=c("stat.star","stat.obs","sd.stat","stat.star.star","sd.stat.star","T.d.star","df.stat"))
    gc(reset = TRUE)

    return(rt)

}

