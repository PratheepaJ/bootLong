#' ComputeMSE
#'
#' Compute MSE in computing K = two-sided probability with different block sizes
#'
#' @param ps Observed \code{phyloseq} class object.
#' @param Wj subset of repeated observations for j-th subject
#' @param qj number of repeated observations for j-th subject
#' @param b numeric, block lenght to account for dependence within-subject.
#' @param R Number of block bootstrap realization.
#' @param RR Number of double block bootstrap realization.
#' @param factors vector of factor variable(s) in the sample data of ps.
#' @param Khat.obs second element of the output of \code{ComputeK} evaluated using the initial block length and full data
#' @param time Time variable at repeated observations.
#' @param T.obs.full If observed statistic is already computed
#' @param compute_stat
#'
#' @return list of MSE computed with block length b, Khat is all subsamples, Khat with initial block length
#'
#'
#' @export
ComputeMSE <- function(ps,qj,Wj,b,R,RR,factors,time,Khat.obs=NULL,T.obs.full=NULL){

    if(is.null(Khat.obs)){stop("User needs to run ComputeK() function with initial block length ")}
    if(is.null(T.obs.full)){stop("User needs to provide observed test statistic")}

    #names(sample_data(ps))[names(sample_data(ps))==time] <- "Time"
    #   Create many ((max(qj)-max(Wj)+1) number of sub-seires) phyloseq with sub-series to compute MSE with block size "b"
    samdf <- data.frame(sample_data(ps))
    samdf <- split(samdf,samdf$SubjectID)
    num.sub.sam <- max(qj)-max(Wj)+1

    samdf.q.W <- mapply(samdf,as.list(qj),as.list(Wj),FUN=list,SIMPLIFY = FALSE)

    samdf.q.W.or <- lapply(samdf.q.W,function(x){
        if(!(is.numeric(x[[1]][,time]))){x[[1]][,time] <- as.numeric(x[[1]][,time])}
        x[[1]] <- dplyr::arrange_(x[[1]],time)
        return(x)
    })


    if(num.sub.sam<5){stop("decrease omega")}

    ps.sub <- list()
    for(i in 1:num.sub.sam){
        sub.sam.i <- lapply(samdf.q.W.or,function(x){
            xd <- x[[1]]
            W <- x[[3]]
            ss <- data.frame(dplyr::slice(x[[1]],i:(W+i-1)))
            return(ss)
        })
        sub.sam.i <- do.call("rbind",sub.sam.i)
        subsam.id <- sub.sam.i$SampleID
        ps.sub[[i]] <- subset_samples(ps,SampleID%in%subsam.id)
        #names(sample_data(ps.sub[[i]]))[names(sample_data(ps.sub[[i]]))=="Time"] <- time
    }

    Khat <- lapply(ps.sub,function(x){
        k.hat <- ComputeK(x,b=b,R=R,RR=RR,factors=factors,time=time,T.obs.full=T.obs.full)
        k.hat <- k.hat[[1]]
        return(k.hat)
    })

    # Khat <- lapply(subsam.num,function(i){
    #
    #     sub.samdf <- lapply(samdf,function(x){x[min(qj[which(names(qj)%in%unique(x$SubjectID))]-Wj[which(names(Wj)%in%unique(x$SubjectID))]+1,i):min(qj[which(names(qj)%in%unique(x$SubjectID))],Wj[which(names(Wj)%in%unique(x$SubjectID))]+i-1),]})
    #
    #     sub.samdf <- do.call("rbind",sub.samdf)
    #
    #     sub.SampleID <- as.character(sub.samdf$SampleID)
    #
    #     ps.m <- prune_samples(sub.SampleID,ps)
    #
    #     k.hat <- ComputeK(ps.m,b=b,R=R,RR=RR,factors=factors,time=time,T.obs.full=T.obs.full)
    #     k.hat <- k.hat[[1]]
    #
    #     rm(sub.samdf)
    #     rm(sub.SampleID)
    #     rm(ps.m)
    #
    #    return(k.hat)
    # })


    rm(ps)
    rm(samdf)
    rm(samdf.q.W)
    rm(samdf.q.W.or)
    rm(ps.sub)


    #deviation squared between two-tailed probability with block lengths lI and lC
    Khat.squared.diff <- lapply(Khat,FUN=function(w){(w-Khat.obs)^2})

    #   dataframe -
    Khat.squared.diff.df <- do.call("cbind",Khat.squared.diff)

    #   MSE over all possible subsamples
    MSE_i <- apply(Khat.squared.diff.df,1,FUN=function(x){mean(x)})

    rm(Khat.squared.diff.df)

    rt <- list(MSE_i=MSE_i,Khat=Khat,Khat.obs=Khat.obs)
    #   free up memory
    gc(reset = TRUE)

    #return(MSE_i)
    return(rt)
}
