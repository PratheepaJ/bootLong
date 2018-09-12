#' bootLongMSEPsi
#'
#' Compute MSE in computing K = two-sided probability with different block sizes
#'
#' @param Wj subset of repeated observations for j-th subject
#' @param qj number of repeated observations for j-th subject
#' @param Khat.obs second element of the output of \code{bootLongPsi} evaluated using the initial block length and full data
#' @inheritParams bootLongPsi
#'
#' @return list of MSE computed with block length b, Khat is all subsamples, Khat with initial block length
#'
#'
#' @export
bootLongMSEPsi <- function(ps,
                           main_factor,
                           time_var,
                           subjectID_var,
                           b,
                           R,
                           RR,
                           qj,
                           Wj,
                           Khat.obs = NULL,
                           T.obs.full = NULL){

        doParallel::registerDoParallel(parallel::detectCores())
        BiocParallel::register(BiocParallel::DoparParam())

        if(is.null(Khat.obs)){
            stop("User needs to run bootLongPsi() function with an initial block length ")
            }
        if(is.null(T.obs.full)){
            stop("User needs to provide observed test statistic")
            }


        samdf <- sample_data(ps) %>% data.frame
        if(!is.numeric(samdf[,time_var])){
            samdf[,time_var] <- as.numeric(samdf[,time_var])
            }
        g <- samdf[,subjectID_var]
        samdf <- split(samdf,g)
        num.sub.sam <- max(qj)-max(Wj)+1

        samdf.q.W <- mapply(samdf,as.list(qj), as.list(Wj), FUN=list,SIMPLIFY = FALSE)

        samdf.q.W.or <- lapply(samdf.q.W,function(x){
            x[[1]] <- dplyr::arrange_(x[[1]],time_var)
            return(x)
        })


        if(num.sub.sam < 5){
            stop("decrease omega")
            }

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
            subsam.id <- as.character(subsam.id)
            ps.sub[[i]] <- prune_samples(subsam.id,ps)

        }

        Khat <- BiocParallel::bplapply(ps.sub,function(x){
            k.hat <- bootLongPsi(x,
                                 main_factor=main_factor,
                                 time_var=time_var,
                                 subjectID_var = subjectID_var,
                                 b=b,
                                 R=R,
                                 RR=RR,
                                 T.obs.full=T.obs.full)
            k.hat <- k.hat[[1]]
            return(k.hat)
        })

        rm(ps)
        rm(samdf)
        rm(samdf.q.W)
        rm(samdf.q.W.or)
        rm(ps.sub)

        Khat.squared.diff <- lapply(Khat, FUN=function(w){
            (w-Khat.obs)^2
            })

        Khat.squared.diff.df <- do.call("cbind", Khat.squared.diff)

        MSE_i <- apply(Khat.squared.diff.df, 1 , FUN=function(x){
            mean(x)
            })

        rm(Khat.squared.diff.df)

        rt <- list(MSE_i=MSE_i,Khat=Khat,Khat.obs=Khat.obs)

        gc(reset = TRUE)

        return(rt)
}
