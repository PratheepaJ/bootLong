#' LMic_correloram
#'
#' Autocorrelation plot
#'
#' @param ps \code{phyloseq} object
#' @param factors factor variable in the sample data of ps.
#' @param time character, time variable at repeated observations.
#' @param taxon numeric, index of taxon to get the variogram
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot2} object of correlogram for one taxon according to the taxon index.
#' @export
#'

LMic_correloram <- function(ps,factors,time,taxon,taxlevel="Species"){
        taxon_name <- tax_table(ps)[taxon,taxlevel]
        df.taxa <- data.frame(sample_data(ps),otu=as.numeric(t(otu_table(ps)[taxon,])))
        names(df.taxa)[names(df.taxa)==factors] <- "Group"
        names(df.taxa)[names(df.taxa)==time] <- "Time"
        #   split by the levels of a factor
        df.taxa.sep <- split(df.taxa,df.taxa$Group)
        #   compute correlogram for each level
        res.sep <- lapply(df.taxa.sep,function(m){
            m$Time <- as.factor(m$Time)
            m <- m %>% group_by(Time) %>% summarise(meant=mean(otu))
            acf.res <- acf(m$meant,plot=F)
            return(acf.res)
        })

        df.sep <- lapply(res.sep,function(m){
            data.frame(lag=m$lag,acf=m$acf)
            })

        Group <- names(df.sep)
        for(i in 1:length(df.sep)){
            df.sep[[i]]$Group <- Group[i]
        }

        df.sep <- do.call("rbind",df.sep)
        df.sep$Group <- as.factor(df.sep$Group)
        df.sep$lag <- as.factor(df.sep$lag)

        p <- ggplot(df.sep,aes(x=lag,y=acf,group=Group,fill=Group))+
            geom_bar(stat="identity", width=.5, position = "dodge")+
            xlab("h")+
            ylab(expression(paste(rho,"(h)")))+
            ggtitle(taxon_name)+
            theme(plot.title = element_text(hjust = 0.5,size=8),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),legend.text = element_text(size=8),legend.title = element_text(size=8))
        # geom_hline(aes(yintercept = 0)) +
        #     geom_segment(mapping = aes(xend = lag, yend = 0),size=.8)+
        #     scale_alpha_discrete(range = c(0.6, 0.4),guide = 'none')+
        return(p)
}


#' LMic_correloram_multaxa
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param time character, time variable at repeated observations.
#' @param starttaxa numeric, first index of taxa ordered by total reads
#' @param endtaxa numeric, the last index of taxa ordered by total reads
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot2} object of correlogram for multiple taxa.
#' @export
#'
LMic_correloram_multaxa <- function(ps,factors,time,starttaxa=1,endtaxa=4,taxlevel="Species"){
    taxa_order <- sort(taxa_sums(ps),decreasing = T)
    ind <- which(taxa_names(ps)%in%names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa,function(x){LMic_correloram(ps=ps,factors=factors,time=time,taxon=x,taxlevel = taxlevel)})
    return(p.all)
}


#' LMic_lagged
#'
#' @param ps ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param time character, time variable at repeated observations.
#' @param taxon numeric, index of taxon to get the variogram
#' @param lags numeric (vector) how many lags to plot
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot2} object of lagged plot for speciec taxon.
#' @export
#'
LMic_lagged <- function(ps,factors,time,taxon,lags,taxlevel="Species"){
    taxon_name <- tax_table(ps)[taxon,taxlevel]
    df.taxa <- data.frame(sample_data(ps),otu=as.numeric(t(otu_table(ps)[taxon,])))
    names(df.taxa)[names(df.taxa)==factors] <- "Group"
    names(df.taxa)[names(df.taxa)==time] <- "Time"
    #   split by the levels of a factor
    df.taxa.sep <- split(df.taxa,df.taxa$Group)
    #   compute correlogram for each level
    res.sep <- lapply(df.taxa.sep,function(m){
        m$Time <- as.factor(m$Time)
        m <- m %>% group_by(Time) %>% summarise(meant=mean(otu))
        xax <- m$meant[-c(1:lags)]
        yax <- m$meant[1:(dim(m)[1]-lags)]
        df.res <- data.frame(xax=xax,yax=yax)
        return(df.res)
    })

    df.sep <- res.sep

    Group <- names(df.sep)
    for(i in 1:length(df.sep)){
        df.sep[[i]]$Group <- Group[i]
    }

    df.sep <- do.call("rbind",df.sep)
    df.sep$Group <- as.factor(df.sep$Group)

    p <- ggplot(df.sep,aes(x=xax,y=yax,col=Group,group=Group))+
        geom_point()+
        ggtitle(paste("lag",lags))+
        xlab("t")+
        ylab(paste("t+",lags))+
        theme(plot.title = element_text(hjust = 0.5,size=8),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),legend.text = element_text(size=8),legend.title = element_text(size=8))
    # scale_alpha_discrete(range = c(0.6, 0.4))+
    return(p)
}

