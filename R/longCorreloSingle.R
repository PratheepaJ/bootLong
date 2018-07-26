#' longCorreloSingle
#'
#' This function plots the correlogram given the taxon index.
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

longCorreloSingle <- function(ps,factors,time,taxon,taxlevel="Species"){
        taxon_name <- tax_table(ps)[taxon,taxlevel]
        df.taxa <- data.frame(sample_data(ps),otu=as.numeric(t(otu_table(ps)[taxon,])))
        names(df.taxa)[names(df.taxa)==factors] <- "Group"
        names(df.taxa)[names(df.taxa)==time] <- "Time"

        df.taxa.sep <- split(df.taxa,df.taxa$Group)

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

        return(p)
}
