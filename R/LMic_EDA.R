#' ps_trans
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#'
#' @return \code{phyloseq} object with transformed \code{otu_table}
#' @export
#' @import "joineR"
ps_trans <- function(ps,factors){
    ot <- as.matrix(otu_table(ps))
    anno <- data.frame(tax_table(ps))
    samdf <- data.frame(sample_data(ps))
    dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samdf)

    #   setting up the model
    des <- as.formula(paste("~", paste(factors, collapse="+")))
    mm <- model.matrix(des,data=samdf)
    pse <- ps
    otu_table(pse) <- otu_table(ps)+1
    pDE <- suppressMessages(phyloseq_to_deseq2(pse,design=des))
    rm(pse)
    pDE <- DESeq2::estimateSizeFactors(pDE)
    sj <- DESeq2::sizeFactors(pDE)
    #   NOTE: if we want to use gene specific size factor, then
    #sij <- normalizationFactors(pDE) # matrix

    v <- voom_arcsine(counts=dgeList, design=mm, lib.size=sj,plot = F)

    transformed.ot <- data.frame(v$E)
    colnames(transformed.ot) <- sample_names(ps)

    pstr <- merge_phyloseq(otu_table(transformed.ot,taxa_are_rows = TRUE),sample_data(ps),tax_table(ps))
    return(pstr)
}


#' LMic_vario
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param taxon numeric, index of taxon to get the variogram
#' @param time character, time variable at repeated observations.
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot} object of variogram
#' @export
#'
LMic_vario <- function(ps,factors,time,taxon,point=FALSE,taxlevel="Species"){
    #ps.tr <- ps_trans(ps,factors="Preterm")
    taxon_name <- tax_table(ps)[taxon,taxlevel]
    df.taxa <- data.frame(sample_data(ps),otu=as.numeric(t(otu_table(ps)[taxon,])))
    names(df.taxa)[names(df.taxa)==factors] <- "Group"
    names(df.taxa)[names(df.taxa)==time] <- "Time"
    #   split by the levels of a factor
    df.taxa.sep <- split(df.taxa,df.taxa$Group)
    #   compute variogram for each level
    res.sep <- lapply(df.taxa.sep,function(m){
        joineR::variogram(m$SubjectID,time=m$Time,Y=m$otu)
    })

    df.sep <- lapply(res.sep,function(m){data.frame(m$svar)})
    Group <- names(df.sep)
    sig2 <- lapply(res.sep,function(m){m$sigma2})
    for(i in 1:length(df.sep)){
        df.sep[[i]]$Group <- Group[i]
        df.sep[[i]]$sigma2 <- sig2[[i]]
    }

    # sig2 <- lapply(res.sep,function(m){m$sigma2})
    sig2 <- unlist(sig2)
    df.sep <- do.call("rbind",df.sep)
    df.sep$Group <- as.factor(df.sep$Group)
    yl <- max(sig2)
    if(point==TRUE){
        p <- ggplot(df.sep)+
            geom_point(aes(x=vt,y=vv,col=Group,group=Group),size=1)+
            geom_smooth(aes(x=vt,y=vv,col=Group,group=Group),method = "loess")+
            geom_hline(aes(yintercept = sigma2,group=Group,col=Group),linetype=2,size=1)+
            xlab("u")+
            ylab("v(u)")+
            ggtitle(taxon_name)+
            theme(plot.title = element_text(hjust = 0.5,size=8),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),legend.text = element_text(size=8),legend.title = element_text(size=8))
        #     stat_summary(aes(x=vt,y=vv,col=Group,group=Group),fun.y="mean",geom="smooth")
    }else{
        p <- ggplot(df.sep)+
            geom_point(aes(x=vt,y=vv,col=Group,group=Group),size=1,col="white")+
            geom_smooth(aes(x=vt,y=vv,col=Group,group=Group),method = "loess")+
            geom_hline(aes(yintercept = sigma2,group=Group,col=Group),linetype=2,size=1)+
            xlab("u")+
            ylab("v(u)")+
            ggtitle(taxon_name)+
            theme(plot.title = element_text(hjust = 0.5,size=8),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),legend.text = element_text(size=8),legend.title = element_text(size=8))
        #     stat_summary(aes(x=vt,y=vv,col=Group,group=Group),fun.y="mean",geom="smooth")
    }

    return(p)
}



#' LMic_vario_multaxa
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#' @param starttaxa numeric, starting index of taxon: taxa will be ordered inside the function
#' @param endtaxa numeric, last index of taxon: taxa will be ordered inside the function
#' @param time character, time variable at repeated observations.
#' @param taxlevel character, taxonomy level to put the title
#'
#' @return \code{ggplot} object of variogram for starttaxa:endtaxa taxa
#' @export
#'
LMic_vario_multaxa <- function(ps,factors,time,starttaxa=1,endtaxa=4,point=FALSE,taxlevel="Species"){
    # ps.tr <- ps_trans(ps,factors=factors)
    #   order tax by taxa sum
    #taxa_order <- names(sort(taxa_sums(ps),decreasing = TRUE))
    taxa_order <- sort(taxa_sums(ps),decreasing = T)
    ind <- which(taxa_names(ps)%in%names(taxa_order)[starttaxa:endtaxa])
    taxa <- as.list(ind)
    p.all <- lapply(taxa,function(x){LMic_vario(ps=ps,factors=factors,time=time,taxon=x,point=point,taxlevel = taxlevel)})
    return(p.all)
}
