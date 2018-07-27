#' psTransform
#' This function takes the \code{phyloseq} object and apply
#' the arcsinh transformation.
#'
#' @param ps \code{phyloseq} object
#' @param factors factor, variable in the sample data of ps.
#'
#' @return \code{phyloseq} object with transformed \code{otu_table}
#' @export
#' @import "joineR"
psTransform <- function(ps,factors){
        ot <- as.matrix(round(otu_table(ps),digits = 0))
        samd <- data.frame(sample_data(ps))
        anno <- data.frame(tax_table(ps))
        dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)

        des <- as.formula(paste("~", paste(factors, collapse="+")))
        mm <- model.matrix(des,data=samd)
        des2 <- as.formula(paste("otu","~", paste(factors, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))

        geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
        sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})

        v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)
        we <- v$weights

        nt <- as.list(seq(1,ntaxa(ps)))

        com.res <- function(ind,samd,ot,sj,we,des2){
                otu <- as.numeric(ot[ind,])
                sj <- as.numeric(sj)
                we <- as.numeric(we[ind,])
                dff <- samd
                dff <- cbind(samd,otu=otu,sj=sj,we=we)

                glmft.tx <- MASS::glm.nb(des2,data = dff,weights = we,method = "glm.fit",link = arcsinhLink())
                return(glmft.tx$residuals)
        }

        resi <- lapply(nt,function(x){com.res(x,samd,(ot+1),sj,we,des2)})
        resi <- data.frame(do.call("rbind",resi))

        colnames(resi) <- sample_names(ps)
        rownames(resi) <- taxa_names(ps)

        ps_res <- phyloseq(otu_table(resi,taxa_are_rows = T),sample_data(ps),tax_table(ps))

            v <- arcsinhTransform(counts=dgeList, design=mm, lib.size=sj,plot = F)

            transformed.ot <- data.frame(v$E)
            colnames(transformed.ot) <- sample_names(ps)

            pstr <- phyloseq(otu_table(transformed.ot,taxa_are_rows = TRUE),sample_data(ps),tax_table(ps))

            rt <- list(pstr,ps_res)
            return(rt)
}

