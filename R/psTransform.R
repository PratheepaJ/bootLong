#' psTransform
#'
#' This function ctakes the \code{phyloseq} object, apply
#' the arcsinh transformation, computes the mean-variance relationship for arcsinh transformed counts, assign a weight for each observation based on the predicted variance.
#'
#' @param ps \code{phyloseq} object
#' @param span the smoother span. see \code{lowess}
#'
#' @return \code{phyloseq} object with transformed \code{otu_table}
#' @export
#' @import "joineR"
#' @importFrom MASS glm.nb
psTransform = function(ps, main_factor, span = 0.5){

    inv_arcs = function(x) {
        y = 0.5*exp(-x)*(exp(2*x)-1)
        return(y)
    }

    if(dim(otu_table(ps))[2]!=nsamples(ps)){
        otu_table(ps)=t(otu_table(ps))
    }

    ot = otu_table(ps) %>% data.frame %>% as.matrix

    if(!isTRUE(all(ot == floor(ot)))){
        stop("otu_table entries must be integer")
    }

    #   variance stabilization for negative binomial distribution accounting for library sizes.
    geo_mean = function(x){
        val = exp(mean(log(x[x>0])))
    }

    geom_mean_row = apply(ot, 1, FUN = geo_mean)

    median_ratios = function(x, geom_mean_row){
        rat = x/geom_mean_row
        median_rat <- median(rat)
        return(median_rat)
    }

    sj = apply(ot, 2, FUN = median_ratios, geom_mean_row = geom_mean_row)

    #   matrix/vector or matrix*vector - rowwise division or rowwise multiplication
    ot_trans = t(asinh(t(ot)*sj))
    lib.size = colSums(ot)*sj
    #   compute weights and residuals
    samdf = sample_data(ps) %>% data.frame
    des = as.formula(paste("otu","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))

    residulas_response <- function(ind, samdf, ot, sj, des){
        otu = as.numeric(ot[ind,])
        sj = as.numeric(sj)
        dff = mutate(samdf, otu=otu, sj=sj)
        glmft = MASS::glm.nb(des, data = dff, method = "glm.fit", link = arcsinhLink())
        res_residuals = resid(glmft,"response")
        return(res_residuals)
    }


    lapply(seq_len(dim(ot_trans)[1]), function(x){fit <- lm(ot_trans[x,]~1); fit$coefficients %*% t(fit$design)})

    #   compute weight for each observation using the mean-variance relationship
    mean_asv = apply(ot_trans, 1, FUN = mean)
    sd_asv = sqrt(apply(ot_trans, 1, FUN = sd))

    allzero = (rowSums(ot) == 0)

    if (any(allzero)) {
        mean_asv = mean_asv[!allzero]
        sd_asv = sd_asv[!allzero]
    }

    l <- lowess(mean_asv, sd_asv, f = span)

    f <- approxfun(l, rule = 2)




        ot = as.matrix(otu_table(ps))
        samd <- data.frame(sample_data(ps))
        anno <- data.frame(tax_table(ps))
        dgeList <- edgeR::DGEList(counts=ot, genes=anno, samples = samd)

        des <- as.formula(paste("~", paste(main_factor, collapse="+")))
        mm <- model.matrix(des,data=samd)
        des2 <- as.formula(paste("otu","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))

        geo.mean.row <- apply((ot+1),1,function(x){exp(sum(log(x))/length(x))})
        sj <- apply((ot+1),2,function(x){median(x/geo.mean.row)})

        v <- arcsinhTransform(counts=dgeList, design=mm,plot = F)
        we <- v$weights

        nt <- as.list(seq(1,ntaxa(ps)))

        com.res <- function(ind,samd,ot,sj,we,des2){
                otu <- as.numeric(ot[ind,])
                sj <- as.numeric(sj)
                otu_sj <-otu*sj
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

