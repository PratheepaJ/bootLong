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

    inv_arcs = function(x){
        y = 0.5*exp(-x)*(exp(2*x)-1)
        return(y)
    }

    if(dim(otu_table(ps))[2] != nsamples(ps)){
        otu_table(ps) = t(otu_table(ps))
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

    # library size normalization factor
    sj = apply(ot, 2, FUN = median_ratios, geom_mean_row = geom_mean_row)

    #   normalized library size
    #lib.size = colSums(ot)*sj

    #   matrix/vector or matrix*vector - rowwise division or rowwise multiplication
    ot_trans = t(asinh(t(ot)*sj))

    ##   compute residuals and weights
    samdf = sample_data(ps) %>% data.frame
    des = as.formula(paste("otu","~", paste(main_factor, collapse="+"),"+","offset(arcsinhLink()$linkfun(sj))"))

    response_residulas_fitted <- function(ind, samdf, ot, sj, des){
        otu = as.numeric(ot[ind,])
        sj = as.numeric(sj)
        dff = mutate(samdf, otu=otu, sj=sj)
        glmft = MASS::glm.nb(des, data = dff, method = "glm.fit", link = arcsinhLink())
        res_residuals = resid(glmft,"response")
        res_fitted = fitted(glmft, "response")
        rt = list(res_residuals, res_fitted)
        names(rt) <- c("response_residuals", "response_fitted")
        return(rt)
    }

    resi_fitted = lapply(seq_len(ntaxa(ps)),function(x){
        response_residulas_fitted(x,samdf = samdf,ot = (ot+1), sj = sj,des = des)
        })

    resi = lapply(resi_fitted, "[[", 1)
    resi = data.frame(do.call("rbind",resi))
    colnames(resi) = sample_names(ps)
    rownames(resi) = taxa_names(ps)

    resi_trans = t(asinh(t(resi)*sj))

    #   compute weight for each observation using the mean-variance relationship of asinh transformed observations accounting for library size
    mean_asv = apply(ot_trans, 1, FUN = mean)
    sd_asv = sqrt(apply(ot_trans, 1, FUN = sd))

    allzero = (rowSums(ot) == 0)

    if (any(allzero)) {
        mean_asv = mean_asv[!allzero]
        sd_asv = sd_asv[!allzero]
    }

    l = lowess(mean_asv, sd_asv, f = span)

    f = approxfun(l, rule = 2)

    fited = lapply(resi_fitted, "[[", 2)
    fited = data.frame(do.call("rbind",fited))
    colnames(fited) = sample_names(ps)
    rownames(fited) = taxa_names(ps)

    fited_trans = t(asinh(t(fited)*sj))

    w = 1/f(fited_trans)^4
    dim(w) = dim(fited_trans)

    ps_resid_asinh = phyloseq(otu_table(resi_trans,taxa_are_rows = T),
                       sample_data(ps),
                       tax_table(ps))

    ps_ot_asinh = phyloseq(otu_table(ot_trans,taxa_are_rows = TRUE),
                        sample_data(ps),
                        tax_table(ps))

    rt <- list(ps_ot_asinh, ps_resid_asinh)
    names(rt) <- c("asinh_transformed_counts", "asinh_transformed_residulas")
    return(rt)
}

