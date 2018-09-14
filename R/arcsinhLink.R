#' asinh link function for glm.nb
#'
#' @return An object of class "link-glm"
#' @export
arcsinhLink = function(){
    linkfun = function(y){log(y + sqrt(y ^ 2 + 1))}

    linkinv = function(eta){0.5*exp(-eta)*(exp(2*eta)-1)}

    mu.eta = function(eta) {.5*(exp(eta)+exp(-eta))}

    valideta = function(eta){TRUE}

    link = "asinh"

    structure(list(linkfun = linkfun,
                   linkinv = linkinv,
                   mu.eta = mu.eta,
                   valideta = valideta,
                   name = link),
              class = "link-glm")
}
