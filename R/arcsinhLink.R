#' arcsinhLink
#' arcsinh link function for GLM
#'
#' @return An object of class "link-glm"
#' @export
arcsinhLink <- function(){
    ## link
    linkfun <- function(y){log(y + sqrt(y ^ 2 + 1))}
    ## inverse link
    linkinv <- function(eta){pmax((0.5*exp(-eta)*(exp(2*eta)-1)), .Machine$double.eps)}
    ## derivative of invlink wrt eta
    mu.eta <- function(eta) {pmax((.5*(exp(eta)+exp(-eta))), .Machine$double.eps)}
    valideta <- function(eta){TRUE}
    link <- "asinh"
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta,
                   name = link),
              class = "link-glm")
}

