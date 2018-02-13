asinh_link <- function() {
    ## link
    linkfun <- function(y) log(y + sqrt(y ^ 2 + 1))
    ## inverse link
    linkinv <- function(eta)  0.5*exp(-eta)*(exp(2*eta)-1)
    ## derivative of invlink wrt eta
    mu.eta <- function(eta) { .5*(exp(eta)+exp(-eta))}
    valideta <- function(eta) TRUE
    link <- "asinh"
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta,
                   name = link),
              class = "link-glm")
}

# #   Basic checks
# vv <- asinh_link()
# vv$linkfun(vv$linkinv(27))  ## check invertibility
# library("numDeriv")
# all.equal(grad(vv$linkinv,2),vv$mu.eta(2))
#
# #   Example
# set.seed(101)
# n <- 1000
# x <- runif(n)
# sh <- 2
# y <- rgamma(n,scale=vv$linkinv(2+3*x)/sh,shape=sh)
# glm(y~x,family=Gamma(link=vv))



#log(E(Y)-1))

# vlog <- function() {
#     ## link
#     linkfun <- function(y) log(exp(y)-1)
#     ## inverse link
#     linkinv <- function(eta)  log(exp(eta)+1)
#     ## derivative of invlink wrt eta
#     mu.eta <- function(eta) { 1/(exp(-eta) + 1) }
#     valideta <- function(eta) TRUE
#     link <- "log(exp(y)-1)"
#     structure(list(linkfun = linkfun, linkinv = linkinv,
#                    mu.eta = mu.eta, valideta = valideta,
#                    name = link),
#               class = "link-glm")
# }
#
# #   Basic checks
# vv <- vlog()
# vv$linkfun(vv$linkinv(27))  ## check invertibility
# library("numDeriv")
# all.equal(grad(vv$linkinv,2),vv$mu.eta(2))
#
# #   Example
# set.seed(101)
# n <- 1000
# x <- runif(n)
# sh <- 2
# y <- rgamma(n,scale=vv$linkinv(2+3*x)/sh,shape=sh)
# glm(y~x,family=Gamma(link=vv))
