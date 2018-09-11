#' link function to be used in the negative binomial regression fitting
#'
#' @return A list of linkFun (link function), varFun(variance function),
#' InvLink (inverse of link function), and InvLinkDeriv (derivative of inverse of the link function)
#' @export

arcsinhlstLink <- function(){
        LinkFun = function(y){log(y + sqrt(y ^ 2 + 1))}
        VarFun = function(y){y*(1+y/theta)}
        InvLink = function(eta){0.5*exp(-eta)*(exp(2*eta)-1)}
        InvLinkDeriv = function(eta){.5*(exp(eta)+exp(-eta))}
        FunList = list(LinkFun,VarFun,InvLink,InvLinkDeriv)
        return(FunList)
}
