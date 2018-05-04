arcsinhlstLink <- function(){
        LinkFun <- function(y) log(y + sqrt(y ^ 2 + 1))
        VarFun <- function(y){y+theta*y^2}
        InvLink <- function(eta)  0.5*exp(-eta)*(exp(2*eta)-1)
        InvLinkDeriv <- function(eta) {.5*(exp(eta)+exp(-eta))}
        FunList <- list(LinkFun,VarFun,InvLink,InvLinkDeriv)
        return(FunList)
}
