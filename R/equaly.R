.packageName <- "bise"

.First.lib <- function(lib, pkg) {
              library.dynam("bise", pkg, lib)
              print("test version of bise loaded")
}

equaly <- function(obs1,year1,obs2,year2) {

    n1 <- length(obs1)
    n2 <- length(obs2)
    count <- 0
    obsday1 <- vector(mode="numeric", length=51)
    obsday2 <- vector(mode="numeric", length=51)
    res <- .Fortran("equaly", n1=as.integer(n1), n2=as.integer(n2), obs1=as.integer(obs1), year1=as.integer(year1), obs2=as.integer(obs2), year2=as.integer(year2), obsday1=as.integer(obsday1), obsday2=as.integer(obsday2), count=as.integer(count), PACKAGE="bise")
    attach(res)
    return(res)
}
