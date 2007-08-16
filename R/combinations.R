.packageName <- "bise"

.First.lib <- function(lib, pkg) {
              library.dynam("bise", pkg, lib)
              print("test version of bise loaded")
}


combinations <- function(id) {
    n <- length(id)
    ncomb <- (n^2-n) / 2
    stat1 <- vector(mode="numeric", length=ncomb)
    stat2 <- vector(mode="numeric", length=ncomb)
    res <- .Fortran("combinations", id=as.integer(id), n=as.integer(n), ncomb=as.integer(ncomb), stat1=as.integer(stat1), stat2=as.integer(stat2), PACKAGE="bise")
    attach(res)
    return(c(stat1,stat2))
}
