.packageName <- "bise"

.First.lib <- function(lib, pkg) {
              library.dynam("bise", pkg, lib)
	      print("test version of bise loaded")
}

# Computes the altitudinal gradient of budburst
# input: vectors containing latitude, longitude, altitude
# and julian day of budburst, e.g. First of May = 120

altGrad <- function(x,y,alt,obsday) {
    n <- length(id)
    delBB <- 0.0
    res <- .Fortran("altdep", x=as.double(x), y=as.double(y), alt=as.integer(alt), obsday=as.integer(obsday), n=as.integer(n), delBB=as.double(delBB) , PACKAGE="bise")
    attach(res)
    return(delBB)
}
