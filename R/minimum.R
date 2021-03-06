.packageName <- "bise"

.First.lib <- function(lib,pkg) {
    library.dynam("bise", pkg, lib)
}

minimum <- function (ndvi)
{
    ndvi <- ifelse(ndvi < 0 | ndvi >= 1 | is.na(ndvi) == TRUE,
        0, ndvi)
    ll <- 366
    db <- vector(mode = "integer", length = 366)
    nb <- ni <- nirm <- vector(mode = "numeric", length = 366)
    day <- 1:366
    SP <- 40
    CUT <- 2
    check <- sort(ndvi, decreasing = TRUE)
    if (check[30] < 0.1) {
        mini <- NA
    }
    else {
        mini <- 0
    }
    if (is.na(mini) == FALSE) {
        res <- .Call("bise", d0 = as.integer(day), n0 = as.numeric(ndvi), SP = as.integer(SP), CUT = as.numeric(CUT), ll = as.integer(ll), db = as.integer(db), nb = as.numeric(nb), ni = as.numeric(ni), nirm = as.numeric(nirm), PACKAGE="bise")
        mini <- min(res[1:order(res, decreasing = TRUE)[1] ] )
    }
    return(mini)
}
