\name{minimum}

\alias{minimum}

\title{minimum of Best Index Slope Extraction}

\description{
    Computes of the upper envelope of daily NDVI values over one year the minimum until max(ndvi).}

\usage{minimum(ndvi)}

\arguments{
  \item{ndvi}{Double(real), observed NDVI values between 0 and 1}
}

\details{to be followed}

\value{NDVI values extracted out of very noise satellite observations. These are regarded to be
       the 'TRUE' observations.}

\references{Viovy et. al}

\author{Daniel Doktor, Alberte Bondeau}

\keyword{utilities}

\examples{\dontrun{
    load("1991-400")
    ndvi <- result[,400]
    ndvi <- ifelse(ndvi < 0 | ndvi >= 1, 0, ndvi)
    ndvi <- ifelse(is.na(ndvi)==TRUE, 0, ndvi)
    minimum(ndvi)}}
