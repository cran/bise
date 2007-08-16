\name{noMix}

\alias{noMix}

\title{Function to dermine the number of underlying distributions}

\description{Uses Gaussian Mixture Modelling based on optimizer and AIC to dermine best model fit}

\usage{noMix(obsday, count)}

\arguments{
  \item{obsday}{Julian day of observed budburst, e.g. 1st of May = 120}
  \item{count}{The number of observed budbursts on that day} }

\details{optimizer, gaussian Mixtures Modelling}

\value{predicted mixtures distribution, AIC, number of Gaussian Mixtures}

\author{Daniel Doktor}

\keyword{models}

\examples{\dontrun{no.mix(obsday, count)}}
