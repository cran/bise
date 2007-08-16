\name{altGrad}

\alias{altGrad}

\title{Altitudinal gradient of budburst }

\description{
    Computes the altitudinal gradient of budburst. Due to on average decreasing temperatures with
    increasing altitude budburst occurs later in mountainous regions than in lowlands. This function
    computes the delayed budburst in days per 1 [m] increase in altitude.}

\usage{altGrad(x, y, alt, obsday)}

\arguments{
  \item{x}{Double(real), Longitude of observations in decimal degrees}
  \item{y}{Doube(real), Latitude of observations in decimal degrees}
  \item{alt}{Integer, Altitude of observations in [m]}
  \item{obsday}{Julian day of observed budburst, e.g. 1st of May = 120} }
  
\details{ gradient on altitude is \deqn{g_{ha} = g^{ha}} for example. }

\value{Delayed budburst in days per 1 [m] increase in altitude}

\references{Doktor D, Badeck FW, Hattermann F, Schaber J, McAllister M (2004) 
            'Analysis and modelling of spatially and temporally varying phenological
            phases'. Geostatistics for Environmental Applications. Proceedings of the Fifth
	    European Conference on Geostatistics for Environmental Applications}

\author{Daniel Doktor}

\keyword{utilities}

\examples{\dontrun{altGrad(x,y,alt,obsday)}}
