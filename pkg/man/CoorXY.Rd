\name{CoorXY}
\alias{CoorXY}
\docType{data}
\title{ Coordinates for plotting the BIOMOD example maps}
\description{
  This dataset contains the coordinates of the data points for the purpose of making plots
  of the different datasets produced with BIOMOD.
}
\usage{data(CoorXY)}
\format{
  A data frame with 2264 observations on the following 2 variables.
  \describe{
    \item{\code{X}}{and}
    \item{\code{Y}}{ : numeric vectors giving respectively the latitude and longitude of the data points}
  }
}
\details{
  These data correspond to a 50x50 km grid for Europe 
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  
}
\examples{
data(CoorXY)
plot(CoorXY[,1]~CoorXY[,2])
}
\keyword{datasets}
