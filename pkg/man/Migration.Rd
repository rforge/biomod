\name{Migration}
\alias{Migration}

\title{ Using migration dynamics in species modelling }
\description{
  This function allows the inclusion of a very simple migration process when projecting species distributions
  into the future. Two datasets are needed : the current species distributions and the future (assuming by
  default unlimited migration).
}
\usage{
Migration(CurrentPred = NULL, FutureProj = NULL, X = NULL, Y = NULL, MaxMigr = NULL, Pred.Save = NULL)
}

\arguments{
  \item{CurrentPred}{ a matrix of n (number of species) columns giving the current state of the species in binary  }
  \item{FutureProj}{ a matrix of n (number of species) columns giving the future predictions of the species in binary }
  \item{X}{ , }
  \item{Y}{ X and Y : the latitude and longitude giving the coordinates of the data points given in CurrentPred and FutureProj }
  \item{MaxMigr}{ a single value or vector, of length the number of species modelled, giving the maximum distance the species 
  could migrate according to the time slice modelled (see the details section). }
  \item{Pred.Save}{ the name where the projections using limited migration will be stored }
}
\details{
  To be able to use this function, both current and future datasets must be ordered in the
 same way and have the same coordinates and resolution.

 Note that the rate of migration should be given in degrees. For instance a species could only have a
 maximum migration of 1 minute (1.6km) per decade. If we project its distribution in 50 years, the rate is
 1 x 0.16667 x 5 (where 0.1667 is the conversion from minute to degree).
 For projection in 2080 : Rate = 1 x 0.16667 x 8.
 For a maximum rate of 3 minutes per decade (4.8km) in 2080: Rate = 3 x 0.16667 x 8
}
\value{ 



}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{Biomod.RangeSize}}, \code{\link{Biomod.Turnover}} }
\examples{


}

\keyword{ models }
\keyword{ optimize }
