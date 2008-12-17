\name{Biomod.Turnover}
\alias{Biomod.Turnover}

\title{ Estimation of the Turnover for a given location and resolution}
\description{
  This function allows to estimate the species lost and gained by pixel for the time slice and scenario considered, 
  and the turnover induced by these (predicted) modifications.
}
\usage{
Biomod.Turnover(CurrentPred = NULL, FutureProj = NULL, Turnover.Save = NULL)
}

\arguments{
  \item{CurrentPred}{ a matrix of n (number of species) columns giving the current state of the species in binary  }
  \item{FutureProj}{ a matrix of n (number of species) columns giving the future predictions of the species in binary }
  \item{Turnover.Save}{ the name given to the new object storing the results }
}
\details{
  
}
\value{
  A 10 columns database is created (and saved under the name you gave it).
The first four columns are relative numbers: 
- Disa represents the number of species predicted to disappear from the given pixel. 
- Stable0 is the number of species which are currently not in the given pixel and not predicted to migrate into it. 
- Stable1 represents the number of species currently occurring in the given pixel, and predicted to remains into the future. 
- Gain represents the number of species which are currently absent but predicted to migrate in the given pixel.

PercLoss, PercGain and Turnover are the relative percentages estimated as the following:
- PercLoss = 100 x L/(SR)
- PercGain = 100 x G/(SR)
- Turnover = 100 x (L+G)/(SR+G)

Where SR is the current species richness, L the species lost and G those gained.

- CurrentSR represents the current modelled species richness in the given pixel.
- FutureSR0Disp represents the future modelled species richness assuming no migration of species
- FutureSR1Disp represents the future modelled species richness assuming migration (depending on
  the datasets given in input, if \code{\link{Migration}} has been used or not).

}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Biomod.RangeSize}}, \code{\link{Migration}} }
\examples{


}

\keyword{ IO }