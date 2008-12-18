\name{Biomod.RangeSize}
\alias{Biomod.RangeSize}

\title{ Analysis of the range size changes }
\description{
  This function allows to estimate the proportion and relative number of pixels (or habitat) lost, gained
and stable for the time slice considered in species-climate modelling under future scenarios.
}
\usage{
Biomod.RangeSize(CurrentPred = NULL, FutureProj = NULL, SpChange.Save = NULL)
}

\arguments{
  \item{CurrentPred}{ a matrix of n (number of species) columns giving the current state of the species in binary  }
  \item{FutureProj}{ a matrix of n (number of species) columns giving the future projections of the species in binary according to a future scenario}
  \item{SpChange.Save}{ the name given to the new object storing the results }
}
\details{
  Note that this function is only relevant if you make projections on the same area with the same resolution and
  for a different time slice (past or future) than for the current data.
}
\value{
  A list of two datasets is created: Compt.By.Species and Diff.By.Pixel

  \item{Compt.By.Species}{ stores the summary of range change for each species (sorted by rows).
  The first four columns are absolute values whereas the next 3 ones are relative values: 
  \item{Disa}{ represents the number of pixels predicted to be lost by the given species.} 
  \item{Stable0}{ is the number of pixels which are not currently occupied by the given species and not predicted to be.} 
  \item{Stable1}{ represents the number of pixels currently occupied by the given species, and predicted to remain occupied into the future.}
  \item{Gain}{ represents the number of pixels which are currently not occupied by the given species but predicted to be into the future.}
  \item{PercLoss}{corresponds to the percentage of currently occupied sites to be lost (Disa/(Disa+Stable1)}
  \item{PercGain}{corresponds to the percentage of new sites considering the species' current distribution size (Gain/(Disa+Stable1). 
  For example, if the there are 30 sites currently occupied and 15 new sites are projected to be occupied in future, it makes PercGain=+50(\%).}
  \item{SpeciesRangeChange}{ it is the overall projection outcome, equal to PercGain-PercLoss. It does not assess for any migration shifts as it 
  strictly compares the range sizes between current and future states.}
  \item{CurrentRangeSize}{ represents the modelled current range size (number of pixels occupied) of the given species.}
  \item{FutureRangeSize0Disp}{ represents the future modelled range size assuming no migration of the given species.}
  \item{FutureRangeSize1Disp}{ represents the future modelled range size assuming migration of the given species (depending on the datasets given in input, if Migration has been used or not).}}
  
  \item{Diff.By.Pixel}{the summary of range change for each species (sorted by columns and with the pixel
  in rows). For each species, a pixel could have four different values :
  
  -2 if the given pixel is predicted to be lost by the species. 
  -1 if the given pixel is predicted to be stable for the species.
   0 is the given pixel was not occupied, and will not be in the future.
   1 if the given pixel was not occupied, and is predicted to be into the future.
  
  This table could be easily plotted into GIS software in order to represent the pattern of change for the selected species (or even with the level.plot() function).}

}

\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Biomod.Turnover}}, \code{\link{Migration}}}
\examples{



}

\keyword{ IO }
