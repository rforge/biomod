\name{ProbDensFunc}
\alias{ProbDensFunc}

\title{ Probability Density Function }
\description{
  Using a variety of parameters in modelling will inevitably bring variability in predictions, especially when it comes to making future predictions.
  This function enables an overall viewing of the future predictions range per species and gives the likelihood of range shift estimations.
}
\usage{
ProbDensFunc(initial, projections, plothist = T, distrib = T, cvsn = T, groups = NULL, resolution = 5)
}

\arguments{
  \item{initial}{ a vector in a binary format (ones and zeros) representing the current distribution of a species which will be used as a reference for the range change calculations }
  \item{projections}{ a matrix grouping all the predictions where each column is a single prediction. Make sure you keep projections in the same order as the initial vector (line1=site1, line2=site2, etc.). }
  \item{plothist}{ set to TRUE to plot the range change histogram }
  \item{distrib}{ if true, the optimal way for condensing 50, 75, 90 and 95\% of the data will be calculated and shown on the histogram }
  \item{cvsn}{ stands for current vs new. If true, the range change calculations will be of two types: the percentage of cells currently 
  occupied by the species to be lost, and the relative percentage of cells currently unoccupied but projected to be, namely 'new' cells, compared to current surface range.}
  \item{groups}{ an option for ungrouping the projections enabling a separated visualisation of the prediction range per given group. A 
  matrix is expected where each column is a single prediction and each line is giving details of one parameter (See the examples section). }
  \item{resolution}{ the step used for classes of prediction in graphics. The default value is 5 }
}

\details{

  The future range changes are calculated as a percentage of the species' present state. For example, if a species currently occupies 100 cells and 
  is estimated by a model to cover 120 cells in the future, the range change will be + 20\%.

  Resolution : Note that modifying the resolution will directly influence the probability scale. Bigger classes will
  cumulate a greater number of predictions and therefore represent a greater fraction of the total
  predictions. The probability is in fact that of the class and not of isolated events.
}
\value{
  

}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{  }
\examples{

}

\keyword{ optimize }
\keyword{ distribution }