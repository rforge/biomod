\name{sre}
\alias{sre}

\title{ Surface Range Envelope }
\description{
  Run a rectilinear surface range envelop (equivalent to BIOCLIM) using the percentile 0.025 or 0.05 as recommended by Nix or Busby.
  The SRE performs a simple analysis of within which range of each variable the data is recorded to render predictions.
}
\usage{
sre(Response = NULL, Explanatory = NULL, NewData = NULL, Perc025 = F, Perc05 = F)
}

\arguments{
  \item{Response}{ a vector giving your species' presences and absences data }
  \item{Explanatory}{ a matrix containing the environmental variables for the sites given in Response. It must
  have as many rows as there are elements in Response. }
  \item{NewData}{  }
  \item{Perc025}{ if True, the 2.5\% most extreme values for each variable will not be taken into account
  for determining the tolerance boundaries for the considered species.  }
  \item{Perc05}{ same as for  Perc025 but the 5\% most extreme values will be avoided }
}
\details{

  The more variables you put in, the more restrictive your model will be (if non-colinear variables).

  This method is very much influenced by the data input, and more specifically by the extremes.
  Where a linear model can discriminate the extreme values from the main tendency, the SRE considers it
  equal as any other data point which leads to notable changes in predictions.
  Note that, as a consequence of its functionning, the predictions are directly given in binary, a site being
  either potentially suitable for all the variables, either out of bounds for at least one variable and therefore
  considered unsuitable. 
}
\value{
  Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  
}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{  }
\examples{

}

\keyword{ models }
\keyword{ multivariate }