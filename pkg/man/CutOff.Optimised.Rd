\name{CutOff.Optimised}
\alias{CutOff.Optimised}

\title{ Inner Function }
\description{ Inner function for the BIOMOD package, not for independent usage purposes
}
\usage{
CutOff.Optimised(Obs, Fit)
}

\arguments{
  \item{Obs}{  }
  \item{Fit}{  }
}
\details{
Function which takes a binary vector (i.e. observed species presence and absence)
and a vector of probabilities (i.e. probability of presence from a model) and search 
for the threshold which will jointly maximise the percentage of presence correctly predicted (sensitivity)
and the percentage of absence correctly predicted (specificity).  
This threshold could be further used to transform the vector a probabilities into
a binary vector (presence and absence)  
}
\value{
A vector with three values: threshold, sensitivity and specificity  
}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{  


}
\seealso{ \code{\link{help}}}
\examples{

}

\keyword{ models }