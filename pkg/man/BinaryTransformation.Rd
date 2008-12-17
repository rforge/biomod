\name{BinaryTransformation}
\alias{BinaryTransformation}

\title{ Transforming data into a binary format }
\description{ Function which transform probabilities of success or occurence into a binary 
vector according to a selected threshold.
}
\usage{
BinaryTransformation(ProbData, CutOffdata)
}

\arguments{
  \item{ProbData}{ a vector or matrix containing the data to be converted }
  \item{CutOffdata}{ a numeric value (if ProbData is a vector) or a vector (if ProbData is a matrix) containing the threshold(s) 
  to be used for converting ProbData. A value needs to be given for each column of ProbData (if a matrix) }
}

\details{
 
  
}
\value{
A binary (0 or 1) vector or matrix .   
}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{  


}
\seealso{ \code{\link{help}}}
\examples{

}

\keyword{ models }