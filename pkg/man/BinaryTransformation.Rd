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
  \item{ProbData}{ a vector containing the data to be converted }
  \item{CutOffdata}{ a numeric value containing the threshold to be used for converting ProbData.}
}


\value{
A binary (0 or 1) vector.   
}


\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{help}}}
\examples{

xx <- rnorm(50,10)
yy <- BinaryTransformation(xx, 10)

cbind(xx,yy)

}

\keyword{ models }