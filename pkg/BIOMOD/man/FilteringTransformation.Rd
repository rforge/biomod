\name{FilteringTransformation}
\alias{FilteringTransformation}


\title{ Transforming data into a filtered format }
\description{ Function which transform probabilities of success or occurence into a filtered 
vector according to a selected threshold.
}
\usage{
FilteringTransformation(ProbData, CutOffdata)
}

\arguments{
  \item{ProbData}{ a vector containing the data to be converted }
  \item{CutOffdata}{ a numeric value containing the threshold to be used for converting ProbData.}
}

\details{
The values higher than the threshold are not converted to a "1" but kept while the values
lower than the threshold are transformed into a zero.
  
}
\value{
A vector of zeros and original values of ProbData.
 
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{help}}}
\examples{

xx <- rnorm(50,10)
yy <- FilteringTransformation(xx, 10)

cbind(xx,yy)

}

\keyword{ models }

