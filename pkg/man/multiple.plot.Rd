\name{multiple.plot}
\alias{multiple.plot}

\title{ Plot and compare prediction maps within BIOMOD }
\description{
  This function allows direct comparisons of the predictions across models. The
  maps produced are the same than with the \code{\link{level.plot}} function 
}
\usage{
multiple.plot(Data, coor, color.gradient='red', plots.per.window=9, cex=1, save.file="no", name="multiple plot")
}

\arguments{
  \item{Data}{ the data you want to visualise }
  \item{coor}{ a 2 columns matrix of the same length as data.in giving the coordinates of the data points }  
  \item{color.gradient}{ available : 'red', 'grey' and 'blue' }
  \item{plots.per.window}{ the number of plots you wish to visualise per window }
  \item{cex}{ to change the point size : >1 will increase size, <1 will decrease it }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ The image size for JPEG and TIFF files saved to disk. Available : 'small', 'standard' and 'large' }
}

\details{
}
\value{
}

\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{level.plot}} }
\examples{


multiple.plot(Sp.Env[,4:18], CoorXY, plots.per.window=16, cex=0.8)

#just the set of environmental variables 
multiple.plot(Sp.Env[,4:10], CoorXY, cex=0.8)

}

\keyword{ aplot }
