\name{level.plot}
\alias{level.plot}

\title{ Plotting 2-dimensional data }
\description{
  Enables to plot data with 2-dimensional geographical coordinates
}
\usage{
level.plot(data.in, XY, color.gradient = "red", cex = 1, level.range = c(min(data.in), max(data.in)), show.scale = TRUE, title = "level plot", SRC=FALSE, save.file="no", ImageSize="small", AddPresAbs=NULL, PresAbsSymbol=c(cex*0.8,16,4))
}

\arguments{
  \item{data.in}{ the data you want to visualise }
  \item{XY}{ a 2 columns matrix of the same length as data.in giving the coordinates of the data points}
  \item{color.gradient}{ available : red, grey and blue }
  \item{cex}{ to change the point size : >1 will increase size, <1 will decrease it }
  \item{level.range}{ the range of values for which you want the color gradient to be used, usefull to increase the resolution
  of the graph especially if you have extreme values (see examples section) }
  \item{show.scale}{ a feature for just keeping the graph without the scale (if set to False) }
  \item{title}{ the title wanted for the graph }
  \item{SRC}{ if TRUE, the function recognizes the values as being 'Species Range Change' values and associates the appropriate colors to them (see \code{\link{Biomod.RangeSize}})}
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{ImageSize}{ The image size for JPEG and TIFF files saved to disk. Available : 'small', 'standard' and 'large' }
  \item{AddPresAbs}{ Optional: adds the presences and absences used for calibration to the plot. Data must be entered as a matrix/dataframe with 3 columns (in this order): X-coordinate, Y-coordinate, Presence(1) or Absence(0). X and Y coordinates must be in the same system as the plot. }
  \item{PresAbsSymbol}{ Optional: a 3 element vector giving the symbols to be used by the AddPresAbs argument for plotting. The elements of the vector must be in this order: size of presence/absence symbols given as a multiplication factor of the 'cex' value entered in the function (e.g. a value of 0.5 means that the symbols will be drawn at a size = 0.5*cex value entered in the function), symbol (in PCH code) to be used for presences, symbol (in PCH code) to be used for absences. An example of input vector for this parameter is 'c(0.4,16,4)' }
}

\author{ Bruno Lafourcade }

\seealso{ \code{\link{multiple.plot}} }
\examples{
  
data(Sp.Env)
data(CoorXY) 

level.plot(Sp.Env[,5], CoorXY, title='var 3')
#x11()
level.plot(Sp.Env[,5], CoorXY, color.gradient='grey', title='var 3') 

#add in some extreme values
Sp.Env[15,5] <- 20000   
Sp.Env[16,5] <- -5000 

#x11()
level.plot(Sp.Env[,5], CoorXY, title='var 3')  

#the extreme values imply a non satisfying view of our data. The contrasts are not clear enough 
#x11()
hist(Sp.Env[,5], breaks=50, main='variable 3')
#Let's narrow the range :
#x11()
level.plot(Sp.Env[,5], CoorXY, level.range=c(0,5000), title='var 3')
 
}

\keyword{ aplot }
