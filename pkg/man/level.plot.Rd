\name{level.plot}
\alias{level.plot}

\title{ Plotting 2-dimensional data }
\description{
  Enables to plot data with 2-dimensional geographical coordinates
}
\usage{
level.plot(data.in, XY, color.gradient = "red", cex = 1, level.range = c(min(data.in), max(data.in)), show.scale = T, title = "level plot")
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
}

\details{}

\value{}

\author{ Bruno Lafourcade }

\seealso{ \code{\link{map.plot}} }
\examples{
  
data(Sp.Env)
data(CoorXY) 

level.plot(Sp.Env[,3], CoorXY, title='var 3')
x11()
level.plot(Sp.Env[,3], CoorXY, color.gradient='grey', title='var 3') 

#add in some extreme values
Sp.Env[15,3] <- 20000   
Sp.Env[16,3] <- -5000 

x11()
level.plot(Sp.Env[,3], CoorXY, title='var 3')  

#the extreme values imply a non satisfying view of our data. The contrasts are not clear enough 
x11()
hist(Sp.Env[,3], breaks=50, main='variable 3')
#Let's narrow the range :
x11()
level.plot(Sp.Env[,3], CoorXY, level.range=c(0,5000), title='var 3')
 
}

\keyword{ aplot }
