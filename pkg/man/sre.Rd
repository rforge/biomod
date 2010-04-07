\name{sre}
\alias{sre}

\title{ Surface Range Envelope }
\description{
  Run a rectilinear surface range envelop (equivalent to BIOCLIM) using the extreme percentiles as recommended by Nix or Busby.
  The SRE performs a simple analysis of within which range of each variable the data is recorded and renders predictions.
}
\usage{
sre(Response = NULL, Explanatory = NULL, NewData = NULL, quants=0.05)
}

\arguments{
  \item{Response}{ a vector giving your species' presences and absences data }
  \item{Explanatory}{ a matrix containing the environmental variables for the sites given in Response. It must
   have as many rows as there are elements in Response. }
  \item{NewData}{ The data for which you want to render predictions with the sre. It must be a matrix of the same type as
   the one given in Explanatory and with precisely the same variable names. }
  \item{quants}{ the value defines the most extreme values for each variable not to be taken into account
   for determining the tolerance boundaries for the considered species.  }
}
\details{

  The more variables you put in, the more restrictive your model will be (if non-colinear variables).

  This method is very much influenced by the data input, and more specifically by the extremes.
  Where a linear model can discriminate the extreme values from the main tendency, the SRE considers it
  equal as any other data point which leads to notable changes in predictions.
  Note that, as a consequence of its functionning, the predictions are directly given in binary, a site being
  either potentially suitable for all the variables, either out of bounds for at least one variable and therefore
  considered unsuitable. 
  
  The quants argument determines the threshold at which the data will be taken into account for calibration : the default
  of 0.05 induces that the 5\% most extreme values will be avoided for each variable on each side of its distribution along
  the gradient. So it in fact takes 5\% away at each end of the variables distribution, giving a total of 10\% of data not considered.
}
\value{
  A vector of the same length as there are rows in NewData giving the prediction in binary (1=presence, 0=absence)
}

\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Models}} }
\examples{

data(Sp.Env)

g <- sre(Sp.Env[,10], Sp.Env[,2:8], NewData=Sp.Env[,2:8])
gg <- sre(Sp.Env[,10], Sp.Env[,2:8], NewData=Sp.Env[,2:8], quants=0.025)
ggg <- sre(Sp.Env[,10], Sp.Env[,2:8], NewData=Sp.Env[,2:8], quants=0.05)


par(mfrow=c(1,4))
level.plot(Sp.Env[,10], XY=CoorXY[,1:2], show.scale=FALSE, title="original data")
level.plot(g, XY=CoorXY[,1:2], show.scale=FALSE, title="full data calibration")
level.plot(gg, XY=CoorXY[,1:2], show.scale=FALSE, title="Perc025=TRUE")
level.plot(ggg, XY=CoorXY[,1:2], show.scale=FALSE, title="Perc05=TRUE")


}

\keyword{ models }
\keyword{ multivariate }