\name{Migration}
\alias{Migration}

\title{ Using migration dynamics in species modelling }
\description{
  This function allows the inclusion of a very simple migration process when projecting species distributions
  into the future. Two datasets are needed : the current species distributions and the future (assuming by
  default unlimited migration).
}
\usage{
Migration(CurrentPred = NULL, FutureProj = NULL, X = NULL, Y = NULL, MaxMigr = NULL, Pred.Save = NULL)
}

\arguments{
  \item{CurrentPred}{ a matrix of n (number of species) columns giving the current state of the species in binary  }
  \item{FutureProj}{ a matrix of n (number of species) columns giving the future predictions of the species in binary }
  \item{X}{ , }
  \item{Y}{ X and Y : the latitude and longitude giving the coordinates of the data points given in CurrentPred and FutureProj }
  \item{MaxMigr}{ a single value or vector, of length the number of species modelled, giving the maximum distance the species 
  could migrate according to the time slice modelled (see the details section). }
  \item{Pred.Save}{ the name where the projections using limited migration will be stored }
}
\details{
  To be able to use this function, both current and future datasets must be ordered in the
 same way and have the same coordinates and resolution.

 Note that the rate of migration should be given in degrees. For instance a species could only have a
 maximum migration of 1 minute (1.6km) per decade. If we project its distribution in 50 years, the rate is
 1 x 0.16667 x 5 (where 0.1667 is the conversion from minute to degree).
 For projection in 2080 : Rate = 1 x 0.16667 x 8.
 For a maximum rate of 3 minutes per decade (4.8km) in 2080: Rate = 3 x 0.16667 x 8
}
\value{ 




}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{Biomod.RangeSize}}, \code{\link{Biomod.Turnover}} }
\examples{



data(Sp.Env)

Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM = T, TypeGLM = "quad", Test = "AIC", MDA = T, RF = T, GBM = F, GAM = F, CTA = F, ANN = F, SRE = F, MARS = F, NbRunEval = 0, DataSplit = 80,
Roc=T, Optimized.Threshold.Roc=T, Kappa=T, TSS=T)


data(Future1)
Projection( Proj = Future1[,2:8], Proj.name="Future1", GLM=T, MDA=T, RF=T, BinTSS=T)


#We are going to run the Turnover function on the models selected by the TSS for each species
PredictionBestModel(method="Roc")
ProjectionBestModel(Proj.name="Future1", method="Roc", Bin.trans=T)



load("pred/PredBestModelByRoc_Bin")
load("proj.Future1/Proj_Future1_BestModelByRoc_Bin")

data(CoorXY)
Migration(CurrentPred = PredBestModelByRoc_Bin, FutureProj = Proj_Future1_BestModelByRoc_Bin,
  X=CoorXY[,1], Y=CoorXY[,2], MaxMigr=5*0.1666*8, Pred.Save="Future1.Migration.BestRoc")


cbind(Proj_Future1_BestModelByRoc_Bin[,1], Future1.Migration.BestRoc[,1], Proj_Future1_BestModelByRoc_Bin[,2], Future1.Migration.BestRoc[,2])[700:730,]


par(mfrow=c(1,2))
level.plot(Proj_Future1_BestModelByRoc_Bin[,1], XY=CoorXY, show.scale=F, title="unlimited migration")
level.plot(Future1.Migration.BestRoc[,1], XY=CoorXY, show.scale=F, title="limited migration")



}

\keyword{ models }
\keyword{ optimize }
