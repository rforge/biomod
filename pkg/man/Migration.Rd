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


\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{Biomod.RangeSize}}, \code{\link{Biomod.Turnover}} }
\examples{


data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,13:14], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


Projection(Proj = Sp.Env[, 4:10], Proj.name="Current", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=TRUE)

#Render Future projections under climate change scenario
data(Future1)
Projection(Proj = Future1[,4:10], Proj.name="Future1", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=TRUE)

Ensemble.Forecasting(Proj.name= "Current", weight.method='Roc', PCA.median=TRUE, 
binary=TRUE, bin.method='Roc', Test=FALSE, decay=1.6, repetition.models=TRUE)

Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=TRUE, 
binary=TRUE, bin.method='Roc', Test=FALSE, decay=1.6, repetition.models=TRUE)

load("proj.Future1/Total_consensus_Future1_Bin")
load("proj.Current/Total_consensus_Current_Bin")



Migration(CurrentPred = Total_consensus_Current_Bin[,,1], FutureProj = Total_consensus_Future1_Bin[,,1],
X=CoorXY[,1], Y=CoorXY[,2], MaxMigr=1, Pred.Save="Future1.Migration")
  

par(mfrow=c(1,3))
level.plot(Sp.Env[,13], XY=CoorXY, show.scale=FALSE, title="current distribution", cex=0.8)
level.plot(Future1.Migration[,1], XY=CoorXY, show.scale=FALSE, title="limited migration", cex=0.8)
level.plot(Total_consensus_Future1_Bin[,1,1], XY=CoorXY, show.scale=FALSE, title="unlimited migration", cex=0.8)

}

\keyword{ models }
\keyword{ optimize }
