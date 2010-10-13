\name{Biomod.Turnover}
\alias{Biomod.Turnover}
\alias{Biomod.Turnover.raster}

\title{ Estimation of the Turnover for a given location and resolution}
\description{
  This function allows to estimate the species lost and gained by pixel for the time slice and scenario considered, 
  and the turnover induced by these (predicted) modifications.
}
\usage{
Biomod.Turnover(CurrentPred = NULL, FutureProj = NULL, Turnover.Save = NULL)
Biomod.Turnover.raster(CurrentPred = NULL, FutureProj = NULL, Turnover.Save = NULL)
}

\arguments{
  \item{CurrentPred}{ a matrix of n (number of species) columns OR a rasterstack giving the current state of the species in binary  }
  \item{FutureProj}{ a matrix of n (number of species) columns OR a rasterstack giving the future predictions of the species in binary }
  \item{Turnover.Save}{ the name given to the new object storing the results }
}
\details{
  
}
\value{
  A 10 columns database or a 10 layers rasterstack is created (and saved under the name you gave it).
The first four are relative numbers: 
- Disa represents the number of species predicted to disappear from the given pixel. 
- Stable0 is the number of species which are currently not in the given pixel and not predicted to migrate into it. 
- Stable1 represents the number of species currently occurring in the given pixel, and predicted to remains into the future. 
- Gain represents the number of species which are currently absent but predicted to migrate in the given pixel.

PercLoss, PercGain and Turnover are the relative percentages estimated as the following:
- PercLoss = 100 x L/(SR)
- PercGain = 100 x G/(SR)
- Turnover = 100 x (L+G)/(SR+G)

Where SR is the current species richness, L the number of species lost and G those gained.

- CurrentSR represents the current modelled species richness in the given pixel.
- FutureSR0Disp represents the future modelled species richness assuming no migration of species
- FutureSR1Disp represents the future modelled species richness assuming migration (depending on
  the datasets given in input, if \code{\link{Migration}} has been used or not).

}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Biomod.RangeSize}}, \code{\link{Migration}} }
\examples{

 
data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,11:18], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. This will hence take several minutes.
Models(GLM = FALSE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


Projection(Proj = Sp.Env[, 4:10], Proj.name="Current", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=T)

#Render Future projections under climate change scenario
data(Future1)
Projection(Proj = Future1[,4:10], Proj.name="Future1", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=T)

Ensemble.Forecasting(Proj.name= "Current", weight.method='Roc', PCA.median=T, 
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=T, 
binary=T, bin.method='Roc', Test=F, decay=1.6, repetition.models=T)

load("proj.Future1/Total_consensus_Future1_Bin")
load("proj.Current/Total_consensus_Current_Bin")


Biomod.Turnover(CurrentPred = Total_consensus_Current_Bin[,,1], FutureProj = Total_consensus_Future1_Bin[,,1], Turnover.Save="Future1.Turnover.consensus")


#Visualise the results
Future1.Turnover.consensus[1150:1200,]


}

\keyword{ IO }