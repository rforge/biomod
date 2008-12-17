\name{ProjectionBestModel}
\alias{ProjectionBestModel}

\title{ Optimising the projections using the evaluation scores }
\description{ Selects the best model per species within Biomod according to its predictive performance
  evaluated by Roc, Kappa and/or TSS.
  
}
\usage{
ProjectionBestModel(Proj.name, Bin.trans = F, Filt.trans = F, method)
}

\arguments{
  \item{Proj.name}{ The name of your projection dataset. Note that the Projection() function has therefore to have run on this particular data.
   It is needed to load back the projections in probabilities from the 'proj' folder.}
  \item{Bin.trans}{ type True to have the projections also converted into binary values }
  \item{Filt.trans}{ type True to have the projections also converted into filtered values }
  \item{method}{ the evaluation method you want to use, can be 'Roc', 'Kappa', 'TSS' or 'all'. If 'all' is selected, then
  the best model will be chosen independently for each evaluation method }
  
}
\details{
  
}
\value{

  No values are returned but objects are produced in the "proj" folder (see \code{\link{Models}} for details).
  Example of ProjBestModelByKappa (but all follow the same logic) : a matrix with a column per species
  containing the predictions in probabilities using the Kappa scores to select the best model. The predictions
  can also be obtained in binary and filtered values (additionnal matrices are thus produced). 
  
}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{PredictionBestModel}}}
\examples{

data(Sp.Env)
Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM = T, TypeGLM = "quad", Test = "AIC", MDA = T, RF = T, GBM = F, GAM = F, CTA = F, ANN = F, SRE = F, MARS = F, NbRunEval = 0, DataSplit = 80,
ROC=T, Optimized.Threshold.ROC=T, Kappa=T, TSS=T)

PredictionBestModel(GLM=T, MDA=T, RF=T, method='all', Bin.trans = T, Filt.trans = T)

data(Future1)
Projection( Proj = Future1[,2:8], Proj.name="Future1", GLM=T, MDA=T, RF=T, BinKappa=T, FiltKappa=T)


#Determine the best projection according to each evaluation method
ProjectionBestModel(Proj.name="Future1", method='all', Bin.trans = T, Filt.trans = T)

#visualise the objects produced, example for TSS
load("pred/BestModelByTSS")
BestModelByTSS

load("proj/Proj_Future1_BestModelByTSS_Bin")
Proj_Future1_BestModelByTSS_Bin[1:100,]

load("proj/Proj_Future1_BestModelByTSS_Filt")
Proj_Future1_BestModelByTSS_Filt[1:100,]


}

\keyword{ models }
\keyword{ optimize }
