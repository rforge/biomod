\name{PredictionBestModel}
\alias{PredictionBestModel}

\title{ Optimising the predictions using the evaluation scores for BIOMOD}
\description{
  Selects the best model per species within Biomod according to its predictive performance
  evaluated by Roc, Kappa or TSS.
}
\usage{
PredictionBestModel(ANN = T, CTA = T, GAM = T, GBM = T, GLM = T, MARS = T, MDA = T, RF = T, SRE = T, Bin.trans = T, Filt.trans = T, method = "all")
}

\arguments{
  \item{ANN}{,} 
  \item{CTA}{,} 
  \item{GAM}{,} 
  \item{GBM}{,} 
  \item{GLM}{,}
  \item{MARS}{,} 
  \item{MDA}{,}
  \item{RF}{and} 
  \item{SRE}{ : select the models by typing True or False to be taken into account for selecting the best model (setting to False is usefull 
  if you don't like a model that has high evaluation scores) }
  \item{Bin.trans}{ type True to have the predictions also converted into binary values}
  \item{Filt.trans}{ type True to have the predictions also converted into filtered values }
  \item{method}{ the evaluation method you want to use, can be 'Roc', 'Kappa', 'TSS' or 'all'. If 'all' is selected, then
  the best model will be chosen independently for each evaluation method}
}
\details{
  
}
\value{
  No values are returned but objects are produced in the "pred" folder (see \code{\link{Models}}).
  Example of PredBestModelByKappa (but all follow the same logic) : a matrix with a column per species
  containing the predictions in probabilities using the Kappa scores to select the best model. The predictions
  can also be obtained in binary and filtered values (additionnal matrices are thus produced).  
}
\references{
Thuiller, W. 2003 BIOMOD: Optimising predictions of species distributions and projecting potential future shifts under global change. Global Change Biology 9, 1353-1362.
}

\author{ Wilfried Thuiller, Bruno Lafourcade }


\seealso{ \code{\link{Models}}, \code{\link{CurrentPred}} , \code{\link{ProjectionBestModel}}}
\examples{

data(Sp.Env)
Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM = T, TypeGLM = "quad", Test = "AIC", MDA = T, RF = T, GBM = F, GAM = F, CTA = F, ANN = F, SRE = F, MARS = F, NbRunEval = 0, DataSplit = 80,
Roc=T, Optimized.Threshold.Roc=T, Kappa=T, TSS=T)

#you can do the evaluation techique individually
PredictionBestModel(GLM=T,GBM=T, GAM=T, CTA=T, ANN=T, MDA=T, MARS=T, RF=T, SRE=T, method='Kappa', Bin.trans = T, Filt.trans = T)
PredictionBestModel(GLM=T,GBM=T, GAM=T, CTA=T, ANN=T, MDA=T, MARS=T, RF=T, SRE=T, method='TSS', Bin.trans = T, Filt.trans = T)

#or all at once
PredictionBestModel(GLM=T,GBM=T, GAM=T, CTA=T, ANN=T, MDA=T, MARS=T, RF=T, SRE=T, method='all', Bin.trans = T, Filt.trans = T)

#visualise the objects produced, example for Roc
load("pred/BestModelByRoc")
BestModelByRoc

load("pred/PredBestModelByRoc_Filt")
PredBestModelByRoc_Filt[50:100,]


}

\keyword{ models }
\keyword{ optimize }