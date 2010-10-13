\name{PredictionBestModel}
\alias{PredictionBestModel}

\title{ Optimising the predictions using the evaluation scores for BIOMOD}
\description{
  Selects the best model per species within Biomod according to its predictive performance
  evaluated by Roc, Kappa or TSS.
}
\usage{
PredictionBestModel(ANN = TRUE, CTA = TRUE, GAM = TRUE, GBM = TRUE, GLM = TRUE, MARS = TRUE, FDA = TRUE, RF = TRUE, SRE = TRUE, Bin.trans = TRUE, Filt.trans = TRUE, method = "all")
}                                                                                                                   

\arguments{
  \item{ANN}{type True for that model to be taken into account for selecting the best model (setting to False is usefull if you don't like a model that has high evaluation scores)} 
  \item{CTA}{type True for that model to be taken into account for selecting the best model} 
  \item{GAM}{type True for that model to be taken into account for selecting the best model} 
  \item{GBM}{type True for that model to be taken into account for selecting the best model} 
  \item{GLM}{type True for that model to be taken into account for selecting the best model}
  \item{MARS}{type True for that model to be taken into account for selecting the best model} 
  \item{FDA}{type True for that model to be taken into account for selecting the best model}
  \item{RF}{type True for that model to be taken into account for selecting the best model} 
  \item{SRE}{ type True for that model to be taken into account for selecting the best model }
  \item{Bin.trans}{ type True to have the predictions also converted into binary values}
  \item{Filt.trans}{ type True to have the predictions also converted into filtered values }
  \item{method}{ the evaluation method you want to use, can be 'Roc', 'Kappa', 'TSS' or 'all'. If 'all', the best model will be chosen independently for each evaluation method}
}
\details{
  The models selected but that have not been run in Models are automatically switched off. In the same line, the binary and 
  filtered procedures are only done for the evaluation methods that have been selected in Models()
}
\value{
  No values are returned but objects are produced in the "pred" folder (see \code{\link{Models}}).
  Example of PredBestModelByKappa.txt (but all follow the same logic) : a matrix with a column per run and per species
  containing the predictions in probabilities using the Kappa scores to select the best model each spieces.
  Example of PredBestModelByKappa : a 3 dimension array where the first dimention correspond to observations selected
  for each run (nb of selected observations), the second dimention refers to run and repetitions. The third coresponds
  to the spieces inice.
  The predictions can also be obtained in binary and filtered values (additionnal matrices are thus produced).
}
\references{
Thuiller, W. 2003 BIOMOD: Optimising predictions of species distributions and projecting potential future shifts under global change. Global Change Biology 9, 1353-1362.
}

\author{ Wilfried Thuiller, Bruno Lafourcade }


\seealso{ \code{\link{Models}}, \code{\link{CurrentPred}} , \code{\link{ProjectionBestModel}}}
\examples{

data(Sp.Env)
data(CoorXY)

Initial.State(Response=Sp.Env[,14:15], Explanatory=Sp.Env[,4:10],
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models. Here we will have 12 runs per species.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2,
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=2, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)

#you can do the evaluation techique individually
PredictionBestModel(method='Kappa', Bin.trans=TRUE, Filt.trans=TRUE)
PredictionBestModel(method='TSS', Bin.trans=TRUE, Filt.trans=TRUE)

#or all at once
PredictionBestModel(method='all', Bin.trans = TRUE, Filt.trans = TRUE)

#visualise the objects produced, example for Roc
#The results appear for each run that has been done
load("pred/BestModelByRoc")
BestModelByRoc

#The prediction data is stored individually for each species
#each evaluation method and each format produced (in contrary to the
#outputs of ProjectionBestModel).
load("pred/PredBestModelByRoc_Filt")
str(PredBestModelByRoc_Filt)
PredBestModelByRoc_Filt[50:100,,1] #for the first spiece

load("pred/PredBestModelByKappa_Filt")
str(PredBestModelByKappa_Filt)
PredBestModelByKappa_Filt[50:100,,2] #for the first spiece

}

\keyword{ models }
\keyword{ optimize }