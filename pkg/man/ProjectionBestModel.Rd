\name{ProjectionBestModel}
\alias{ProjectionBestModel}

\title{ Optimising the projections using the evaluation scores }
\description{ Selects the best model per species within Biomod according to its predictive performance
  evaluated by Roc, Kappa and/or TSS.
  
}
\usage{
ProjectionBestModel(Proj.name, Bin.trans = TRUE, Filt.trans = TRUE, method='all')
}

\arguments{
  \item{Proj.name}{ The name of your projection dataset. Note that the Projection() function has therefore to have run on this particular data.
   It is needed to load back the projections in probabilities from the appropriate directory.}
  \item{Bin.trans}{ type True to have the projections also converted into binary values }
  \item{Filt.trans}{ type True to have the projections also converted into filtered values }
  \item{method}{ the evaluation method you want to use, can be 'Roc', 'Kappa', 'TSS' or 'all'. If 'all' is selected, then
  the best model will be chosen independently for each evaluation method }
  
}
\details{
  The Proj.name is dependent of the one already given in \code{\link{Projection}}. It will be used to restore some
  data in the correspondant directory.
}

\value{
  No values are returned but objects are produced in the appropriate folder (see \code{\link{Projection}} for details).
  Example of ProjBestModelByKappa (but all follow the same logic) : a matrix with a column per species
  containing the projections in probabilities using the Kappa scores to select the best model. The projections
  can also be obtained in binary and filtered values (additionnal matrices are thus produced). 
}

\references{ 
Thuiller, W. 2003 BIOMOD: Optimising predictions of species distributions and projecting potential future shifts under global change. Global Change Biology 9, 1353-1362.
}
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{PredictionBestModel}}}
\examples{

data(Sp.Env)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. Here we will have 12 runs per species. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=2, strategy="sre", coor=NULL, distance=0, nb.absences=1000)


#Render Future projections under climate change scenario
data(Future1)
Projection( Proj = Future1[,4:10], Proj.name="Future1", GLM=TRUE, MDA=TRUE, RF=TRUE, BinKappa=TRUE, FiltKappa=TRUE)


#This function needs to be run before hand in order to produce the reference files for ProjectionBestModel
PredictionBestModel(method='all', Bin.trans = TRUE, Filt.trans = TRUE)

#Determine the best projection according to each evaluation method
ProjectionBestModel(Proj.name="Future1", method='all', Bin.trans = TRUE, Filt.trans = TRUE)

#visualise the objects produced, example for TSS
load("pred/BestModelByTSS")
BestModelByTSS

#In contrast with PredictionBestModel(), the results are stored in one array per method with
#the third dimension being the number of species modelled.

load("proj.Future1/Proj_Future1_BestModelByTSS_Bin")
dim(Proj_Future1_BestModelByTSS_Bin)
dimnames(Proj_Future1_BestModelByTSS_Bin)[-1]
Proj_Future1_BestModelByTSS_Bin[1:10,,]

load("proj.Future1/Proj_Future1_BestModelByTSS_Filt")
Proj_Future1_BestModelByTSS_Filt[1:10,,]

}

\keyword{ models }
\keyword{ optimize }
