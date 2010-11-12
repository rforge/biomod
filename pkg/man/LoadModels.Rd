\name{LoadModels}
\alias{LoadModels}

\title{ Loading models produced by BIOMOD }
\description{ The function will load all the models wanted from the hardrive to R's memory for any usage.
}

\usage{
LoadModels(Sp=1, PA='all', rep='all', models='all')
}

\arguments{
  \item{Sp}{ a vector giving the species for which you want the models to be loaded. Numbers are wanted, not names. Typing 'all' selects all the species modelled.}
  \item{PA}{ a vector giving the PA runs for which you want the models to be loaded. Typing 'all' selects all the pseudo-absences runs. }
  \item{rep}{ a vector giving the repetitions for which you want the models to be loaded.  Typing 'all' selects all the repetition runs. }
  \item{models}{ a vector giving the models for which you want the models to be loaded. names are wanted. Typing 'all' selects all the models available. }
  
}

\details{
rep : '1' is the 100\% model and '2' the first repetitions
  
  
  BEWARE : storing the models outside of R is made to avoid overloading of the memory. Loading
  back a lot of models will tend towards saturation, and an exceeding number will crash R, or make it unable to run properly.
   
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{help}} }

\examples{


data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,11:18], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


ls()

LoadModels(Sp=1)
ls()
LoadModels(Sp=3, models=c("GLM", "ANN"))
ls()
LoadModels(Sp='all')
ls()

}

\keyword{ models }
