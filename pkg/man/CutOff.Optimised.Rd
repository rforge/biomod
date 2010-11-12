\name{CutOff.Optimised}
\alias{CutOff.Optimised}

\title{ Treshold estimation minimising difference between sensitivity and specificity }
\description{ Function that calculates the threshold minimising the absolute difference between sensitivity and specificity
}
\usage{
CutOff.Optimised(Obs, Fit)
}

\arguments{
  \item{Obs}{ Observed data usually presence and absence of a given species }
  \item{Fit}{ Probability of occurrence from a species distribution model }
}
\details{
Function which takes a binary vector (i.e. observed species presence and absence)
and a vector of probabilities (i.e. probability of presence from a model) and search 
for the threshold which will jointly maximise the percentage of presence correctly predicted (sensitivity)
and the percentage of absence correctly predicted (specificity). In other words, this is the treshold minimising the difference between sensitivity and specificity.   
This threshold could be further used to transform the vector of probabilities into
a binary vector (presence and absence)  
}
\value{
A vector with three values: threshold (Cutoff), sensitivity (se) and specificity (sp).   
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{help}}}
\examples{

data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. Here we will have 12 runs per species. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", NbRunEval = 0, DataSplit = 100, Roc=TRUE, Optimized.Threshold.Roc = FALSE)

#Load the predictions
load("pred/Pred_Sp277")	
	
## Cutoff estimation estimation

CutOff.Optimised(Sp.Env$Sp277, Pred_Sp277[,"GLM",1,1])

# TSS estimated inside the Models function
Evaluation.results.Roc	


}

\keyword{ models }