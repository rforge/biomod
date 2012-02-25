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
\dontrun{
data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env$Sp281, Explanatory=Sp.Env[,4:10], sp.name='Sp281',
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are full data run . This will hence take several minutes.
Models(RF=TRUE, SRE=TRUE, GLM = FALSE, TypeGLM = "quad", Test = "AIC", CTA = FALSE, CV.tree = 50, ANN = FALSE, CV.ann = 2,
   NbRunEval = 1, DataSplit = 100, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)
   
#Load the predictions
load("pred/Pred_Sp281")	
	
## Cutoff estimation estimation

CutOff.Optimised(Sp.Env$Sp281, Pred_Sp281[,"RF",1,1])

# TSS estimated inside the Models function
Evaluation.results.Roc	
}

}

\keyword{ models }