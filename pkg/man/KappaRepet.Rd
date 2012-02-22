\name{KappaRepet}
\alias{KappaRepet}


\title{ Maximising Kappa and True Skill Statistics for incremental threshold }
\description{ Find the threshold maximising either Kappa or True Skill Statistic
}

\usage{
KappaRepet(Obs, Fit, TSS = FALSE)
}

\arguments{
  \item{Obs}{ vector of observed data, coded in 0 (absent) and 1 (present) }
  \item{Fit}{ vector of predicted data from distribution models, usually ranging from 0 to 1 }
  \item{TSS}{ logical. Should TSS be calculated instead of Kappa? }
}

\details{
  Build a mis-classification matrix for incremental thresholds, estimate the Kappa or TSS and find the threshold that maximise on the statistics.
}
\value{
  A vector with Kappa or TSS, the selected threshold (CutOff), the true positive (TP), the sensitivity (se) that is the percentage of true positive, the true negative (TN) and the specificity (sp) that is the percentage of true false positive.
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{help}}}
\examples{
	
data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env$Sp281, Explanatory=Sp.Env[,4:10], sp.name='Sp281',
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done a lone full run. This will hence take several minutes.
Models(RF=TRUE, SRE=TRUE, GLM = FALSE, TypeGLM = "quad", Test = "AIC", CTA = FALSE, CV.tree = 50, ANN = FALSE, CV.ann = 2,
   NbRunEval = 1, DataSplit = 100, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)

#Load the predictions
load("pred/Pred_Sp281")	
	
## TSS estimation

MyTSS <- KappaRepet(Sp.Env$Sp281, Pred_Sp281[,"RF",1,1], TSS = TRUE)
MyTSS

ls()
# TSS estimated inside the Models function
Evaluation.results.TSS	
	
}

\keyword{ models }