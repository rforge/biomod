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
Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. Here we will have 12 runs per species. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", NbRunEval = 0, DataSplit = 100, TSS=TRUE)

#Load the predictions
load("pred/Pred_Sp277")	
	
## TSS estimation

MyTSS <- KappaRepet(Sp.Env$Sp277, Pred_Sp277[,"GLM",1,1], TSS = TRUE)
MyTSS

ls()
# TSS estimated inside the Models function
Evaluation.results.TSS	
	
}

\keyword{ models }