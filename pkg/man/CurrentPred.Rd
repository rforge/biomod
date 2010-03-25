\name{CurrentPred}
\alias{CurrentPred}

\title{ Produce binary and filtered data within Biomod }
\description{
  This function enables to create binary and/or filtered datasets from the probability data and the thresholds produced
  by the Models() function of Biomod for current data.
}
\usage{
CurrentPred(GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, FDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE)
}

\arguments{
  \item{GLM}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{GAM}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{GBM}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{CTA}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{RF}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{FDA}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{SRE}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{MARS}{ type True or False to produce binary and/or filtered prediction data for that model }
  \item{ANN}{  type True or False to produce binary and/or filtered prediction data for that model }
  \item{BinRoc}{ set to True if you want the predictions converted from probabilities to binary data using Roc }
  \item{BinKappa}{ set to True if you want the predictions converted from probabilities to binary data using Kappa }
  \item{BinTSS}{ set to True if you want the predictions converted from probabilities to binary data using TSS }
  \item{FiltRoc}{ set to True if you want the predictions converted from probabilities to filtered data using Roc }
  \item{FiltKappa}{ set to True if you want the predictions converted from probabilities to filtered data using Kappa }
  \item{FiltTSS}{ set to True if you want the predictions converted from probabilities to filtered data using TSS }
}
\details{
  The thresholds used are those stored in the Evaluation.results.Roc/Kappa/TSS objects produced by the Models() function.
  As a general reminder, do think about reloading these objects stored in the .RData file also produced by the Models() function
  from one R session to another .
}
\value{
  No values are returned but a series of objects are produced in the "pred" folder (see \code{\link{Models}} for further details).
  Like for the predictions, the files are arrays of 4 diensions and produced per species. They contain the predictions in binary or 
  filtered predictions data. See the structure in the examples below.
 
}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{PredictionBestModel}}}
\examples{

data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. Here we will have 12 runs per species. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=1, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


#The results of the Models() function are in probabilities (scaled between 0 to 1000). The CurrentPred() function
#enables to produce binary or filtered results. Here we decide to create binary data using the Roc and TSS thresholds, 
#and filtered data using only the TSS threshold.
CurrentPred(BinRoc=TRUE, BinKappa=FALSE, BinTSS=TRUE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=TRUE) 


#let's view some examples of them
load("pred/Pred_Sp277_BinRoc")
dim(Pred_Sp277_BinRoc)
dimnames(Pred_Sp277_BinRoc)[-1]
Pred_Sp277_BinRoc[120:140,,1,1]


load("pred/Pred_Sp277_FiltTSS")
Pred_Sp277_FiltTSS[120:140,,1,1]

#Note that the level.plot function will not work directly if NbRepPA > 0
#But will work on full calibrated models

load("pred/Pred_Sp277")
data(CoorXY)
nrow(Pred_Sp277)
nrow(CoorXY)
level.plot(Pred_Sp277[,"GLM",1,1], CoorXY, show.scale=FALSE)

#Make full calibrations
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=1000)


CurrentPred(BinRoc=TRUE, BinKappa=FALSE, BinTSS=TRUE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=TRUE) 


load("pred/Pred_Sp277")
load("pred/Pred_Sp277_BinRoc")
load("pred/Pred_Sp277_FiltTSS")

#plot results compared to the original data (here we plot the results for the CTA)   
#For the 100 percent calibrated model 
par(mfrow=c(1,4))
level.plot(Sp.Env[,13], CoorXY, show.scale=FALSE, title="original data")
level.plot(Pred_Sp277[,"CTA",1,1], CoorXY, show.scale=FALSE, title="probabilities with CTA full")
level.plot(Pred_Sp277_BinRoc[,"CTA",1,1], CoorXY, show.scale=FALSE, title="binary by Roc full")
level.plot(Pred_Sp277_FiltTSS[,"CTA",1,1], CoorXY, show.scale=FALSE, title="filtered by TSS full")

#For the second repetitions
x11()
par(mfrow=c(1,4))
level.plot(Sp.Env[,13], CoorXY, show.scale=FALSE, title="original data")
level.plot(Pred_Sp277[,"CTA",3,1], CoorXY, show.scale=FALSE, title="probabilities with CTA rep2")
level.plot(Pred_Sp277_BinRoc[,"CTA",3,1], CoorXY, show.scale=FALSE, title="binary by Roc rep2")
level.plot(Pred_Sp277_FiltTSS[,"CTA",3,1], CoorXY, show.scale=FALSE, title="filtered by TSS rep2")


}
\keyword{ models }
\keyword{ optimize }
