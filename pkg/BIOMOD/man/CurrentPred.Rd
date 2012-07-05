\name{CurrentPred}
\alias{CurrentPred}

\title{ Produce binary and filtered data within Biomod }
\description{
  This function enables to create binary and/or filtered datasets from the probability data and the thresholds produced
  by the Models() function of Biomod for current data.
}
\usage{
CurrentPred(GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, FDA=TRUE, MARS=TRUE, RF=TRUE,
            BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE)
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

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{PredictionBestModel}}}
\examples{
\dontrun{
data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env$Sp277, Explanatory=Sp.Env[,4:10], sp.name='Sp277',
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done a Full Calib. This will hence take several minutes.
Models(RF=TRUE, SRE=TRUE, GLM = FALSE, TypeGLM = "quad", Test = "AIC", CTA = FALSE, CV.tree = 50, ANN = FALSE, CV.ann = 2,
   NbRunEval = 1, DataSplit = 100, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


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

level.plot(Pred_Sp277_FiltTSS[,"RF",1,1], CoorXY, show.scale=FALSE)
}
}
\keyword{ models }
\keyword{ optimize }
