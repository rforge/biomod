\name{CurrentPred}
\alias{CurrentPred}

\title{ Produce binary and filtered data within Biomod }
\description{
  This function enables to create binary and/or filtered datasets from the probability data and the thresholds produced
  by the Models() function of Biomod for current data.
}
\usage{
CurrentPred(GLM, GBM, GAM, CTA, ANN, SRE, MDA, MARS, RF, BinRoc, BinKappa, BinTSS, FiltRoc, FiltKappa, FiltTSS)
}

\arguments{
  \item{GLM}{,}
  \item{GAM}{,}
  \item{GBM}{,}
  \item{CTA}{,}
  \item{RF}{,}
  \item{MDA}{,}
  \item{SRE}{,}
  \item{MARS}{and}
  \item{ANN}{ : select the models by typing True or False for which you want the predictions to be produced in binary and/or filtered format}
  \item{BinRoc}{ set to True if you want the predictions converted from probabilities to binary data using Roc }
  \item{BinKappa}{ set to True if you want the predictions converted from probabilities to binary data using Kappa }
  \item{BinTSS}{ set to True if you want the predictions converted from probabilities to binary data using TSS }
  \item{FiltRoc}{ set to True if you want the predictions converted from probabilities to filtered data using Roc }
  \item{FiltKappa}{ set to True if you want the predictions converted from probabilities to filtered data using Kappa }
  \item{FiltTSS}{ set to True if you want the predictions converted from probabilities to filtered data using TSS }
}
\details{
  The thresholds used are those stored in the Evaluation.results.Roc/Kappa/TSS objects produced by the Models() function.
  As a general reminder, do think about reloading these objects stored in the .RData file also produced by the Models() function from one R session to another .
}
\value{
  No values are returned but a series of objects are produced in the "pred" folder (see \code{\link{Models}} for further details).
  The files are produced per species and contain the predictions in binary or filtered format with a column per model selected. See the
  structure in the examples below.
 
}
\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{PredictionBestModel}}}
\examples{

data(Sp.Env)

Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM = T, TypeGLM = "quad", Test = "AIC", MDA = T, RF = T, GBM = F, GAM = F, CTA = F, ANN = F, SRE = F, MARS = F, NbRunEval = 0, DataSplit = 80,
ROC=T, Optimized.Threshold.ROC=T, Kappa=T, TSS=T)


#The results of the Models() function are in probabilities (scaled between 0 to 1000). The CurrentPred() function
#enables to produce binary or filtered results. Here we decide to create binary data using the Roc and TSS thresholds, 
#and filtered data using only the TSS threshold.
CurrentPred(BinRoc=T, BinKappa=F, BinTSS=T, FiltRoc=F, FiltKappa=F, FiltTSS=T) 


#let's view some examples of them
load("pred/Pred_Species2_BinRoc")
Pred_Species2_BinRoc[120:140,]      #the head() function reads the first six lines of a dataframe


load("pred/Pred_Species2_FiltTSS")
Pred_Species2_FiltTSS[120:140,]        #the head() function reads the first six lines of a dataframe


#and even plot them compared to the original data (here we plot the results for the GLM)
load("pred/Pred_Species2")
data(CoorXY)


par(mfrow=c(1,4))
level.plot(Sp.Env[,10], CoorXY, show.scale=F, title="original data")
level.plot(Pred_Species2[,"GLM"], CoorXY, show.scale=F, title="probabilities with GLM")
level.plot(Pred_Species2_BinRoc[,"GLM"], CoorXY, show.scale=F, title="binary by Roc")
level.plot(Pred_Species2_FiltTSS[,"GLM"], CoorXY, show.scale=F, title="filtered by TSS")


}
\keyword{ models }
\keyword{ optimize }
