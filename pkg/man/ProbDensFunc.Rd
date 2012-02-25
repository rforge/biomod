\name{ProbDensFunc}
\alias{ProbDensFunc}

\title{ Probability Density Function }
\description{
  Using a variety of parameters in modelling will inevitably bring variability in predictions, especially when it comes to making future predictions.
  This function enables an overall viewing of the future predictions range per species and gives the likelihood of range shift estimations. It will 
  calculate the optimal way for condensing 50, 75, 90 and 95\% of the data. 
}
\usage{
ProbDensFunc(initial, projections, plothist = TRUE, cvsn = TRUE, groups = NULL, resolution = 5, save.file="no", name="ProbDF_plot")
}

\arguments{
  \item{initial}{ a vector in a binary format (ones and zeros) representing the current distribution of a species which will be used as a reference for the range change calculations }
  \item{projections}{ a matrix grouping all the predictions where each column is a single prediction. Make sure you keep projections in the same order as the initial vector (line1=site1, line2=site2, etc.). }
  \item{plothist}{ set to TRUE to plot the range change histogram }
  \item{cvsn}{ stands for "current vs new". If true, the range change calculations will be of two types: the percentage of cells currently 
  occupied by the species to be lost, and the relative percentage of cells currently unoccupied but projected to be, namely 'new' cells, compared to current surface range 
  (Please read on the BIOMOD manual for further explanation).}
  \item{groups}{ an option for ungrouping the projections enabling a separated visualisation of the prediction range per given group. A 
  matrix is expected where each column is a single prediction and each line is giving details of one parameter (See the examples section). }
  \item{resolution}{ the step used for classes of prediction in graphics. The default value is 5 }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
}

\details{
  The future range changes are calculated as a percentage of the species' present state. For example, if a species currently occupies 100 cells and 
  is estimated by a model to cover 120 cells in the future, the range change will be + 20\%.

  Resolution : Note that modifying the resolution will directly influence the probability scale. Bigger classes will
  cumulate a greater number of predictions and therefore represent a greater fraction of the total
  predictions. The probability is in fact that of the class and not of isolated events.
}
\value{
  This is a plotting function, no objects are returned or created.
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Projection}}, \code{\link{Ensemble.Forecasting}} }
\examples{

#We will do an example with 20 repetitions of the ANN and the RF to assess the variability
#in prediction making when the calibration of the model is done on partial data. Only Species7
#will be done.

\dontrun{
data(Sp.Env)

Initial.State(Response = Sp.Env[15], Explanatory = Sp.Env[,4:10], 
IndependentResponse = NULL, IndependentExplanatory = NULL)

Models(ANN=TRUE, RF=TRUE, CTA=TRUE, FDA=TRUE, GLM=TRUE, NbRunEval=20, DataSplit=60, NbRepPA=0, 
   Roc = TRUE, Optimized.Threshold.Roc = TRUE, Kappa = TRUE, TSS=TRUE)  

data(Future1)

Projection(Proj = Future1[,4:10], Proj.name="rep50_Future1", 
   ANN=TRUE, RF=TRUE, CTA=TRUE, FDA=TRUE, GLM=TRUE, BinRoc = TRUE, BinKappa = TRUE, BinTSS = TRUE, repetition.models=TRUE)              


#we will also compare with the consensus results.
Ensemble.Forecasting(Proj.name= "rep50_Future1", weight.method='Roc', PCA.median=TRUE, 
binary=TRUE, bin.method='Roc', Test=FALSE, decay=1.6, repetition.models=TRUE)


#constitute the dataframe that will be used by the ProbDensFunc 
load("proj.rep50_Future1/Proj_rep50_Future1_Sp163_BinRoc")                                                                     
load("proj.rep50_Future1/Proj_rep50_Future1_Sp163_BinKappa")
load("proj.rep50_Future1/Proj_rep50_Future1_Sp163_BinTSS")

load("proj.rep50_Future1/Total_consensus_rep50_Future1")
#But we need binary data
Bin_consensus <- BinaryTransformation(Total_consensus_rep50_Future1[,1,1], mean(consensus_rep50_Future1_results[["Sp163"]][["thresholds"]][1,]))


#Build the dataframes that will be used by the function.
#One for the projections, one for the grouping information.

models <- c("ANN", "CTA", "GLM", "FDA", "RF")
evaluation <- c("Roc", "Kappa", "TSS")

DataFrame <- matrix(NA, 2264, 2)
Groups <- matrix(NA, 2, 316)
Groups[1,] <- c(rep(rep(models, each=21), 3), "consensus")
Groups[2,] <- c(rep(evaluation, each=105), "consensus")

for(ev in evaluation){
   for(m in models){
      add.data <- eval(parse(text=paste("Proj_rep50_Future1_Sp163_Bin", ev, sep="")))
      DataFrame <- cbind(DataFrame, add.data[, m,, 1])
   }
}
DataFrame <- cbind(DataFrame, Bin_consensus)
DataFrame <- DataFrame[,-1][,-1]

ProbDensFunc(initial=Sp.Env[,15], projections=DataFrame, cvsn=TRUE, groups=Groups, resolution=5)
}
}

\keyword{ optimize }
\keyword{ distribution }