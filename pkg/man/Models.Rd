\name{Models}
\alias{Models}

\title{ Running the species-climate envelop models available in BIOMOD }
\description{
  This function allows to train various species-climate modelling algorithms. It enables also
  to run three evaluation techniques (TSS, Kappa, Roc curve) to estimate the predictive
  performance of the models produced. 
}
\usage{
Models(GLM = FALSE, TypeGLM = "simple", Test = "AIC", GBM = FALSE, No.trees = 5000, GAM = FALSE, Spline = 3, CTA = FALSE, CV.tree = 50, ANN = FALSE, CV.ann = 5, SRE = FALSE, quant=0.025, FDA = FALSE, MARS = FALSE, RF = FALSE, NbRunEval = 1, DataSplit = 100, NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=NULL, Yweights = NULL, VarImport = 0, Roc = FALSE, Optimized.Threshold.Roc = FALSE, Kappa = FALSE, TSS = FALSE, KeepPredIndependent = FALSE)
}

\arguments{
  \item{GLM}{ type True to run the Generalised Linear Model, or False to switch it off }
  \item{TypeGLM}{ the complexity of the terms of the GLM : 'simple', 'poly', 'quad' standing for linear, polynomial or quadratic respectively (see manual for further details). }
  \item{Test}{ the stepwise procedure for the GLM, either use the AIC or BIC criteria. }
  \item{GBM}{ type True to run the Generalised Boosting Model, or False to switch it off}
  \item{No.trees}{ the maximum number of trees for the GBM }
  \item{GAM}{ type True to run the Generalised Additive Model, or False to switch it off }
  \item{Spline}{ the degree of smoothing of the spline function for the GAM (3 by default)}
  \item{CTA}{ type True to run the Classification and regression Tree Analysis, or False to switch it off }
  \item{CV.tree}{ number of cross-validation in CTA to estimate the length of the optimal tree }
  \item{ANN}{ type True to run the Artificial Neural Network, or False to switch it off }
  \item{CV.ann}{ number of cross validation for the ANN }
  \item{SRE}{ type True to run the Surface Range Envelope, or False to switch it off }
  \item{quant}{ the value defines the extreme percentiles of the response data not to be used by the SRE for calibration (see \code{\link{sre}} for more details) }
  \item{FDA}{ type True to run the Flexible Discriminant Analysis, or False to switch it off  }
  \item{MARS}{ type True to run the Multiple Adaptive Regression Splines, or False to switch it off }
  \item{RF}{ type True to run the Random Forest, or False to switch it off }
  \item{NbRunEval}{ the number of evaluation runs to proceed before the final 100 percent data model is run (see the manual or the details section below for further explanations) }
  \item{DataSplit}{ to determine the ratio used for splitting the original database in calibration and evaluation subsets. The value given (between 0 and 100) is the percentage allocated to calibration }
  \item{NbRepPA}{ The number of repetitions for the selection of pseudo absences}
  \item{strategy}{ the strategy to use for selecting pseudo absences. Can be either "circles", "squares", "per", "random" or "sre". See \code{\link{pseudo.abs}} for details.}
  \item{coor}{ a two columned matrix giving the coordinates of your data points. It is needed for the "per", "circles" and "squares" strategies}
  \item{distance}{ a value giving the distance to use for the "per", "circles" and "squares" strategy of the pseudo absences selection. See \code{\link{pseudo.abs}} for details.}
  \item{nb.absences}{ the number of pseudo absences wanted to run the models with. They are randomly selected from the pool of pseudo absences available selected by the given strategy.}
  \item{Yweights}{ a N-columns matrix, N being the number of species, giving the weights that the user can set for the response variables. 
   This is similar to an index of detectability for each site, which allows users to give stronger weights to more reliable presences or absences }
  \item{VarImport}{ if True, this parameter enables a direct comparison of the explanatory variable importance across models (see the manual for further explanations) }
  \item{Roc}{ If True, the Area Under the ROC curve score will be evaluated for all the models selected and for each species }
  \item{Optimized.Threshold.Roc}{ Type True to determine a threshold with the Roc method optimising the percentage of presences
   and absences correctly predicted (this is not automatic as the Roc curve is a threshold independent method). }
  \item{Kappa}{ If True, the Kappa score will be evaluated and a threshold produced for all the models selected and for each species }
  \item{TSS}{ If True, the True Skill Statistic score will be evaluated and a threshold produced for all the models selected and for each species }
  \item{KeepPredIndependent}{ If true, the predictions on the independent data will be saved. Otherwise, only the predictive accuracy on the evaluation set are conserved (NbRunEval) }
}

\details{
  The dataset used is DataBIOMOD produced by the \code{\link{Initial.State}} function. A prior run of this function is therefore necessary. The
  two arguments explained below are crucial as they determine the number of runs that will be done for each species. Be carefull not to over
  multiply thevtotal number of runs in order to keep a reasonable running time. 
  
  NbRepPA : BIOMOD now enables to calibrate the models with only a limited number of absences from the full dataset given in \code{\link{Initial.State}}.
  If NbRepPA is equal to 0 the calibration procedure will be the standard one, i.e. taking all the data available. Else (NbRepPA > 0), an inner run of the 
  \code{\link{pseudo.abs}} function will be done for each species, with the strategy chosen, to determine the pool of pseudo-absences available for calibration.
  The number of absences taken for calibration is determined by the nb.absences argument. For example, if there are 50000 pseudo-absences available and 
  nb.absences=2000, those 2000 absences will be randomly selected from the preselected 50000 pseudo-absences. This random selection is repeated NbRepPA-times 
  to constitute the runs PA1, PA2, and so forth. 
  
  NbRunEval : as already explained in the \code{\link{Initial.State}} help file, the common trend is to split the original dataset into two subsets for the calibration 
  and evaluation steps of the models construction. Here we provide the possibility to do multiple runs with this standard procedure (NbRunEval times) to have an 
  evaluation on independent-like data, whilst building the final model on 100 percent of the given data. Since march 2009, all the models produced by BIOMOD and their 
  related informations are saved on the hard disk. This is to enable the possibility of assessing the stochastic effect of the splitting procedure of the calibration step
  by being able to render projections with those models.
  
  Yweights : When NbRepPA is set to a value higher than zero, some weights are awarded automatically. For all the runs across all the species, the weights are determined
  to obtain a prevalence of 0.5, meaning that the presences will have the same importance than the absences in the calibration process of the models. 

}
\value{
  No values are returned but several objects are created : 
  
  \item{Evaluation.results.Kappa}{ }
  \item{Evaluation.results.Roc}{ }
  \item{Evaluation.results.TSS}{a list of matrices (one for each run) storing the result of the evaluation procedures for Kappa / TSS / Roc for each model : 
  the score of the model on the evaluation and calibration data, the threshold to be used for converting the probabilities into binary or filtered data, the 
  sensitivity and specificity of the model-threshold combination }
  \item{VarImportance}{ a list of matrices (one per species) containing the results of the importance of the variables analysis on the final model. The highest the value,
  the more influence the variable has on the model. A value of this 0 assumes no influence of that variable on the model. Note that this technique does not account for 
  interactions between the variables.}
  \item{Models.information}{ Crucial information for the models for rendering projections (but not of any interest for the user)}
  \item{Biomod.PA.data}{ if NbRepPA>0, it stores the available data for each species after the inner \code{\link{pseudo.abs}} run.}  
  \item{Biomod.PA.sample}{ if NbRepPA>0, it stores the data that has been selected for each run from Biomod.PA.data for calibrating the models.} 
  
  
  Additionnal objects are stored out of R in two different directories for memory storage purposes. They are created by the function directly on the root
  of your working directory set in R ("models" and "pred" directories). The first one contains the models, the second one the predictions on current data.
  Other objects will later be produced in the "pred" directory by other functions. Additionnal directories will be created by the the \code{\link{Projection}} 
  function. 
  
  The models are currently stored as objects to be read exclusively in R. To load them back (the same stands for all objects stored on the hard disk) 
  use the \code{\link{load}} function (see examples section below). 
  
  The predictions are stored as arrays (see the examples) with 4 dimensions. Please, see the "Biomod Manual.pdf" (inside the BIOMOD package) or the 
  examples below for further details.  
  %text files -> for a facilitated exportation into other softwares.
  
}


\author{ Wilfried Thuiller, Bruno Lafourcade, Robin Engler }
\seealso{ \code{\link{Initial.State}}, \code{\link{Projection}}, \code{\link{CurrentPred}}, \code{\link{PredictionBestModel}}, \code{\link{response.plot}}, \code{\link{level.plot}} }
\examples{


data(Sp.Env)
data(CoorXY)
#use fix(Sp.Env) to visualise the dataset

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,12:14], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. Here we will have 36 runs per species. This will hence take several minutes.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", GBM = TRUE, No.trees = 3000, GAM = TRUE, CTA = TRUE, CV.tree = 50, 
   ANN = TRUE, CV.ann = 2, SRE = TRUE, quant=0.05, FDA = TRUE, MARS = TRUE, RF = TRUE, NbRunEval = 2, DataSplit = 80,
   Yweights=NULL, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, KeepPredIndependent = FALSE, VarImport=5,
   NbRepPA=2, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


#view all the produced objects with ls()
ls()

#Take a look at the format in which the predictions are stored
load("pred/Pred_Sp290")
dim(Pred_Sp290)
dimnames(Pred_Sp290)[-1]

Pred_Sp290[1:20,,,1]
#please check the UserGuide for a full explanation of the structure of the array.


#You can load back a model and visualise it with the load() function :
load("models/Sp290_ANN_PA1")
Sp290_ANN_PA1

load("models/Sp290_ANN_PA1")
Sp290_ANN_PA1
str(Sp290_ANN_PA1)

#to see the results of the evaluation procedures
Evaluation.results.Kappa


#check which variable has the most impact on the models (only for the 100 percent model)
VarImportance


#This next function enables to plot the response curves of the model
load("models/Sp290_GLM_PA1")
response.plot(Sp290_GLM_PA1, Sp.Env[4:10])


#Another series of functions enable to perform further steps, see ?CurrentPred and ?PredictionBestModel.

}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }