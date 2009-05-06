\name{Models}
\alias{Models}

\title{ Running the species-climate envelop models available in BIOMOD }
\description{
  This function allows to train various species-climate modelling algorithms. It enables also
  to run three evaluation techniques (TSS, Kappa, Roc curve) to estimate the predictive
  performance of the models produced. 
}
\usage{
Models(GLM = F, TypeGLM = "simple", Test = "AIC", GBM = F, No.trees = 2000, GAM = F, Spline = 3, CTA = F, CV.tree = 50, ANN = F, CV.ann = 5, SRE = F, Perc025 = F, Perc05 = F, MDA = F, MARS = F, RF = F, NbRunEval = 1, DataSplit = 100, Yweights = NULL, Roc = F, Optimized.Threshold.Roc = F, Kappa = F, TSS = F, KeepPredIndependent = F, VarImport = 0, NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=NULL)
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
  \item{Perc025}{ if True, the 0.025 extreme percentiles of the response data will not be used by the SRE for calibration (see \code{\link{sre}} for further details)}
  \item{Perc05}{ if True, the 0.05 extreme percentiles of the response data will not be used by the SRE for calibration (see \code{\link{sre}} for further details)}
  \item{MDA}{ type True to run the Mixture Discriminant Analysis, or False to switch it off  }
  \item{MARS}{ type True to run the Multiple Adaptive Regression Splines, or False to switch it off }
  \item{RF}{ type True to run the Random Forest, or False to switch it off }
  \item{NbRunEval}{ the number of evaluation runs to proceed before the final 100 percent data model is run (see the manual or the details section below for further explanations) }
  \item{DataSplit}{ to determine the ratio used for splitting the original database in calibration and evaluation subsets. The value given (between 0 and 100) is the percentage allocated to calibration }
  \item{Yweights}{ a N-columns matrix, N being the number of species, giving the weights that the user can set for the response variables. 
   This is similar to an index of detectability for each site, which allows users to give stronger weights to more reliable presences or absences }
  \item{Roc}{ If True, the Area Under the ROC curve score will be evaluated for all the models selected and for each species }
  \item{Optimized.Threshold.Roc}{ Type True to determine a threshold with the Roc method optimising the percentage of presences
   and absences correctly predicted (this is not automatic as the Roc curve is a threshold independent method). }
  \item{Kappa}{ If True, the Kappa score will be evaluated and a threshold produced for all the models selected and for each species }
  \item{TSS}{ If True, the True Skill Statistic score will be evaluated and a threshold produced for all the models selected and for each species }
  \item{KeepPredIndependent}{ If true, the predictions on the independent data will be saved. Otherwise, only the predictive accuracy on the evaluation set are conserved (NbRunEval) }
  \item{VarImport}{ if True, this parameter enables a direct comparison of the explanatory variable importance across models (see the manual for further explanations) }
  \item{NbRepPA}{ The number of repetitions for the selection of pseudo absences}
  \item{strategy}{ the strategy to use for selecting pseudo absences. Can be either "circles", "squares", "per", "random" or "sre". See \code{\link{pseudo.abs}} for details.}
  \item{coor}{ a two columned matrix giving the coordinates of your data points. It is needed for the "per", "circles" and "squares" strategies}
  \item{distance}{ a value giving the distance to use for the "per", "circles" and "squares" strategy of the pseudo absences selection. See \code{\link{pseudo.abs}} for details.}
  \item{nb.absences}{ the number of pseudo absences wanted to run the models with. They are randomly selected from the pool of pseudo absences available selected by the given strategy.}
}

\details{
  The dataset used is DataBIOMOD produced by the \code{\link{Initial.State}} function. A prior run of this function is therefore necessary.
  
  NbRunEval : the evaluation is as follows : as already explained in the \code{\link{Initial.State}} help file, the common trend is to split the original dataset
  into two subsets for the calibration and evaluation processes of the models. Here we provide the possibility to do a series of standard runs with this procedure
  with the NbRunEval argument to have an evaluation on independent-like data, whilst building the final model on 100 percent of the available data. This is to avoid the 
  stochastic effect of splitting data before the calibration step.

}
\value{
  No values are returned but several objects are created : 
  
  \item{Evaluation.results.Kappa}{ }
  \item{Evaluation.results.Roc}{ }
  \item{Evaluation.results.TSS}{a list of matrices (one per species) storing the result of the evaluation procedures for Kappa / TSS / Roc for each model : the score
  of the evaluation runs (if Nb RunEval>0) , the threshold to be used for converting the probabilities into binary or filtered data, the sensitivity and specificity of 
  the model-threshold combination }
  \item{VarImportance}{ a list of matrices (one per species) containing the results of the evaluation of the importance of the variables for each model selected. The highest the value,
  the more influence the variable has on the model. A value of this 0 assumes no influence of that variable on the model. Note that this technique does not account for interactions
  between the variables.}
  \item{Models.information}{ Crucial information for the models (but not of any interest for the user)}
  
  Additionnal objects are stored out of R in three different directories for memory storage purposes. Two of them are created
  by BIOMOD with the Models function on the root of your working directory set in R ("models" and "pred" directories), and the last one by
  the Projection() function ("proj" directory) . The first one contains the models, the second one the predictions on current data, the third 
  contains the projections on other sapce or time scales.
  
  The models are currently stored as objects to be read exclusively in R. To load them back (the same stands for all objects in the same format) 
  use the \code{\link{load}} function (see examples section below))
  The predictions are stored as R objects and in text files for a facilitated exportation into other softwares.

}


\references{ }


\author{ Wilfried Thuiller, Bruno Lafourcade, Robin Engler }
\note{

}
\seealso{ \code{\link{Initial.State}}, \code{\link{Projection}}, \code{\link{CurrentPred}}, \code{\link{PredictionBestModel}}, \code{\link{response.plot}}, \code{\link{level.plot}}, \code{\link{map.plot}} }
\examples{


data(Sp.Env)
#use fix(Sp.Env) to visualise the dataset

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,10:12], Explanatory=Sp.Env[,2:8], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#this will take several minutes to run
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", GBM = TRUE, No.trees = 3000, GAM = TRUE, CTA = TRUE, CV.tree = 50, 
ANN = TRUE, CV.ann = 2, SRE = TRUE, Perc025=TRUE, Perc05=FALSE, MDA = TRUE, MARS = TRUE, RF = TRUE, NbRunEval = 1, DataSplit = 80,
Yweights=NULL, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, KeepPredIndependent = FALSE, VarImport=5)

#view all the produced objects with ls()
ls()

#check which variable has the most impact on the models
VarImportance

#You can load back a model and visualise it with the load() function :
load("models/Species2_ANN")
Species2_ANN

#Other functions of the BIOMOD package are very usefull to visualise the different outputs produced. Here are only a few examples :

response.plot("GLM", Sp=1)

data(CoorXY)  #this data is needed for the next two plotting functions
              #Note that, in our case, the data loaded already has the right name for working
              #with the map.plot() function which requires the coordinate data to be named "CoorXY". 

load("pred/Pred_Species2")
x11()
level.plot(Pred_Species2[,"GLM",1,1], CoorXY)

for(i in 1:3)
map.plot(Sp=i, coor=CoorXY, method='Roc', format.type='probs', wanted='prediction' , color.gradient='red')


#Another series of functions enable to perform further steps, see ?CurrentPred and ?PredictionBestModel.


}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }