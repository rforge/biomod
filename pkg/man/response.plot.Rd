\name{response.plot}
\alias{response.plot}

\title{ Analysis of the response curves of a model within Biomod}
\description{ Adaptation of the Evalution Strip proposed by Elith et al.(2005). This function enables to plot the response curves of a model independently of the algorithm used for
building the model. It therefore permits a direct comparisons of models built using different statistical approaches on the same data. 
  
}
\usage{
response.plot(model, Data, save.file="no", name="response_curve")
}

\arguments{
  \item{model}{ the model for which you want the response curves to be plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA and MARS.}
  \item{Data}{ the variables for which you want the response curves to be plotted. A data frame is wanted with one column per vriable. They have to have the same names 
   as the ones used to calibrate the model.}
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
}
\details{
  For building such response curves, n-1 variables are set to their median value and only the one
  of interest is varying accross its whole range. The variations observed and the curve thus obtained shows
  the sensibility of the model to that specific variable. This method does therefore not account for interactions between variables.
}
\value{
}

\author{ Wilfried Thuiller, Bruno Lafourcade }

\references{ 
Elith, J., Ferrier, S., Huettmann, F. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting 
predicted responses from species distribution models. Ecological Modelling 186, 280-289.
}

\seealso{ \code{\link{Models}}  }
\examples{

data(Sp.Env)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,12:13], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Let's run only a few models : GLM, FDA and RF. This will take a few moments to run.
#Do note that the response.plot() function works for all the models available. Feel free to switch them on or off.
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", FDA = TRUE, CV.tree = 50, RF = TRUE, CV.ann = 2, 
   NbRunEval = 0, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=1000)


#The response plots for the first species modelled 
load("models/Sp290_GLM_full")
response.plot(Sp290_GLM_full, Sp.Env[,4:10])

load("models/Sp290_FDA_full")
response.plot(Sp290_FDA_full, Sp.Env[,4:10])

load("models/Sp290_RF_full")
response.plot(Sp290_RF_full, Sp.Env[,4:10])



#for the second species modelled
load("models/Sp277_GLM_full")
response.plot(Sp277_GLM_full, Sp.Env[,4:10])

load("models/Sp277_FDA_full")
response.plot(Sp277_FDA_full, Sp.Env[,4:10])

load("models/Sp277_RF_full")
response.plot(Sp277_RF_full, Sp.Env[,4:10])

}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }
