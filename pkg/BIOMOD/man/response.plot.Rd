\name{response.plot}
\alias{response.plot}

\title{ Analysis of the response curves of a model within Biomod}
\description{ Adaptation of the Evaluation Strip proposed by Elith et al.(2005). This function enables to plot the response curves of a model independently of the algorithm used for
building the model. It therefore permits a direct comparisons of models built using different statistical approaches on the same data. 
  
}
\usage{
response.plot(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE)
}

\arguments{
  \item{model}{ the model for which you want the response curves to be plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA and MARS.}
  \item{Data}{ the variables for which you want the response curves to be plotted. A data frame is wanted with one column per vriable. They have to have the same names 
   as the ones used to calibrate the model.}
  \item{show.variables}{ give in the column numbers of 'Data' for selecting the variables that are wanted for plotting }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ the size of the image in pixels if save.file is different to "no". Affects "jpeg" and "tiff" outputs only. Default if 480 pixels which is the R default. }
  \item{plot}{ if TRUE (the default) then a plot is produced. If not, an array containing predictions is returned (see details)}
}
\details{
  For building such response curves, n-1 variables are set to their median value and only the one
  of interest is varying across its whole range. The variations observed and the curve thus obtained shows
  the sensibility of the model to that specific variable. This method does therefore not account for interactions between variables.
  
  plot : if FALSE, a 3D array is returned. It contains the necessary outputs to produce the plots. This is useful to plot the response curves of all calibrated models onto the same plot (see example).  
  Array returned structure : 
    First dimension: the same number of rows than the original data
    Second dimension: 2 columns: The first one is the explanatory variable to plot, the second one, the probability of occurrence to plot
    Third dimension: The set of environmental variables for which the response.plot was asked to run.
  
  
  
}

\author{ Wilfried Thuiller}

\references{ 
Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting 
predicted responses from species distribution models. Ecological Modelling 186, 280-289.
}

\seealso{ \code{\link{Models}}  }
\examples{
\dontrun{
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
response.plot(Sp290_GLM_full, Sp.Env[,4:10], plot=TRUE)

load("models/Sp290_FDA_full")
response.plot(Sp290_FDA_full, Sp.Env[,4:10], plot=TRUE)

load("models/Sp290_RF_full")
response.plot(Sp290_RF_full, Sp.Env[,4:10], plot=TRUE)

### Save the object
GLM.RC <- response.plot(Sp290_GLM_full, Sp.Env[,4:10], plot=F)
FDA.RC <- response.plot(Sp290_FDA_full, Sp.Env[,4:10], plot=F)
RF.RC <- response.plot(Sp290_RF_full, Sp.Env[,4:10], plot=F)

### plot the response for GLM for the first variable (the 1 for the third dimension)
plot(GLM.RC[,1,1], GLM.RC[,2,1],  type="l", ylim=c(0,1), xlab=colnames(GLM.RC[,1,])[1], ylab="Probability of occurrence")
### then add the lines for the two other models with different colors
lines(FDA.RC[,1,1], FDA.RC[,2,1], col='blue')
lines(RF.RC[,1,1], RF.RC[,2,1], col='red')
library(plots)
smartlegend(x="left",y="top", inset = 0, c("GLM", "FDA", "RF"), fill = c("black", "blue", "red"))


#### the same but with all variables plotted in the same plot window. 
library(gplots)
par(mfrow=c(3,3))
for(i in 1:ncol(Sp.Env[,4:10])){
	
	plot(GLM.RC[,1,i], GLM.RC[,2,i],  type="l", ylim=c(0,1), xlab=colnames(GLM.RC[,1,])[i], ylab="Probability of occurrence")
	lines(FDA.RC[,1,i], FDA.RC[,2,i], col='blue')
	lines(RF.RC[,1,i], RF.RC[,2,i], col='red')
	smartlegend(x="left",y="top", inset = 0, c("GLM", "FDA", "RF"), fill = c("black", "blue", "red"),cex=0.6)
	}
}
}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }
