\name{response.plot}
\alias{response.plot}

\title{ Analysis of the response curves of a model within Biomod}
\description{ Adaptation of the Evaluation Strip method proposed by Elith et al.(2005). This function enables to plot the response curves of a model independently of the algorithm used for
building the model. It therefore permits a direct comparisons of different statistical approaches on the same data. 
  
}
\usage{
response.plot(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE)
}

\arguments{
  \item{model}{ the model for which you want the response curves to be plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA and MARS.}
  \item{Data}{ the variables for which you want the response curves to be plotted. A data frame is wanted with one column per variable. They must have the same names 
   as the ones used to calibrate the model.}
  \item{show.variables}{ give in the column numbers of 'Data' for selecting the variables that are wanted to be plotted }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ the size of the image in pixels if save.file is different to "no". Affects "jpeg" and "tiff" outputs only. Default if 480 pixels which is the R default. }
  \item{plot}{ if TRUE (the default) then a plot is produced. If not, an array containing predictions is returned (see details)}
}
\details{
  For building such response curves, n-1 variables are set to their median value and only the one
  of interest is varying across its whole range. The variations observed and the curve thus obtained shows
  the sensibility of the model to that specific variable. This method does not account for interactions between variables.
  
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

\seealso{ \code{\link{BIOMOD_Modeling}}  }
\examples{
# 0. Load data & Selecting Data

# species occurances
species_occ <- read.csv(system.file("external/species/species_occ.csv",package="biomod2"))

# we consider only presences of MyocastorCoypus species
myRespName <- 'MyocastorCoypus'
myRespCoord <- species_occ[which(!is.na(species_occ[,myRespName])),c('x','y')]
myResp <- as.numeric(na.omit(species_occ[,myRespName]))

# Environemental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack(system.file("external/climat/current/bio3.grd",package="biomod2"),
                       system.file("external/climat/current/bio4.grd",package="biomod2"),
                       system.file("external/climat/current/bio7.grd",package="biomod2"),
                       system.file("external/climat/current/bio11.grd",package="biomod2"),
                       system.file("external/climat/current/bio12.grd",package="biomod2"))

# 1. Formating Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('GLM','RF','FDA'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=3, 
                                       models.eval.meth = c('TSS','ROC'),
                                       SaveObj = TRUE )


# 4. Plot response curves

#The response plots for the first species modelled 
load("MyocastorCoypus/models/MyocastorCoypus_PA1_Full_GLM")
response.plot(MyocastorCoypus_PA1_Full_GLM, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=TRUE)

load("MyocastorCoypus/models/MyocastorCoypus_PA1_Full_RF")
response.plot(MyocastorCoypus_PA1_Full_RF, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=TRUE)

load("MyocastorCoypus/models/MyocastorCoypus_PA1_Full_FDA")
response.plot(MyocastorCoypus_PA1_Full_FDA, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=TRUE)

### Save the object
GLM.RC <- response.plot(MyocastorCoypus_PA1_Full_GLM, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=FALSE)
RF.RC <- response.plot(MyocastorCoypus_PA1_Full_RF, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=FALSE)
FDA.RC <- response.plot(MyocastorCoypus_PA1_Full_FDA, getModelsInputData(myBiomomodModelOut, 'expl.var'), plot=FALSE)

### plot the response for GLM for the first variable (the 1 for the third dimension)
plot(GLM.RC[,1,1], GLM.RC[,2,1],  type="l", ylim=c(0,1), xlab=colnames(GLM.RC[,1,])[1], ylab="Probability of occurrence")
### then add the lines for the two other models with different colors
lines(FDA.RC[,1,1], FDA.RC[,2,1], col='blue')
lines(RF.RC[,1,1], RF.RC[,2,1], col='red')
# library(plots)
# smartlegend(x="left",y="top", inset = 0, c("GLM", "FDA", "RF"), fill = c("black", "blue", "red"))


#### the same but with all variables plotted in the same plot window. 
# library(gplots)
# par(mfrow=c(3,3))
# for(i in 1:ncol(getModelsInputData(myBiomomodModelOut, 'expl.var'))){
#   
#   plot(GLM.RC[,1,i], GLM.RC[,2,i],  type="l", ylim=c(0,1), xlab=colnames(GLM.RC[,1,])[i], ylab="Probability of occurrence")
# 	lines(FDA.RC[,1,i], FDA.RC[,2,i], col='blue')
# 	lines(RF.RC[,1,i], RF.RC[,2,i], col='red')
# 	smartlegend(x="left",y="top", inset = 0, c("GLM", "FDA", "RF"), fill = c("black", "blue", "red"),cex=0.6)
# 	}

}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }
