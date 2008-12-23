\name{response.plot}
\alias{response.plot}

\title{ Analysis of the response curves of a model within Biomod}
\description{ Adaptation of the Evalution Strip proposed by Elith et al.(2005). This function enables to plot the response curves of a model independently of the algorithm used for
building the model. It therefore permits a direct comparisons of models built using different statistical approaches on the same data. 
  
}
\usage{
response.plot(mod, Sp)
}

\arguments{
  \item{mod}{ the model for which you want the response curves to be plotted, should be one of GAM, GBM, GLM, ANN, CTA, SRE, RF, MDA MARS.}
  \item{Sp}{ the species for which you want the response curves to be plotted. A number should be given refering to the column number in the dataset given 
  in the Response argument of the \code{\link{Initial.State}} function }
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

Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#Let's run only a few models : GLM, MDA and RF. This will take a few moments to run.
#Do note that the response.plot() function works for all the models available. Feel free to switch them on or off.
Models(GLM = T, TypeGLM = "quad", Test = "AIC", MDA = T, RF = T, GBM = F, GAM = F, CTA = F, ANN = F, SRE = F, MARS = F, NbRunEval = 0, DataSplit = 80, TSS=T)

#The response plots for the first species modelled
species=1

response.plot("GLM", Sp=species)
response.plot("MDA", Sp=species)
response.plot("RF", Sp=species)	


#for the second species modelled, and so on ...
species=2

response.plot("GLM", Sp=species)
response.plot("MDA", Sp=species)
response.plot("RF", Sp=species)	


}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }
