\name{map.plot}
\alias{map.plot}

\title{ Plot and compare prediction maps within BIOMOD }
\description{
  This function allows direct comparisons of the predictions across models. The
  maps produced are the same than with the \code{\link{level.plot}} function 
}
\usage{
map.plot(Sp, coor, method, format.type, wanted, color.gradient = "red", Proj.name = NULL)
}

\arguments{
  \item{Sp}{ the number of the species you want, as stored in the DataBIOMOD dataframe (see examples)}
  \item{coor}{ a 2 columns matrix of the same length as data.in giving the coordinates of the data points }  
  \item{method}{ this will generate the scores of each model in the bottom right corner of the plot according to the evaluation method chosen }
  \item{format.type}{ which data to visualise, should be either 'probs', 'binary', or 'filtered' }
  \item{wanted}{ specify if you want 'prediction' (current) or 'projection' (future, different resolution) to be plotted }
  \item{color.gradient}{ available : 'red', 'grey' and 'blue' }
  \item{Proj.name}{ the name of the projection data used in Projection() }
}

\details{
The Proj.name argument needs to correspond to a projection that has been made with the Projection() function :
The function loads back the data from the projection directory (created by the Projection() function). 
  
}
\value{

}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{level.plot}} }
\examples{

data(Sp.Env)

Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM = TRUE, TypeGLM = "quad", Test = "AIC", MDA = TRUE, RF = TRUE, GBM = FALSE, GAM = FALSE, CTA = FALSE, ANN = FALSE, SRE = FALSE, MARS = FALSE, NbRunEval = 0, DataSplit = 80,
Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE)

#'Sp' stands for the column of the species in DataBIOMOD without counting the variables.
#For examples, we give Sp=2 which will plot the second species modelled with the Models function (but it is
#the 9th column of the DataBIOMOD dataframe as there are 7 variables)
data(CoorXY)
map.plot(Sp=2, coor=CoorXY, method='Roc', format.type='probs', wanted='prediction' , color.gradient='red')

#The 'method' argument will generate the scores of each model in the bottom right corner of the plot according to the
#evaluation method chosen. This particularly interesting for plotting projection to remember how performed each model

#Render Future projections under climate change scenario
data(Future1)
Projection( Proj = Future1[,2:8], Proj.name="Future1", GLM=TRUE, MDA=TRUE, RF=TRUE, BinKappa=TRUE, FiltKappa=TRUE)

#view the probabilities produced by the models
map.plot(Sp=2, coor=CoorXY, method='TSS', format.type='probs', wanted='projection' , color.gradient='red', Proj.name="Future1")
#view the projections in binary using Kappa
map.plot(Sp=2, coor=CoorXY, method='TSS', format.type='binary', wanted='projection' , color.gradient='red', Proj.name="Future1")


}

\keyword{ aplot }
