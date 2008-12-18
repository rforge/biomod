\name{Ensemble.Forecasting}
\alias{Ensemble.Forecasting}

\title{ Make ensemble forecasts for future projection of species' distributions }
\description{
  The ensemble forecasts approach avoids the selection of one particular technique and is an answer to inter-model variability. 
  Comitee and weighted averages are calculated from the full projections produced by the models of BIOMOD, and means are calculated
  across the binary projections of the models. 
}
\usage{
Ensemble.Forecasting(ANN = T, CTA = T, GAM = T, GBM = T, GLM = T, MARS = T, MDA = T, RF = T, SRE = T, Proj.name, weight.method, decay = 1.6, PCA.median = T, binary = T, bin.method = "Roc", Test = F)
}

\arguments{
  \item{ANN}{,}
  \item{CTA}{,}
  \item{GAM}{,}
  \item{GBM}{,}
  \item{GLM}{,}
  \item{MARS}{,}
  \item{MDA}{,}
  \item{RF}{ and }
  \item{SRE}{ : if TRUE, the model will be taken into account for making ensemble forecasts (setting a model to False is 
  usefull if you don't like or trust a model that has high evaluation scores) }
  \item{Proj.name}{ the name under which the results will be stored }
  \item{weight.method}{ the evaluation method you want to make a hierarchie of the models }
  \item{decay}{ will define the relative importance of the weights. A high value will strongly discriminate the 'good' models from the 'bad' ones (see the details section) }
  \item{PCA.median}{ if True, a single model using a standard PCA will be selected to be considered as the best representation of the predictions' variability.
  Please read on the references for further details of this technique.}
  \item{binary}{ set to True if you also want the results in binary (they are produced in probabilities by default) }
  \item{bin.method}{ the method to convert the probabilities to binary values (if wanted) }
  \item{Test}{ set to True if you want the predictive performance of the ensemble method to be evaluated on current data }
}
\details{
 The decay is the ratio between a weight and the following or prior one. The formula is : W = W(-1) * decay. For example,
 with the default value of 1.6 and 4 weights wanted, the relative importance of the weights will be 1 /1.6/2.56(=1.6*1.6)/4.096(=2.56*1.6)
 from the weakest to the strongest, and gives 0.11/0.17/0.275/0.445  considering that the sum of the weights is equal to one. The lowest the decay, the smoother the differences
 between the weights enhancing a weak discrimination between models (please read the BIOMOD manual for some examples).
  
}
\value{ A list is returned, named after the Proj.name argument (see examples), and several objects are produced in the "proj" folder.
 The list contains the details of the ensemble forecast process for every species modelled. For each species, the information is as follows :

  \item{`stats`}{ contains a 2-columns matrix : if Test=T, the first one gives the assumed predictive performance for each consensus method. These scores are
  obtained by applying the same ensemble computation on the current predictions as on the future forecasts, and compared with the data input for that species 
  (using the Roc evaluation method). If binary=T, the second column gives the threshold obtained by each method for converting the probabilities into presence-absence predictions
  (see the manual for more details).}
  \item{`weights`}{ contains the weights awarded to each model selected for the ensemble forecasts (only concerns the weighted.mean).}
  \item{`PCA.median`} defines the model selected by the median PCA consensus approach. Note that nothing else is produced in regards to this method.

}

\references{ 
Araújo, M. B. & New, M. 2007 Ensemble forecasting of species distributions. Trends in Ecology and Evolution 22, 42-47.
Araújo, M. B., Whittaker, R. J., Ladle, R. & Erhard, M. 2005 Reducing uncertainty in projections of extinction risk from climate change. Global Ecology and Biogeography 14, 529-538.
Marmion, M., Parviainen, M., Luoto, M., Heikkinen, R. K. & Thuiller, W. 2009 Evaluation of consensus methods in predictive species distribution modelling. Diversity and Distributions 15, 59-69.
Thuiller, W. 2004 Patterns and uncertainties of species' range shifts under climate change. Global Change Biology 10, 2020-2027.
}
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ 

}
\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{ProjectionBestModel}} }
\examples{

data(Sp.Env)
Initial.State(Response=Sp.Env[,10:11], Explanatory=Sp.Env[,2:8], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)

#This line will only run a few models
Models(GLM=T, TypeGLM="quad", Test="AIC", MDA=T, RF=T, GBM=F, GAM=F, CTA=F, ANN=F, SRE=F, MARS=F, NbRunEval=1, DataSplit=80,
Roc=T, Optimized.Threshold.Roc=T, Kappa=T, TSS=T)

#Render Future projections under climate change scenario
data(Future1)
Projection( Proj = Future1[,2:8], Proj.name="Future1", GLM=T, MDA=T, RF=T, BinKappa=T, FiltKappa=T)

#Run the ensemble forecasting function on this scenario
Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=T, binary=T, bin.method='Roc', Test=T)

#The list returned contains some of the results (see the value section for details)
#Some objects, named consensus_..., are stored in the proj folder
load("proj/consensus_Future1_mean")

#this one contains the mean projection for all the species modelled
data(CoorXY)
level.plot(consensus_Future1_mean[,1], CoorXY, title="mean projection for Species2") 



}

\keyword{ models }
\keyword{ optimize }