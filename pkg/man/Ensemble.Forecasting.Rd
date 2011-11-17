\name{Ensemble.Forecasting}
\alias{Ensemble.Forecasting}
\alias{Ensemble.Forecasting.raster}

\title{ Make ensemble forecasts for future projection of species' distributions using data frame or rasters }
\description{
  The ensemble forecasts approach avoids the selection of one particular technique and is an answer to inter-model variability. 
  Comitee and weighted averages are calculated from the full projections produced by the models of BIOMOD, and means are calculated
  across the binary projections of the models. 
}
\usage{
Ensemble.Forecasting(ANN = TRUE, CTA = TRUE, GAM = TRUE, GBM = TRUE, GLM = TRUE, MARS = TRUE, FDA = TRUE,
    RF = TRUE, SRE = TRUE, Proj.name, weight.method, decay = 1.6, PCA.median = TRUE, binary = TRUE,
    bin.method = "Roc", Test = FALSE, repetition.models=TRUE, final.model.out=FALSE, qual.th=0, compress="xz")

Ensemble.Forecasting.raster(ANN = TRUE, CTA = TRUE, GAM = TRUE, GBM = TRUE, GLM = TRUE, MARS = TRUE,
   FDA = TRUE, RF = TRUE, SRE = TRUE, Proj.name, weight.method, decay = 1.6, PCA.median = TRUE,
   binary = TRUE, bin.method = "Roc", Test = FALSE, repetition.models = TRUE, final.model.out = FALSE, qual.th=0, compress="xz")
}

\arguments{
  \item{ANN}{ type TRUE if you want this model to be used for building the ensemble forecasts (setting a model to False is usefull if you don't like or trust a model that has high evaluation scores)}
  \item{CTA}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{GAM}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{GBM}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{GLM}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{MARS}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{FDA}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{RF}{ type TRUE if you want this model to be used for building the ensemble forecasts}
  \item{SRE}{ type TRUE if you want this model to be used for building the ensemble forecasts  }
  \item{Proj.name}{ the same name used for rendering the projections with the Projection() function. It will be used for loading the projection result files. 
   It is also the base under which the results will be stored. }
  \item{weight.method}{ the evaluation method you want to make a hierarchie of the models }
  \item{decay}{ will define the relative importance of the weights. A high value will strongly discriminate the 'good' models from the 'bad' ones (see the details section). If the value of this parameter is set to 'proportional' (entered as a string), then the attributed weights are proportional to the chosen 'weight.method' }
  \item{PCA.median}{ if True, a single model using a standard PCA will be selected to be considered as the best representation of the predictions' variability.
  Please read on the references for further details of this technique.}
  \item{binary}{ set to True if you also want the results in binary (they are produced in probabilities by default) }
  \item{bin.method}{ the method to convert the probabilities to binary values (if wanted) }
  \item{Test}{ set to True if you want the predictive performance of the ensemble method to be evaluated on current data }
  \item{repetition.models}{ set to True if you want the ensembles to be build with the repetition models taken into account }
  \item{final.model.out}{ set to True if you want the total ensemble to be build with the final models taken into account }
  \item{qual.th}{ enables to switch off the models under a certain evaluation score. This will be applied to all models on all species. }
   \item{compress}{logical or character string specifying whether saving to a named file is to use compression. FALSE corresponds to no compression, and character strings "gzip", or "xz" specify the type of compression. See ?save for more details. Default is "xz"}
}
\details{
 The decay is the ratio between a weight and the following or prior one. The formula is : W = W(-1) * decay. For example,
 with the default value of 1.6 and 4 weights wanted, the relative importance of the weights will be 1 /1.6/2.56(=1.6*1.6)/4.096(=2.56*1.6)
 from the weakest to the strongest, and gives 0.11/0.17/0.275/0.445  considering that the sum of the weights is equal to one. The lowest the decay, the smoother the differences
 between the weights enhancing a weak discrimination between models (please read the BIOMOD manual for some examples).
 
 The value 'proportional' is also possible for the decay : the weights are awarded for each method proportionaly to their evaluation scores. The advantage
 is that the discrimination is more fair than with the decay. In the latter case, close scores can strongly diverge in the weights they are awarded, when the
 proportional method will consider them as being fairly similar in prediction quality and award them a similar weight. 

Compression :
See help for the save function for more detailed explanation about compression options. Only gzip and xz are available here.
}

\value{ 
  Some objects are stored directly on the hard disk in the same directory as for the projections.

  A list is returned, named after the Proj.name argument (see examples). It contains the details of the ensemble forecast process for every species modelled and
  for every run. For each species, the information is as follows :

  \item{stats}{ contains a 2-columns matrix : if Test=TRUE, the first one gives the assumed predictive performance for each consensus method. These scores are
  obtained by applying the same ensemble computation on the current predictions as on the future forecasts, and compared with the data input for that species 
  (using the Roc evaluation method). If binary=TRUE, the second column gives the threshold obtained by each method for converting the probabilities into presence-absence predictions
  (see the manual for more details).}
  \item{weights}{ contains the weights awarded to each model selected for the ensemble forecasts (only concerns the weighted.mean).}
  \item{PCA.median}{ defines the model selected by the median PCA consensus approach. Note that nothing else is produced in regards to this method.}
}

\references{ 
Araujo, M. B. & New, M. 2007 Ensemble forecasting of species distributions. Trends in Ecology and Evolution 22, 42-47.
Araujo, M. B., Whittaker, R. J., Ladle, R. & Erhard, M. 2005 Reducing uncertainty in projections of extinction risk from climate change. Global Ecology and Biogeography 14, 529-538.
Marmion, M., Parviainen, M., Luoto, M., Heikkinen, R. K. & Thuiller, W. 2009 Evaluation of consensus methods in predictive species distribution modelling. Diversity and Distributions 15, 59-69.
Thuiller, W. 2004 Patterns and uncertainties of species' range shifts under climate change. Global Change Biology 10, 2020-2027.
}
\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}, \code{\link{Projection}}, \code{\link{ProjectionBestModel}} }
\examples{

example(Projection)

#Run the ensemble forecasting function on this scenario
Ensemble.Forecasting(Proj.name= "Future1", weight.method='Roc', PCA.median=TRUE, binary=TRUE, bin.method='Roc', Test=TRUE, repetition.models=TRUE)

#The list returned contains some of the results (see the value section for details)
load("proj.Future1/consensus_Sp277_Future1")

#this next one contains the consensus projections accross all runs for each species and consensus method
load("proj.Future1/Total_consensus_Future1")
dim(Total_consensus_Future1)
dimnames(Total_consensus_Future1)[-1]

#lets plot all these consensus results and compare them to the original data used for calibrating the models
data(CoorXY)
data <- Total_consensus_Future1

par(mfrow=c(2,7))
par(mar=c(0.6,0.6,2,0.6))
level.plot(DataBIOMOD[,8], CoorXY, show.scale=FALSE, title='sp277_original_data', cex=0.5)
level.plot(data[,1,1], CoorXY, show.scale=FALSE, title='sp277_mean', cex=0.5)
level.plot(data[,1,2], CoorXY, show.scale=FALSE, title='sp277_weighted_mean', cex=0.5)
level.plot(data[,1,3], CoorXY, show.scale=FALSE, title='sp277_median', cex=0.5)
level.plot(data[,1,4], CoorXY, show.scale=FALSE, title='sp277_Roc_mean', cex=0.5)
level.plot(data[,1,5], CoorXY, show.scale=FALSE, title='sp277_Kappa_mean', cex=0.5)
level.plot(data[,1,6], CoorXY, show.scale=FALSE, title='sp277_Roc_mean', cex=0.5)

level.plot(DataBIOMOD[,9], CoorXY, show.scale=FALSE, title='sp164_original_data', cex=0.5)
level.plot(data[,2,1], CoorXY, show.scale=FALSE, title='sp164_mean', cex=0.5)
level.plot(data[,2,2], CoorXY, show.scale=FALSE, title='sp164_weighted_mean', cex=0.5)
level.plot(data[,2,3], CoorXY, show.scale=FALSE, title='sp164_median', cex=0.5)
level.plot(data[,2,4], CoorXY, show.scale=FALSE, title='sp164_Roc_mean', cex=0.5)
level.plot(data[,2,5], CoorXY, show.scale=FALSE, title='sp164_Kappa_mean', cex=0.5)
level.plot(data[,2,6], CoorXY, show.scale=FALSE, title='sp164_TSS_mean', cex=0.5)


}

\keyword{ models }
\keyword{ optimize }