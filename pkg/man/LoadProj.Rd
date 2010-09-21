\name{LoadProj}
\alias{LoadProj}

\title{ Loading predictions and projections produced by BIOMOD }
\description{ The function will load all the wanted data from the hardrive to R's memory for any further usage.
}

\usage{
LoadProj(Sp=1, Bin=TRUE, Filt=FALSE, which.pred="pred", raster=FALSE)
}

\arguments{
  \item{Sp}{ a vector giving the species for which you want the models to be loaded. Numbers are wanted, not names. Typing 'all' selects all the species modelled.}
  \item{Bin}{ Logical. If they exist, the binary datasets will be loaded in R. }
  \item{Filt}{ Logical. If they exist, the filtered datasets will be loaded in R.  }
  \item{which.pred}{ a vector giving the names of the prediction/projection used in Projection() for which you want the data to be loaded. 
  Names are wanted. Type "pred" for the predictions on the original data.}
  \item{raster}{ Logical. Define if the projection wanted has been done with rasters. }
}

\details{
  

}
\value{


  
}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{  
}
\seealso{ \code{\link{Models}}, \code{\link{help}} }

\examples{


data(Sp.Env)
data(CoorXY)

#This command is necessary for the run of BIOMOD as a new dataframe is produced for the Models function
Initial.State(Response=Sp.Env[,11:18], Explanatory=Sp.Env[,4:10], 
IndependentResponse=NULL, IndependentExplanatory=NULL)

#Here are done 2 PA runs and 2 repetitions for each run. This will hence take several minutes.
Models(GLM = FALSE, TypeGLM = "quad", Test = "AIC", CTA = TRUE, CV.tree = 50, ANN = TRUE, CV.ann = 2, 
   NbRunEval = 2, DataSplit = 80, Roc=TRUE, Optimized.Threshold.Roc=TRUE, Kappa=TRUE, TSS=TRUE, VarImport=5,
   NbRepPA=0, strategy="circles", coor=CoorXY, distance=2, nb.absences=1000)


Projection(Proj = Sp.Env[, 4:10], Proj.name="Current", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=T)

#Render Future projections under climate change scenario
data(Future1)
Projection(Proj = Future1[,4:10], Proj.name="Future1", GLM=TRUE, CTA=TRUE, ANN=TRUE, BinKappa=TRUE, FiltKappa=TRUE, repetition.models=T)

ls()

LoadProj(Sp='all', which.pred=c("pred", "Current", "Future1")) 
ls()


}

\keyword{ models }
