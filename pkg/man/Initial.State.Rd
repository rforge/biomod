\name{Initial.State}
\alias{Initial.State}

\title{ Initialise the datasets for usage in Biomod }

\description{
  This function rearranges the datasets so they can easily be used within Biomod. A few
  other objects are stored in the memory like the number of species, the length of the datasets, etc.
  This is only a small manipualtion prior to the running of BIOMOD's models.
}
\usage{
Initial.State(Response = NULL, Explanatory = NULL, IndependentResponse = NULL, IndependentExplanatory = NULL, sp.name=NULL)
}

\arguments{
  \item{Response}{ a vector or matrix containing your species data (one column per species) in binary (ones for presences, zeros for absences) }
  \item{Explanatory}{ a matrix containing your environmental variables with which you want to explain the species distributions. Should have the same
  number of rows as in the Response matrix }
  \item{IndependentResponse}{ same format as Response. It has to have the same number of columns (and with the species sorted the same way)
  but can be of a different length. See the section details for more information. }
  \item{IndependentExplanatory}{ same as IndependentResponse but for the explanatory variables. The same variables are needed have to have the same names (if names
  were given in Explanatory)} 
  \item{sp.name}{ only necessary in the case of modelling one species at a time (i.e. a vector given in Response, give in the name of the species modelled between 
  quotation marks (for no overwritting purposes) }
}
\details{
  Giving indenpendent data is an option to evaluate more accurately the predictive performance of the models trained. 
  A model should be tested on independent data for assessing of its predictive performance, that is on data that has not been used for training it to avoid
  an overestimation of its efficiency. Generally speaking in ecology modelling, it is rare to have 2 independent datasets on the same data. Therefore, the common trend
  is to split the original dataset in 2 subsets, one for training the model and the other for testing it. But in this case, the subsets are not exactly independent.
  If you do have such real idependent data, give it in for every species to obtain a more reliable estimation of the models' accuracy.  
}
\value{
  A few objects are produced (but not returned): 
  \item{Biomod.material}{a list containing vital information for BIOMOD to run like the number of variables, the number of species modelled, their names, etc. This object
  will later be used by other functions such as \code{\link{Models}} or \code{\link{Projection}}.}
  \item{DataBIOMOD}{The dataset produced grouping the species and the variables to be used by other Biomod functions such as \code{\link{Models}}}
  \item{DataEvalBIOMOD}{An optional dataset produced by the function if independent data has been given in}

}


\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{Models}}}
\examples{

data(Sp.Env)
data(CoorXY)

#we will just take the first two species of our dataset
Initial.State(Response=Sp.Env[,11:13], Explanatory=Sp.Env[,4:10], 
		IndependentResponse=NULL, IndependentExplanatory=NULL)		


#View the data
level.plot(DataBIOMOD[,3], CoorXY, level.range=c(0,5000), title="variable 3")
x11()
level.plot(DataBIOMOD[,3], CoorXY, level.range=c(0,5000), title="variable 3", color.gradient='grey')
x11()
level.plot(DataBIOMOD[,9], CoorXY, title="distribution of 'Sp290'") 

}
                 
\keyword{ models }
\keyword{ datasets }
