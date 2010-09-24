\name{CreateProjMatrix}
\alias{CreateProjMatrix}

\title{ Creating a dataframe compiling the wanted projections }
\description{
   Creating a dataframe for the convinience of using it within the ProDensFunc() function.
}
\usage{
CreateProjMatrix(ANN=TRUE, CTA=TRUE, GLM=TRUE, GBM=TRUE, GAM=TRUE, SRE=TRUE, FDA=TRUE, MARS=TRUE, RF=TRUE, Sp, Probs=TRUE, Bin=FALSE, Filt=FALSE, which.proj)}

\arguments{
  \item{ANN}{ if TRUE, the model will be used to render projections }
  \item{CTA}{ if TRUE, the model will be used to render projections }
  \item{GLM}{ if TRUE, the model will be used to render projections }
  \item{GBM}{ if TRUE, the model will be used to render projections }
  \item{GAM}{ if TRUE, the model will be used to render projections }
  \item{SRE}{ if TRUE, the model will be used to render projections }
  \item{FDA}{ if TRUE, the model will be used to render projections }
  \item{MARS}{ if TRUE, the model will be used to render projections }
  \item{RF}{ if TRUE, the model will be used to render projections } 
  \item{Sp}{ the name of the species for which you wish to compile the projections }
  \item{Probs}{ set it to TRUE if you want the projections in probabilities to be compiled }
  \item{Bin}{ set it to TRUE if you want the projections in binary to be compiled }
  \item{Filt}{ set it to TRUE if you want the projections in filtered values to be compiled }
  \item{which.proj}{ a collection of projection names that have been run by BIOMOD with the Projection() function }
}


\details{
  Note that the models set to TRUE need to have been run in the Models function to be able to run with Projection(). 
  Else, they simply won't be taken into account.
}

\value{
  No values are returned but 2 types of objects are produced in Rs memory.
  The names under which the objects are created follow the same logic : 
  "SpName_ProjMat" for the probabilities with an extension for bianry and/or filtered.
  This one contains the collection of projections.
  And "SpName_ProjMat_Groups" is a matrix containing various indexes for each line to
  keep track of where the projections come from.
}

\references{ }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{ }
\seealso{ \code{\link{Projection}}, \code{\link{ProbDensFunc}} }
\examples{

}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }

