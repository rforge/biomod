\name{Sp.Env}
\alias{Sp.Env}
\docType{data}
\title{ Species and Environmental data for species-climate modelling within BIOMOD}
\description{
  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Sp.Env)}
\format{
  A data frame with 2264 observations on the following 15 variables.
  \describe{
    \item{\code{Id}}{an Id line for correspondance with the coordinates data}
    \item{\code{Var1}}{an environmental variable}
    \item{\code{Var2}}{an environmental variable}
    \item{\code{Var3}}{an environmental variable}
    \item{\code{Var4}}{an environmental variable}
    \item{\code{Var5}}{an environmental variable}
    \item{\code{Var6}}{an environmental variable}
    \item{\code{Var7}}{an environmental variable}
    \item{\code{Species1}}{binary record of a species geographical localisation}
    \item{\code{Species2}}{binary record of a species geographical localisationr}
    \item{\code{Species3}}{binary record of a species geographical localisation}
    \item{\code{Species4}}{binary record of a species geographical localisation}
    \item{\code{Species5}}{binary record of a species geographical localisation}
    \item{\code{Species6}}{binary record of a species geographical localisation}
    \item{\code{Species7}}{binary record of a species geographical localisation}
  }
}
\details{
  
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Sp.Env)
fix(Sp.Env)

data(CoorXY)
levelplot(Sp.Env[,4], CoorXY[,1:2])
}
\keyword{datasets}
