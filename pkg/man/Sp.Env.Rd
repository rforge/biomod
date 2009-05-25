\name{Sp.Env}
\alias{Sp.Env}
\docType{data}
\title{ Species and Environmental data for species-climate modelling within BIOMOD}
\description{ A dataset comprising the distribution of seven species (binary records where 1=presence and 0=absence)
  and an associated set of explanatory variables for 2264 sites. This data file is an example dataset for running BIOMOD.
}
\usage{data(Sp.Env)}
\format{
  A data frame with 2264 observations on the following 18 variables.
  \describe{
    \item{\code{Idw}}{an Id line for correspondance with the coordinates data}
    \item{\code{X}}{longitude of the sites}
    \item{\code{Y}}{latitude of the sites}    
    \item{\code{Var1}}{an environmental variable}
    \item{\code{Var2}}{an environmental variable}
    \item{\code{Var3}}{an environmental variable}
    \item{\code{Var4}}{an environmental variable}
    \item{\code{Var5}}{an environmental variable}
    \item{\code{Var6}}{an environmental variable}
    \item{\code{Var7}}{an environmental variable}
    \item{\code{Sp281}}{binary record of a species geographical localisation}
    \item{\code{Sp290}}{binary record of a species geographical localisation}
    \item{\code{Sp277}}{binary record of a species geographical localisation}
    \item{\code{Sp164}}{binary record of a species geographical localisation}
    \item{\code{Sp163}}{binary record of a species geographical localisation}
    \item{\code{Sp177}}{binary record of a species geographical localisation}
    \item{\code{Sp185}}{binary record of a species geographical localisation}
    \item{\code{Sp191}}{binary record of a species geographical localisation}
  }
}
\details{
  
}
\source{
  
}
\references{
  
}
\examples{
data(Sp.Env)
fix(Sp.Env)

data(CoorXY)
levelplot(Sp.Env[,4], CoorXY[,1:2])
}
\keyword{datasets}
