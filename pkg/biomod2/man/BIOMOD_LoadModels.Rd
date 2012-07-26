\name{BIOMOD_LoadModels}
\alias{BIOMOD_LoadModels}

\title{ Load models built within BIOMOD_Modeling function}

\description{
  This function was implemented to help \pkg{biomod2} users to load individual models built during the \code{\link[biomod2]{BIOMOD_Modeling}} step.
}
\usage{
  BIOMOD_LoadModels(bm.out, ... )
}

\arguments{
  \item{bm.out}{ a \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} returned by \code{\link[biomod2]{BIOMOD_Modeling}} }
  \item{...}{ additional arguments (see details)  }
}

\details{
  This function is particulary useful when you plane to make some response plot analyses. It will load a set of individual models built within \code{\link[biomod2]{BIOMOD_Modeling}} function in your working space.
  
  Runing this function referencing only \code{bm.out} argument will lead to the loading of all models built. You can make a subselection of models you want to load using the additional arguments.
  
  \section{\bold{Additional arguments (\ldots{})} }
}

\note{
  \code{MAXENT} and \code{SRE} models are not supported yet. They will be automaticly excludes even if they are selected.
}

\value{
  A character vector filled with the loaded models names. 
}


\author{ Damien Georges }

\seealso{ \code{\link{BIOMOD_Modeling}}}
\examples{
## Not Done Yet ##
##
##
}
                 
\keyword{ models }
\keyword{ datasets }
