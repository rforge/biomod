\name{CVnnet}
\alias{CVnnet}

\title{ Inner Function }
\description{ Inner function for the BIOMOD package, not for independent usage purposes
}


\usage{
CVnnet(formula, data, truth, size = c(round(0.75 * Biomod.material[["NbVar"]]), round(0.75 * Biomod.material[["NbVar"]]), round(0.75 * Biomod.material[["NbVar"]]), Biomod.material[["NbVar"]], Biomod.material[["NbVar"]], Biomod.material[["NbVar"]]), lambda = rep(c(0.01, 0.05, 0.1), 2), nreps = 1, nifold = 3, ...)
}

\arguments{
  \item{formula}{  }
  \item{data}{  }
  \item{truth}{  }
  \item{size}{  }
  \item{lambda}{  }
  \item{nreps}{  }
  \item{nifold}{  }
  \item{\dots}{  }
}
\details{
  
}
\value{
  
}
\references{  }
\author{ Wilfried Thuiller, Bruno Lafourcade }
\note{  


}
\seealso{ \code{\link{help}}}
\examples{
}

\keyword{ models }
