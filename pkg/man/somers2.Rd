\name{somers2}
\alias{somers2}

\title{ Inner Function }
\description{ Computes Somers' Dxy rank correlation between a variable x and a binary (0-1) variable y, and the corresponding receiver operating characteristic curve area c. Note that Dxy = 2(c-0.5). somers allows for a weights variable, which specifies frequencies to associate with each observation.
}
\usage{
somers2(x, y, weights=NULL, normwt=FALSE, na.rm=TRUE)
}

\arguments{
  \item{x}{ typically a predictor variable. NAs are allowed. }
  \item{y}{ a numeric outcome variable coded 0-1. NAs are allowed. }
  \item{weights}{ a numeric vector of observation weights (usually frequencies). Omit or specify a zero-length vector to do an unweighted analysis. }
  \item{normwt}{ set to TRUE to make weights sum to the actual number of non-missing observations. }
  \item{na.rm}{ set to FALSE to suppress checking for NAs.  }
}


\details{
The rcorr.cens function, which although slower than somers2 for large sample sizes, can also be used to obtain Dxy for non-censored binary y, and it has the advantage of computing the standard deviation of the correlation index.  
}
\value{
a vector with the named elements C, Dxy, n (number of non-missing pairs), and Missing. Uses the formula C = (mean(rank(x)[y == 1]) - (n1 + 1)/2)/(n - n1), where n1 is the frequency of y=1.   
}
\references{  }
\author{ 
Frank Harrell
Department of Biostatistics
Vanderbilt University School of Medicine
f.harrell@vanderbilt.edu 
}
\note{  
See also rcorr.cens, rank, wtd.rank (functions of the Hmisc package)
}
\seealso{}
\examples{
set.seed(1)
predicted <- runif(200)
dead      <- sample(0:1, 200, TRUE)
roc.area <- somers2(predicted, dead)["C"]

# Check weights
x <- 1:6
y <- c(0,0,1,0,1,1)
f <- c(3,2,2,3,2,1)
somers2(x, y)
somers2(rep(x, f), rep(y, f))
somers2(x, y, f)
}

\keyword{ models }
