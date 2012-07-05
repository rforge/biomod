\name{BIOMOD-class}

\docType{class}

\alias{BIOMOD.formated.data-class}
\alias{BIOMOD.formated.data.PA-class}
\alias{BIOMOD.Model.Options-class}
\alias{BIOMOD.models.out-class}
\alias{BIOMOD.projection.out-class}
\alias{BIOMOD.EnsembleModeling.out-class}


\alias{BIOMOD.stored.data-class}
\alias{BIOMOD.stored.array-class}
\alias{BIOMOD.stored.formated.data-class}
\alias{BIOMOD.stored.models.options-class}


% BIOMOD.formated.data-method
\alias{show,BIOMOD.formated.data-method}
\alias{plot,BIOMOD.formated.data-method}

% BIOMOD.formated.data.PA-method
\alias{show,BIOMOD.formated.data.PA-method}
\alias{plot,BIOMOD.formated.data.PA-method}

% BIOMOD.Model.Options-method
\alias{show,BIOMOD.Model.Options-method}

% BIOMOD.models.out-method
\alias{show,BIOMOD.models.out-method}
\alias{getModelsPrediction,BIOMOD.models.out-method}
\alias{getModelsPredictionEval,BIOMOD.models.out-method}
\alias{getModelsVarImport,BIOMOD.models.out-method}
\alias{getModelsOptions,BIOMOD.models.out-method}
\alias{getModelsInputData,BIOMOD.models.out-method}
\alias{getModelsBuiltModels,BIOMOD.models.out-method}

% BIOMOD.projection.out-method
\alias{show,BIOMOD.projection.out-method}
\alias{plot,BIOMOD.projection.out-method}

% BIOMOD.EnsembleModeling.out-method
\alias{show,BIOMOD.EnsembleModeling.out-method}
\alias{getEMalgos,BIOMOD.EnsembleModeling.out-method}
\alias{getEMkeptModels,BIOMOD.EnsembleModeling.out-method}
\alias{getEMeval,BIOMOD.EnsembleModeling.out-method}
\alias{getEMbuiltModels,BIOMOD.EnsembleModeling.out-method}


\title{ BIOMOD* classes}

\description{
The BIOMOD* objects are created, used and returned by BIOMOD functions.  

  \code{BIOMOD.formated.data*}:{ user input data transformed into biomod2 format with or without Pseudo Absences sampling (see \code{\link[biomod2]{BIOMOD_FormatingData}})}  \\

  \code{BIOMOD.Model.Options}:{ individual internal models options (see \code{\link[biomod2]{BIOMOD_ModelingOptions}})}  \\

  \code{BIOMOD.models.out}:{ modeling outputs (see \code{\link[biomod2]{BIOMOD_Modeling}})} \\
  
  \code{BIOMOD.stored.*}:{ Inner object to store different kind of data into BIOMOD.models.out object} \\


}


\section{Slots}{
BIOMOD.formated.data* objects
  \describe{
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{coord}:}{ "data.frame", xy coordinates}
    \item{\code{data.species}:}{ "numeric", presences, true absences and no information species points}
    \item{\code{data.env.var}:}{ "data.frame", expalnatory variables}
    
    \item{\code{PA}:}{"data.frame", selected pseudo absences (contains also presences and true absences)}
  }
  
BIOMOD.Model.Options objects
  \describe{
    \item{\code{GLM}:}{ "list", GLM options }
    \item{\code{GBM}:}{ "list", GBM options }
    \item{\code{GAM}:}{ "list", GAM options }
    \item{\code{CTA}:}{ "list", CTA options }
    \item{\code{ANN}:}{ "list", ANN options }
    \item{\code{SRE}:}{ "list", SRE options }
    \item{\code{FDA}:}{ "list", FDA options }
    \item{\code{MARS}:}{ "list", MARS options }
    \item{\code{RF}:}{ "list", RF options }
    \item{\code{MAXENT}:}{ "list", MAXENT options }
  }

BIOMOD.models.out objects
  \describe{
    \item{\code{sp.name}:}{ "character", specie name }
    \item{\code{expl.var.names}:}{ "character", explanatory variables name }
    \item{\code{models.computed}:}{ "character", names of computed models }
    \item{\code{models.failed}:}{ "character", names of failed models }
    \item{\code{models.evaluation}:}{ "BIOMOD.stored.array", evaluations of each model computed according to selected evaluation methods }
    \item{\code{variables.importances}:}{ "BIOMOD.stored.array", models variable importances }
    \item{\code{models.prediction}:}{ "BIOMOD.stored.array", predictions of each models on their own calibration + validation dataset }
    \item{\code{formated.input.data}:}{ "BIOMOD.stored.formated.data", input data }
    \item{\code{calib.lines}:}{ "BIOMOD.stored.array", calibrations and evaluation lines selected for each run }
    \item{\code{models.options}:}{ "BIOMOD.stored.models.options", options used to build each model }
  }

BIOMOD.stored.* objects
  \describe{
    \item{\code{inMemory}:}{ "logical", is the object loaded on memory }
    \item{\code{link}:}{ "character", path to the object stored in hard drive }
    \item{\code{val}:}{ "array", "BIOMOD.Model.Options" or "BIOMOD.formated.data*", the object considered }
  }
  
}



\author{ Damien Georges }


\examples{
showClass("BIOMOD.formated.data")
}

\keyword{models}
\keyword{option}

