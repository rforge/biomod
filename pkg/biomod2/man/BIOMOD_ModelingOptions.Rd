% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIOMOD_ModelingOptions.R
\name{BIOMOD_ModelingOptions}
\alias{BIOMOD_ModelingOptions}
\title{Configure the modeling options for each selected model}
\usage{
BIOMOD_ModelingOptions(GLM = NULL,
                         GBM = NULL,
                         GAM = NULL,
                         CTA = NULL,
                         ANN = NULL,
                         SRE = NULL,
                         FDA = NULL,
                         MARS = NULL,
                         RF = NULL,
                         MAXENT.Phillips = NULL,
                         MAXENT.Tsuruoka = NULL)
}
\arguments{
\item{GLM}{list, GLM options}

\item{GBM}{list, GBM options}

\item{GAM}{list, GAM options}

\item{CTA}{list, CTA options}

\item{ANN}{list, ANN options}

\item{SRE}{list, SRE options}

\item{FDA}{list, FDA options}

\item{MARS}{list, MARS options}

\item{RF}{list, RF options}

\item{MAXENT.Phillips}{list, MAXENT.Phillips options}

\item{MAXENT.Tsuruoka}{list, MAXENT.Tsuruoka options}
}
\value{
A \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} object given to \code{\link[biomod2]{BIOMOD_Modeling}}
}
\description{
Parametrize and/or tune biomod's single models options.
}
\details{
The aim of this function is to allow advanced user to change some default parameters of BIOMOD inner models.  
  For each modeling technique, options can be set up.  
  
  Each argument have to be put in a list object.
  
  The best way to use this function is to print defaut models options (\code{\link{Print_Default_ModelingOptions}}) or create a default 'BIOMOD.model.option object' and print it in your console. Then copy the output, change only the required parameters, and paste it as function arguments. (see example)  
  
  Here the detailed list of modifiable parameters. They correspond to the traditional parameters that could be setted out for each modeling technique (e.g. ?GLM)
}
\section{GLM (\code{\link[stats]{glm}})}{

  
  \itemize{
    
    \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
      You can choose to either: 
        \itemize{
          \item{generate automatically the GLM formula by using the type and interaction.level arguments 
            type (default \code{'quadratic'}) : formula given to the model ('simple', 'quadratic' or 'polynomial'). 
            interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM.}
          \item{or construct specific formula}
        }}
    
    \item{\code{test} (default \code{'AIC'}) : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.}
    
    \item{\code{family} (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default.}
    
    \item{\code{control} : a list of parameters for controlling the fitting process. For glm.fit this is passed to \code{\link{glm.control}}.}
    
  }
}

\section{GBM (default \code{\link[gbm]{gbm}})}{

  
  Please refer to \code{\link[gbm]{gbm}} help file to get the meaning of this options.
  \itemize{
    \item{ \code{distribution} (default \code{'bernoulli'})}
    \item{ \code{n.trees} (default \code{2500})}
    \item{ \code{interaction.depth} (default \code{7})}
    \item{ \code{n.minobsinnode} (default \code{5})}
    \item{ \code{shrinkage} (default \code{0.001})}
    \item{ \code{bag.fraction} (default \code{0.5})}
    \item{ \code{train.fraction} (default \code{1})}
    \item{ \code{cv.folds} (default \code{3})}
    \item{ \code{keep.data} (default \code{FALSE})}
    \item{ \code{verbose} (default \code{FALSE})}
    \item{ \code{perf.method} (default \code{'cv'})}
  }
}

\section{GAM (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}})}{

  \itemize{
    
    \item{algo : either "GAM_gam" (default), "GAM_mgcv" or "BAM_mgcv" defining the chosen GAM function (see \code{\link[mgcv]{gam}}, \code{\link[gam]{gam}} resp. \code{\link[mgcv]{bam}} for more details)}
    
    \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
      You can choose to either: 
        \itemize{
          \item{generate automatically the GAM formula by using the type and interaction.level arguments 
            type : the smother used to generate the formula. Only "s_smoother" available at time. 
            interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GAM. Interaction are not considered if you choosed "GAM_gam" algo}
          \item{or construct specific formula}
        }}
    
    \item{k (default \code{-1} or \code{4}): a smooth term in a formula argument to gam (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
    \item{family (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default. }
    \item{control : see \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}}}
    
    \item{some extra "GAM_mgcv" specific options (ignored if algo = "GAM_gam")
      \itemize{
        \item{\code{method} (default \code{'GCV.Cp'})}
        \item{\code{optimizer} (default \code{c('outer','newton')})}
        \item{\code{select} (default \code{FALSE})}
        \item{\code{knots} (default \code{NULL})}
        \item{\code{paramPen} (default \code{NULL})}
      }
      
    }
  }
}

\section{CTA (\code{\link[rpart]{rpart}})}{

  
  Please refer to \code{\link[rpart]{rpart}} help file to get the meaning of the following options.
  \itemize{
    \item{\code{method} (default \code{'class'})}
    \item{\code{parms} (default \code{'default'}) : if \code{'default'}, default \pkg{rpart} parms value are kept}
    \item{\code{cost} (default \code{NULL})}
    \item{\code{control}: see \code{\link[rpart]{rpart.control}}}
  }
  
  NOTE: for method and parms, you can give a 'real' value as described in the rpart help file or 'default' that implies default \code{\link[rpart]{rpart}} values.
}

\section{ANN (\code{\link[nnet]{nnet}})}{

  
  \itemize{
    \item{\code{NbCV} (default \code{5}) : nb of cross validation to find best size and decay parameters}
    \item{\code{size}} (default \code{NULL}) : number of units in the hidden layer. If \code{NULL} then size parameter will be optimised by cross validation based on model AUC (\code{NbCv} cross validation; tested size will be the following c(2,4,6, 8) ). You can also specified a vector of size you want to test. The one giving the best model AUC will be then selected.
    \item{\code{decay}} (default \code{NULL}) : parameter for weight decay. If \code{NULL} then decay parameter will be optimised by cross validation on model AUC (\code{NbCv} cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1) ). You can also specified a vector of decay you want to test. The one giving the best model AUC will be then selected.
    \item{\code{rang} (default \code{0.1}) : Initial random weights on [-rang, rang]}
    \item{\code{maxit} (default \code{200}): maximum number of iterations.}
  }
}

\section{SRE (\code{\link[biomod2]{sre}})}{

  \itemize{
    \item{\code{quant} (default \code{0.025}): quantile of 'extreme environmental variable' removed for selection of species envelops}
  }
}

\section{FDA (\code{\link[mda]{fda}})}{

  
  Please refer to \code{\link[mda]{fda}} help file to get the meaning of these options.
  \itemize{
    \item{\code{method} (default \code{'mars'})}
    \item{\code{add_args} (default \code{NULL}) : additional arguments to \code{method} given as a list of parameters
      (corespond to the \ldots options of fda function) }
  }
}

\section{MARS (\code{\link[earth]{earth}})}{

  
  Please refer to \code{\link[earth]{earth}} help file to get the meaning of these options.
  \itemize{
    \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
      You can choose to either: 
        \itemize{
          \item{generate automatically the GLM formula by using the type and interaction.level arguments 
            type (default \code{'simple'}) : formula given to the model ('simple', 'quadratic' or 'polynomial'). 
            interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM/MARS.}
          \item{or construct specific formula}
        }}
    %    \item{\code{degree} (default \code{2})}
    \item{\code{nk}} (default \code{NULL}) : an optional integer specifying the maximum number of model terms. If NULL is given then default mars function value is used ( i.e max(21, 2 * nb_expl_var + 1) )
    \item{\code{penalty} (default \code{2})}
    \item{\code{thresh} (default \code{0.001})}
    \item{\code{nprune} (default \code{NULL})}
    \item{\code{pmethod} (default \code{"backward"})}
  }
}

\section{RF (\code{\link[randomForest]{randomForest}})}{

  
  \itemize{
    \item{\code{do.classif} (default \code{TRUE}) : if TRUE classification random.forest computed else regression random.forest will be done}
    \item{\code{ntree} (default \code{500})}
    \item{\code{mtry} (default \code{'default'})}
    \item{\code{nodesize} (default \code{5})}
    \item{\code{maxnodes} (default \code{NULL})}
  }
  
  NOTE: for mtry, you can give a 'real' value as described in randomForest help file or 'default' that implies default randomForest values
}

\section{MAXENT.Phillips (\url{www.cs.princeton.edu/~schapire/maxent/}) }{

  \itemize{
    \item{\code{path_to_maxent.jar} : character, the link to \pkg{maxent.jar} file (the working directory by default) }
    \item{\code{memory_allocated} : integer (default \code{512}), the amount of memory (in Mo) reserved for java to run MAXENT.Phillips. should be 64, 128, 256, 512, 1024, 2048... or NULL if you want to use default java memory limitation parameter.}
    \item{\code{background_data_dir} : character, path to a directory where explanatory variables are stored as ASCII files (raster format). 
      If specified MAXENT will generate it's own background data from expalantory variables rasters (as usually done in MAXENT studies). If not 
      set, then MAXENT will use the same pseudo absences than other models (generated within biomod2 at formatting step) as background data.}
    \item{\code{maximumbackground} : integer, the maximum nuber of background data to sample. This parameter will be use only if \code{background_data_dir}
      option has been set to a non default value.}
    \item{\code{maximumiterations} : integer (default \code{200}), maximum iteration done}
    \item{\code{visible} : logical (default \code{FALSE}), make the Maxent user interface visible}
    \item{\code{linear} : logical (default \code{TRUE}), allow linear features to be used}
    \item{\code{quadratic} : logical (default \code{TRUE}), allow quadratic features to be used}
    \item{\code{product} : logical (default \code{TRUE}), allow product features to be used}
    \item{\code{threshold} : logical (default \code{TRUE}), allow threshold features to be used}
    \item{\code{hinge} : logical (default \code{TRUE}), allow hinge features to be used}
    \item{\code{lq2lqptthreshold} : integer (default \code{80}), number of samples at which product and threshold features start being used}
    \item{\code{l2lqthreshold} : integer (default \code{10}), number of samples at which quadratic features start being used}
    \item{\code{hingethreshold} : integer (default \code{15}), number of samples at which hinge features start being used}
    \item{\code{beta_threshold} : numeric (default \code{-1.0}), regularization parameter to be applied to all threshold features; negative value enables automatic setting}
    \item{\code{beta_categorical} : numeric (default \code{-1.0}), regularization parameter to be applied to all categorical features; negative value enables automatic setting}
    \item{\code{beta_lqp} : numeric (default \code{-1.0}), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting}
    \item{\code{beta_hinge} : numeric (default \code{-1.0}), regularization parameter to be applied to all hinge features; negative value enables automatic setting}
    \item{\code{betamultiplier} : numeric (default \code{1}), multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.}
    \item{\code{defaultprevalence} : numeric (default \code{0.5}), default prevalence of the species: probability of presence at ordinary occurrence points}
  }
}

\section{MAXENT.Tsuruoka (\code{\link[maxent]{maxent}})}{

  
  \itemize{
    \item{\code{l1_regularizer} (default \code{0.0}): An numeric turning on L1 regularization and setting the regularization parameter. A value of 0 will disable L1 regularization}
    \item{\code{l2_regularizer} (default \code{0.0}): An numeric turning on L2 regularization and setting the regularization parameter. A value of 0 will disable L2 regularization}
    \item{\code{use_sgd} (default \code{FALSE}): A logical indicating that SGD parameter estimation should be used. Defaults to FALSE}
    \item{\code{set_heldout} (default \code{0}): An integer specifying the number of documents to hold out. Sets a held-out subset of your data to test against and prevent overfitting}
    \item{\code{verbose} (default \code{FALSE}): A logical specifying whether to provide descriptive output about the training process}
  }
  
  NOTE: if you use the \code{set_heldout} parameter then the data that will be held out will be taken in the 
  calibration data pool. It can be penilizing in case of low number of occurences dataset.
}
\examples{
  ## default BIOMOD.model.option object
  myBiomodOptions <- BIOMOD_ModelingOptions()
  
  ## print the object
  myBiomodOptions
  
  ## you can copy a part of the print, change it and custom your options 
  ## here we want to compute quadratic GLM and select best model with 'BIC' criterium
  myBiomodOptions <- BIOMOD_ModelingOptions(
    GLM = list( type = 'quadratic',
                interaction.level = 0,
                myFormula = NULL,
                test = 'BIC',
                family = 'binomial',
                control = glm.control(epsilon = 1e-08, 
                                      maxit = 1000, 
                                      trace = FALSE) ))
  
  ## check changes was done
  myBiomodOptions
  
  ##' you can prefer to establish your own GLM formula
  myBiomodOptions <- BIOMOD_ModelingOptions(
    GLM = list( myFormula = formula("Sp277 ~ bio3 + 
                    log(bio10) + poly(bio16,2) + bio19 + bio3:bio19")))
  
  ## check changes was done
  myBiomodOptions
  
  ##' you also can directly print default parameters and then follow the same processus
  Print_Default_ModelingOptions()
  
}
\author{
Damien Georges, Wilfried Thuiller
}
\keyword{models}
\keyword{options}

