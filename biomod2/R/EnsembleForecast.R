
###############################################################################################################################
### Function description:
### *********************
### 
### The "EnsembleForecast" function returns a  
###
### -> "InMatrix": A matrix of continous projections [0:1000] where each column corresponds to a given modelling technique.
###                The columns must be names after the modelling technique whose value they contain.

### -> "Models.eval": a 3D array resulting of evalution of models quality (model used to do prediction InMatrix)
###               dim1 : Models, 
###               dim2 = "Cross.validation indepdt.data total.score Cutoff   Sensitivity Specificity" 
###               dim3 = evaluation type

### -> "EM.algo": The methods to be used for the ensemble forecast. Must be given as a vector, e.g. c("GAM","GBM","GLM","MARS","FDA","RF")
### -> "WeightMethod": must be one of the following: "Roc", "TSS" or "Kappa"
### -> "Decay": Either "proportional" or a numeric value > 0
### -> "QualThresh": models with evaluation < QualThresh will be ignored in the computing of the ensemble forecast.
### -> "PCA.median": if TRUE, then the function also returns the projections of the "PCA.median" EM method.
### -> "TruthData": optional binary vector (0 or 1) that the user can provied to evaluate the ensemble forecasts.
###                 The lenght of this vector must be equal to the number of rows of "InMatrix".
###
### Author: Robin Engler. September 2011.
###############################################################################################################################


setGeneric( "EnsembleForecast_V2", 
            def = function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
                    standardGeneric( "EnsembleForecast_V2" )
                    } )

setMethod('EnsembleForecast_V2', signature(InMatrix='matrix'),
  function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
  	
    ### Verify user input is correct
  #   if(!any(WeightMethod==c("Roc","Kappa","TSS"))) stop("\n weight.method should be one of 'Roc', 'Kappa' or 'TSS' \n") 
  	if(QualThresh<0 | QualThresh>1) stop("\n 'QualThresh' should be a value in the range [0:1[. \n") 
  	if(is.numeric(Decay)) if(Decay<0) stop("\n 'Decay' should be either 'proportional' or a numeric value > 0. \n") 
  	if(!is.numeric(Decay)) if(Decay!="proportional") stop("\n 'Decay' should be either 'proportional' or a numeric value > 0. \n") 
#   	if(!is.matrix(InMatrix) & !is.data.frame(InMatrix)) stop("\n 'InMatrix' must be a matrix or data frame object. \n") 
#   	if(is.data.frame(InMatrix)) InMatrix <- as.matrix(InMatrix) # if the input data is a dataframe, then we convert it to a matrix.
  	if(!is.null(TruthData)) if(nrow(InMatrix)!=length(as.vector(TruthData)))  stop("\n 'InMatrix' must have the same length as 'TruthData'. \n") 
  
    if(length(dim(Models.eval)) != 3) stop("\n 'Models.eval' sould be a 3d array" )	
  
  	### Assign default value to EM.algo
  	if(is.null(EM.algo) | EM.algo[1]=='all') EM.algo <- unlist(dimnames(Models.eval)[3])
  	EM.algo <- EM.algo[which(EM.algo!="SRE")] # remove SRE technique if it was selected (SRE cannot be used for ensemble forecast).            
  	
  
    cat("\n\t\tmodels selected to do EF are : ", toString(EM.algo))	
   
  	### Create list that will contain the function's output
  	EF <- list()
  	EF$EM <- rep(NA,nrow(InMatrix))
  # 	EF$EM.eval <- matrix(NA,3,4,dimnames=list(c("AUC","TSS","Kappa"),c("Eval","Threshold","Sensitivity","Specificity")))
  	EF$PCA.median <- NA
    EF$EM.name <- name
  	
  	### Verify that the columns of "InMatrix" match with the selected algorithms in EM.algo.
  	### We also use this opportunity to discard any column that contains NA values (typically because a model did not callibrate).
  	if(length(EM.algo)==0){cat("EM could not be completed as no appropriate modeling technique was selected (note that SRE cannot be used for ensemble forecast). Function returned NA values. \n"); return(EF)}
  # 	for(J in length(colnames(InMatrix)):1){
  # 		if(!any(colnames(InMatrix)[J]==EM.algo)) InMatrix <- InMatrix[,-J]
  # 	}
    InMatrix <- InMatrix[,colnames(InMatrix) %in% EM.algo]
  	InMatrix <- InMatrix[,!is.na(apply(InMatrix,2,sum))] # Delete any column that contains a value of NA
  	if(ncol(InMatrix)==0){cat("EM could not be completed. Likely reason is that no model calibrated successfully. Function returned NA values \n"); return(EF)}
  	EM.algo <- colnames(InMatrix) # re-assign EM.algo to be sure that EM.algo and columns of InMatrix are in the same order.

  	####TO MAKE WORK....
  # 	### Get the Evaluation values for each modelling modelling technique that will later be used to weight each model's 
  # 	### contribution to the ensemble forecast.
  # 	if(Biomod.material$NbRunEval>0){EvalCol <- 1                           # select column 1 if we have evaluated the models using split-sampling (i.e., Biomod.material$NbRunEval>0)
  # 	}  else{if(exists("DataEvalBIOMOD")) EvalCol <- 2 else EvalCol <- 3}    # column 2 if independent evaluation was used, and column 3 if not.
  	
   
  	EMweights <- as.numeric(Models.eval[WeightMethod,'Cross.validation',EM.algo])
  	
  	### Set to zero the weight of any model that is < QualThresh. Note that bad models can sometimes produce 
  	### evaluation values < 0 (worse than random), so these are converted to zero at the same occasion.
  	EMweights <- ifelse(EMweights<QualThresh, 0, EMweights)
  	EMweights <- ifelse(is.na(EMweights), 0, EMweights)  # any model that has "NA" as weight has its weight set to 0.
  	if(sum(EMweights)==0){cat("EM could not be completed. Likely reason is that no model calibrated successfully! \n"); return(EF)}
  	### If "Roc" was chosen as weight method, then we rescale the values in the range 0 to 1.
  	if(WeightMethod=="ROC") EMweights[EMweights != 0] <- (EMweights[EMweights != 0]-0.5)*2
  	
  	### If the used has chosen "proportional" as "Decay", then the weights are alredy OK as they are
  	### If the used has chosen a number as "Decay", then the weights are alredy OK, then weights are "decay" times decreased for each subsequent model in model quality order. 
  	if(is.numeric(Decay)){                                                            # weights are "decay" times decreased for each subsequent model in model quality order.                              
  		EMweights <- round(EMweights, 10) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
  		DecayCount <- sum(EMweights>0)
  		WOrder <- order(EMweights, decreasing=T)
  		Dweights <- EMweights
  		for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * Decay
  		#If 2 or more score are identical -> make a mean weight between the ones concerned
  		for(J in 1:length(EMweights)){
  			if(sum(EMweights[J]==EMweights)>1) Dweights[which(EMweights[J]==EMweights)] <- mean(Dweights[which(EMweights[J]==EMweights)])
  		}
  		EMweights <- Dweights
  		rm(Dweights,DecayCount,WOrder)
  	}

  	### Standardise model weights
  	EMweights <- EMweights/sum(EMweights)
  	
  	### Compute ensemble forecast
  	EF$EM <- round(as.vector(InMatrix %*% EMweights))
  	
  	### If the user has provided the TruthData, we Evaluate the EM projections
  	if(!is.null(TruthData)){
      if(length(unlist(dimnames(Models.eval)[1])) > 0){
        cat("\n\t\tEvaluating Model stuff...")
      }
      if(is.null(evalLines) | sum(evalLines)==0) evalLines <- rep(TRUE,length(as.vector(TruthData)))
      # considere here NA as absences
      TruthData <- as.vector(TruthData)
      TruthData[is.na(TruthData)] <- 0
  
      cross.validation <- sapply(unlist(dimnames(Models.eval)[1]),
                                Find.Optim.Stat,
                                Fit = EF$EM[evalLines],
                                Obs = as.vector(TruthData)[evalLines],
                                Pecision = 5)
      rownames(cross.validation) <- c("Cross.validation","Cutoff","Sensitivity", "Specificity")
      EF$EM.eval <- t(round(cross.validation,digits=3))
      rm(cross.validation)
  	}
  	
  	###If the user has asked for it, determine the model selected by the PCA consensus approach
  	if(PCA.median){
  	    if(sum(search()=="package:ade4")==0) library(ade4)  
  	    cons <- dudi.pca(InMatrix, scale=TRUE, scannf = FALSE, nf=2)
  	    EF$PCA.median <- colnames(InMatrix)[which.min(abs(cons$co[,2]))]
  	    #x11() #plotting the pca 
  	    #s.corcircle(cons$co, lab = colnames(InMatrix), full = FALSE, box = FALSE, sub=Biomod.material$species.names[i])
  	    rm(cons)
  	}    
  	
  	### release function output
  	return(EF)
  })

###############################################################################################################################
###############################################################################################################################

setMethod('EnsembleForecast_V2', signature(InMatrix='data.frame'),
  function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
    #just convert data.frame into matrix
    InMatrix <- as.matrix(InMatrix)
    return(EnsembleForecast_V2(InMatrix, Models.eval, EM.algo, WeightMethod, Decay, QualThresh, PCA.median, TruthData, evalLines, name))
  })

###############################################################################################################################
###############################################################################################################################

setMethod('EnsembleForecast_V2', signature(InMatrix='array'),
  function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
    # check Inmatrix and Models.eval dimentions
    if(length(dim(InMatrix)) != 4) stop('InMatrix must be a 4D array')
    if(length(dim(Models.eval)) != 5) stop('Models.eval must be a 5D array')
    if(dim(InMatrix)[4] != dim(Models.eval)[5] | dim(InMatrix)[3] != dim(Models.eval)[4] ) stop('Incompatible dimentions of InMatrix and Models.eval arrays')
    
    # getting dim names
    if(dim(InMatrix)[4] == 1 & length(unlist(strsplit(unlist(dimnames(InMatrix)[4]),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      dataset.names <- unlist(sapply(unlist(dimnames(InMatrix)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
    }
    
    run.eval.names <- dimnames(InMatrix)[3]
    mod.names <- paste('EF.',WeightMethod, sep='')
    
    EF.out <- list()
    # loop on each dimentions
    EF.out <- lapply(unlist(dimnames(InMatrix)[4]), function(d4){ # data.set
                lapply(unlist(dimnames(InMatrix)[3]), function(d3){ # evaluation run
                  
                  InMatrix.mod <- unlist(dimnames(InMatrix)[2])
                  if(length(grep('EF.', InMatrix.mod) > 0)){ # remove EF models if computed
                    InMatrix.mod <- InMatrix.mod[-grep('EF.', InMatrix.mod)]
                  }

                  Models.eval.mod <- unlist(dimnames(Models.eval)[3])
                  if(length(grep('EF.', Models.eval.mod) > 0)){ # remove EF models if computed
                    Models.eval.mod <- Models.eval.mod[-grep('EF.', Models.eval.mod)]
                  }
                                       
                  if( length(InMatrix.mod) > 1 & length(Models.eval.mod) > 1){
                    lapply(WeightMethod,function(w.meth){
                      return(EnsembleForecast_V2(as.matrix(InMatrix[,InMatrix.mod,d3,d4]), Models.eval[,,Models.eval.mod,d3,d4], EM.algo, w.meth, Decay, QualThresh, PCA.median, TruthData, evalLines, name))
                    })
                  } else {
                    cat('\nOnly 1 model given so EF has no sens here!')

                    return(NULL)
                  }
                })
              })
                    
    if(is.null(names(EF.out))) names(EF.out) <- unlist(dimnames(InMatrix)[4])
                    
    ef.out <- list()
    ef.out$pred <- .transform.outputs(EF.out,'EF.prediction', list(dataset.names, run.eval.names, mod.names))

    if(PCA.median){
      ef.out$PCA.median <- .transform.outputs(EF.out,'EF.PCA.median', list(dataset.names, run.eval.names, mod.names))
    }
                
    if(!is.null(TruthData)){
      ef.out$evaluation <- .transform.outputs(EF.out,'EF.evaluation', list(dataset.names, run.eval.names, mod.names))
    }
    
    rm(EF.out)            
             
    return(ef.out)
         
  })

###############################################################################################################################
###############################################################################################################################

setMethod('EnsembleForecast_V2', signature(InMatrix='RasterStack', Models.eval='matrix'),
  function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
    ### Verify user input is correct
  #   if(!any(WeightMethod==c("Roc","Kappa","TSS"))) stop("\n weight.method should be one of 'Roc', 'Kappa' or 'TSS' \n") 
  	if(QualThresh<0 | QualThresh>1) stop("\n 'QualThresh' should be a value in the range [0:1[. \n") 
  	if(is.numeric(Decay)) if(Decay<0) stop("\n 'Decay' should be either 'proportional' or a numeric value > 0. \n") 
  	if(!is.numeric(Decay)) if(Decay!="proportional") stop("\n 'Decay' should be either 'proportional' or a numeric value > 0. \n") 

    if(ncol(Models.eval) != nlayers(InMatrix)) stop("\n You must give as many layers as eval.score \n")
   
   
#   	if(!is.null(TruthData)) if(nrow(InMatrix)!=length(as.vector(TruthData)))  stop("\n 'InMatrix' must have the same length as 'TruthData'. \n") 
  
  	### Assign default value to EM.algo
      if(is.null(EM.algo) | EM.algo[1]=='all') EM.algo <- unlist(dimnames(Models.eval)[2])
  	EM.algo <- EM.algo[which(EM.algo!="SRE")] # remove SRE technique if it was selected (SRE cannot be used for ensemble forecast).            
  	
   if(length(grep("SRE",layerNames(InMatrix))) > 0){
     selected.layers <- 1:nlayers(InMatrix)
     selected.layers <- selected.layers[-grep("SRE",layerNames(InMatrix))]
     InMatrix <- InMatrix[[selected.layers]]
   }
    
    cat("\n\t\tmodels selected to do EF are : ", toString(EM.algo))	
   
  	### Create list that will contain the function's output
  	EF <- list()
  	EF$EM <- InMatrix[[1]]#reclass(InMatrix[[1]], -Inf, Inf, 0)
    EF$EM[!is.na(EF$EM[])] <- 0
  # 	EF$EM.eval <- matrix(NA,3,4,dimnames=list(c("AUC","TSS","Kappa"),c("Eval","Threshold","Sensitivity","Specificity")))
  	EF$PCA.median <- NA
    EF$EM.name <- name
  	### Verify that the columns of "InMatrix" match with the selected algorithms in EM.algo.
  	### We also use this opportunity to discard any column that contains NA values (typically because a model did not callibrate).
  	if(length(EM.algo)==0){cat("EM could not be completed as no appropriate modeling technique was selected (note that SRE cannot be used for ensemble forecast). Function returned NA values. \n"); return(EF)}
  # 	for(J in length(colnames(InMatrix)):1){
  # 		if(!any(colnames(InMatrix)[J]==EM.algo)) InMatrix <- InMatrix[,-J]
  # 	}
#     InMatrix <- InMatrix[[which(layerNames(InMatrix) %in% EM.algo)]]
#   	InMatrix <- InMatrix[,!is.na(apply(InMatrix,2,sum))] # Delete any column that contains a value of NA
#   	if(ncol(InMatrix)==0){cat("EM could not be completed. Likely reason is that no model calibrated successfully. Function returned NA values \n"); return(EF)}
#   	EM.algo <- colnames(InMatrix) # re-assign EM.algo to be sure that EM.algo and columns of InMatrix are in the same order.

  	####TO MAKE WORK....
  # 	### Get the Evaluation values for each modelling modelling technique that will later be used to weight each model's 
  # 	### contribution to the ensemble forecast.
  # 	if(Biomod.material$NbRunEval>0){EvalCol <- 1                           # select column 1 if we have evaluated the models using split-sampling (i.e., Biomod.material$NbRunEval>0)
  # 	}  else{if(exists("DataEvalBIOMOD")) EvalCol <- 2 else EvalCol <- 3}    # column 2 if independent evaluation was used, and column 3 if not.
   
  	EMweights <- as.numeric(Models.eval[which(rownames(Models.eval) %in% WeightMethod),which(colnames(Models.eval) %in% EM.algo)])
  	
  	### Set to zero the weight of any model that is < QualThresh. Note that bad models can sometimes produce 
  	### evaluation values < 0 (worse than random), so these are converted to zero at the same occasion.
  	EMweights <- ifelse(EMweights<QualThresh, 0, EMweights)
  	EMweights <- ifelse(is.na(EMweights), 0, EMweights)  # any model that has "NA" as weight has its weight set to 0.
  	if(sum(EMweights)==0){cat("EM could not be completed. Likely reason is that no model calibrated successfully! \n"); return(EF)}
  	### If "Roc" was chosen as weight method, then we rescale the values in the range 0 to 1.
#   	if(WeightMethod=="ROC") EMweights[EMweights != 0] <- (EMweights[EMweights != 0]-0.5)*2
  	
  	### If the used has chosen "proportional" as "Decay", then the weights are alredy OK as they are
  	### If the used has chosen a number as "Decay", then the weights are alredy OK, then weights are "decay" times decreased for each subsequent model in model quality order. 

  	if(is.numeric(Decay)){                                                            # weights are "decay" times decreased for each subsequent model in model quality order.                              
  		EMweights <- round(EMweights, 10) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
  		DecayCount <- sum(EMweights>0)
  		WOrder <- order(EMweights, decreasing=T)
  		Dweights <- EMweights
  		for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * Decay
  		#If 2 or more score are identical -> make a mean weight between the ones concerned
  		for(J in 1:length(EMweights)){
  			if(sum(EMweights[J]==EMweights)>1) Dweights[which(EMweights[J]==EMweights)] <- mean(Dweights[which(EMweights[J]==EMweights)])
  		}
  		EMweights <- Dweights
  		rm(Dweights,DecayCount,WOrder)
  	}

  	### Standardise model weights
  	EMweights <- EMweights/sum(EMweights)
  	
  	### Compute ensemble forecast
    for( i in 1:length(EMweights)){
      EF$EM <- EF$EM + (InMatrix[[i]] * EMweights[i])
    }
   
   layerNames(EF$EM) <- paste( paste(unlist(strsplit(layerNames(InMatrix)[1],'_'))[-length(unlist(strsplit(layerNames(InMatrix)[1],'_')))], collapse='_'),  paste("EF.",WeightMethod,sep=""),collapse='_', sep='_')
  	
  	### If the user has provided the TruthData, we Evaluate the EM projections
#   	if(!is.null(TruthData)){
#       if(length(unlist(dimnames(Models.eval)[1])) > 0){
#         cat("\n\t\tEvaluating Model stuff...")
#       }
#       if(is.null(evalLines) | sum(evalLines)==0) evalLines <- rep(TRUE,length(as.vector(TruthData)))
#       # considere here NA as absences
#       TruthData <- as.vector(TruthData)
#       TruthData[is.na(TruthData)] <- 0
#   
#       cross.validation <- sapply(unlist(dimnames(Models.eval)[1]),
#                                 Find.Optim.Stat,
#                                 Fit = EF$EM[evalLines],
#                                 Obs = as.vector(TruthData)[evalLines],
#                                 Pecision = 5)
#       rownames(cross.validation) <- c("Cross.validation","Cutoff","Sensitivity", "Specificity")
#       EF$EM.eval <- t(round(cross.validation,digits=3))
#       rm(cross.validation)
#   	}
  	
  	###If the user has asked for it, determine the model selected by the PCA consensus approach
#   	if(PCA.median){
#   	    if(sum(search()=="package:ade4")==0) library(ade4)  
#   	    cons <- dudi.pca(InMatrix, scale=TRUE, scannf = FALSE, nf=2)
#   	    EF$PCA.median <- colnames(InMatrix)[which.min(abs(cons$co[,2]))]
#   	    #x11() #plotting the pca 
#   	    #s.corcircle(cons$co, lab = colnames(InMatrix), full = FALSE, box = FALSE, sub=Biomod.material$species.names[i])
#   	    rm(cons)
#   	}    
  	
  	### release function output
  	return(EF)
  })
          
setMethod('EnsembleForecast_V2', signature(InMatrix='RasterStack', Models.eval='array'),
  function(InMatrix, Models.eval, EM.algo=NULL, WeightMethod="Roc", Decay="proportional", QualThresh=0, PCA.median=TRUE, TruthData=NULL, evalLines=NULL, name=NULL){
    # check Inmatrix and Models.eval dimentions
    if(length(dim(Models.eval)) != 5) stop('Models.eval must be a 5D array')
    
    # getting dim names
    if(dim(Models.eval)[5] == 1 & length(unlist(strsplit(layerNames(InMatrix)[1],'_'))) == 3 ){
      dataset.names <- 'AllData'
    } else{
      dataset.names <- unlist(sapply(unlist(dimnames(Models.eval)[5]), function(name){
        return(tail(unlist(strsplit(name,'_')),1))}))
    }
    
    run.eval.names <- unlist(dimnames(Models.eval)[4])
    mod.names <- paste('EF.',WeightMethod, sep='')
    
    EF.out <- list()
    # loop on each dimentions
    EF.out <- lapply(dataset.names, function(d4){ # data.set
                lapply(run.eval.names, function(d3){ # evaluation run
                  
                  InMatrix.mod <- sapply(layerNames(InMatrix), function(name){ 
                    return(tail(unlist(strsplit(name,'_'))))})
                  
                  if(length(grep('EF.', InMatrix.mod) > 0)){ # remove EF models if computed
                    InMatrix.mod <- InMatrix.mod[-grep('EF.', InMatrix.mod)]
                  }

                  Models.eval.mod <- unlist(dimnames(Models.eval)[3])
                  if(length(grep('EF.', Models.eval.mod) > 0)){ # remove EF models if computed
                    Models.eval.mod <- Models.eval.mod[-grep('EF.', Models.eval.mod)]
                  }
                                       
                  if( length(InMatrix.mod) > 1 & length(Models.eval.mod) > 1 &
                    sum(grepl(d3,layerNames(InMatrix)) & grepl(d4,layerNames(InMatrix))) > 0 ){
                    
                    lapply(WeightMethod,function(w.meth){
                        return(EnsembleForecast_V2(subset(InMatrix, 
                                                          which(grepl(d3,layerNames(InMatrix)) 
                                                                & grepl(d4,layerNames(InMatrix)) ),
                                                          drop=FALSE),
                                                          Models.eval[,'Cross.validation',Models.eval.mod,d3,d4],
#                                                    matrix(Models.eval[,'Cross.validation',Models.eval.mod,d3,d4],
#                                                           ncol = dim(Models.eval)[1],
#                                                           dimnames = list(dimnames(Models.eval)[[1]], Models.eval.mod) ),
                                                                    EM.algo, w.meth, Decay,
                                                                    QualThresh, PCA.median, TruthData, evalLines, name)['EM'])
                      
#                       return(EnsembleForecast_V2(InMatrix[[which(grepl(d3,layerNames(InMatrix)) &
#                         grepl(d4,layerNames(InMatrix)) )]], matrix(Models.eval[,'Cross.validation',Models.eval.mod,d3,d4],
#                                                                    ncol=length(Models.eval.mod),
#                                                                    dimnames = list(dimnames(Models.eval)[[1]], Models.eval.mod) ),
#                                                                     EM.algo, w.meth, Decay,
#                                                                     QualThresh, PCA.median, TruthData, evalLines, name)['EM'])
                    })
                  } else {
#                     cat('\nOnly 1 model given so EF has no sens here!')

                    return(NULL)
                  }
                })
              })
    
    EF.out <- unlist(EF.out)
    names(EF.out) <- unlist(lapply(EF.out,layerNames))
                    
    return(EF.out)
    
  })

