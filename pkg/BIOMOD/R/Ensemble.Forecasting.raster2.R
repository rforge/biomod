`Ensemble.Forecasting.raster2` <-
function(ANN=TRUE,CTA=TRUE,GAM=TRUE,GBM=TRUE,GLM=TRUE,MARS=TRUE,FDA=TRUE,RF=TRUE,SRE=TRUE, Proj.name, weight.method, decay=1.6, PCA.median=TRUE, binary=TRUE, bin.method='Roc', Test=FALSE, repetition.models=TRUE, final.model.out=FALSE, qual.th=0, compress="xz", return.stack=TRUE)
{
  library('raster', quietly=TRUE)
  ### TMP PARAMS THAT WILL HAVE TO ADD IN FCT ARGS...
  ef.assembling <- c('')
  
  ### check Parameters
  
  args <- .check.Ensemble.Forecasting.raster(ANN,CTA,GAM,GBM,GLM,MARS,FDA,RF,SRE, Proj.name, weight.method, decay, PCA.median, binary, bin.method, Test, repetition.models, final.model.out, qual.th, compress, return.stack)

  compress <- args$compress
  repetition.models <- args$repetition.models
  ProjRunEval <- args$ProjRunEval
  final.model.out <- args$final.model.out
  ens.choice <- args$ens.choice
  work.with.stack <- args$work.with.stack
  ef.algo <- args$ef.algo
  prob.mean.weigth.decay <- args$prob.mean.weigth.decay
  return.stack <- args$return.stack

  rm('args')
    
  #create storing list of outputs
  list.out <- vector('list', Biomod.material$NbSpecies)
  names(list.out) <- Biomod.material$species.names                                                                                                                                 
    
  #--------- start species loop ---------       
  for(i in 1:Biomod.material$NbSpecies){
    work.space.species.loop <- ls()
    
    cat('\n\n-=-=-= ',Biomod.material$species.names[i])
    
    # 1. get all attributs we need to compute EF -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    
    # considering the number of PA runs and reps that were actually done 
    NbPA <- max(Biomod.material$NbRepPA, 1)

    nbrep <- 1
    if(repetition.models) nbrep <- nbrep + Biomod.material$NbRunEval
#     nbrep <- nbrep + Biomod.material$NbRunEval
    
    # getting the scores and the thresholds of chosen eval method
    if(!is.null(Biomod.material$Independent.data.set)){
      if(Biomod.material$Independent.data.set){
        scoresCol = "indepdt.data"
      } else if(nbrep > 1 ){
        scoresCol = "Cross.validation"
      } else { scoresCol = "total.score" }
      
    } else if(nbrep > 1 ){
      scoresCol = "Cross.validation"
    } else { scoresCol = "total.score" }
    
    scores.Roc <- scores.TSS <- scores.Kappa <- NULL
    
    if(!is.na(Evaluation.results.Roc[1])){
      scores.Roc <- array(NA, dim=c(2,sum(ens.choice), nbrep, NbPA))
      
      if(nbrep<=1){
        dim2.names <- 'allData'
      } else {
        dim2.names <- c('allData', paste('rep',1:(nbrep-1), sep=""))
      }
      
      if(Biomod.material$NbRepPA<1){
        dim3.names <- 'full'
      } else {
        dim3.names <- paste("PA", 1:NbPA, sep="")
      }

      dimnames(scores.Roc) <- list(c('score', 'threshold'),
                                   names(ens.choice)[which(ens.choice == TRUE)],
                                   dim2.names,
                                   dim3.names)
                                   
      for(mod.name in grep(Biomod.material$species.names[i] ,names(Evaluation.results.Roc), value=TRUE ) ){
        # do not take repetition into acount if they are not selected
        if(!repetition.models){
          if(grepl('_rep', mod.name)){
            next
          }
        }
        
        # get information on run with table name
        mod.name.split <- unlist(strsplit(mod.name,'_'))
        if(grepl("full", tail(mod.name.split,1)) | grepl("PA", tail(mod.name.split,1))){
          mod.name.split <- c(mod.name.split, "allData")
        }

        # fill corectly our array
        # don't take the cv score for models calibrates with all data
        if(tail(mod.name.split,1) == "allData" & scoresCol == "Cross.validation"){
          scores.Roc['score', , "allData", tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Roc[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "total.score"])
        } else {
          scores.Roc['score', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Roc[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], scoresCol])
        }
        scores.Roc['threshold', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Roc[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "Cutoff"])
      }
    }
    
    if(!is.na(Evaluation.results.TSS[1])){
      scores.TSS <- array(NA, dim=c(2,sum(ens.choice), nbrep, NbPA))
      
      if(nbrep<=1){
        dim2.names <- 'allData'
      } else {
        dim2.names <- c('allData', paste('rep',1:(nbrep-1), sep=""))
      }
      
      if(Biomod.material$NbRepPA<1){
        dim3.names <- 'full'
      } else {
        dim3.names <- paste("PA", 1:NbPA, sep="")
      }
      
      dimnames(scores.TSS) <- list(c('score', 'threshold'),
                                   names(ens.choice)[which(ens.choice == TRUE)],
                                   dim2.names,
                                   dim3.names)
      
      for(mod.name in grep(Biomod.material$species.names[i] ,names(Evaluation.results.TSS), value=TRUE ) ){
        # do not take repetition into acount if they are not selected
        if(!repetition.models){
          if(grepl('_rep', mod.name)){
            next
          }
        }
        
        # get information on run with table name
        mod.name.split <- unlist(strsplit(mod.name,'_'))
        if(grepl("full", tail(mod.name.split,1)) | grepl("PA", tail(mod.name.split,1))){
          mod.name.split <- c(mod.name.split, "allData")
        }
        
        # fill corectly our array
        # don't take the cv score for models calibrates with all data
        if(tail(mod.name.split,1) == "allData" & scoresCol == "Cross.validation"){
          scores.TSS['score', , "allData", tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.TSS[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "total.score"])
        } else {
          scores.TSS['score', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.TSS[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], scoresCol])
        }
        
        scores.TSS['threshold', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.TSS[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "Cutoff"])
      }      

    }
    
    if(!is.na(Evaluation.results.Kappa[1])){
      scores.Kappa <- array(NA, dim=c(2,sum(ens.choice), nbrep, NbPA))
      
      if(nbrep<=1){
        dim2.names <- 'allData'
      } else {
        dim2.names <- c('allData', paste('rep',1:(nbrep-1), sep=""))
      }
      
      if(Biomod.material$NbRepPA<1){
        dim3.names <- 'full'
      } else {
        dim3.names <- paste("PA", 1:NbPA, sep="")
      }
      
      dimnames(scores.Kappa) <- list(c('score', 'threshold'),
                                   names(ens.choice)[which(ens.choice == TRUE)],
                                   dim2.names,
                                   dim3.names)
      
      for(mod.name in grep(Biomod.material$species.names[i] ,names(Evaluation.results.Kappa), value=TRUE ) ){
        # do not take repetition into acount if they are not selected
        if(!repetition.models){
          if(grepl('_rep', mod.name)){
            next
          }
        }
        
        # get information on run with table name
        mod.name.split <- unlist(strsplit(mod.name,'_'))
        if(grepl("full", tail(mod.name.split,1)) | grepl("PA", tail(mod.name.split,1))){
          mod.name.split <- c(mod.name.split, "allData")
        }
        
        # fill corectly our array
        # don't take the cv score for models calibrates with all data
        if(tail(mod.name.split,1) == "allData" & scoresCol == "Cross.validation"){
          scores.Kappa['score', , "allData", tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Kappa[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "total.score"])
        } else {
          scores.Kappa['score', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Kappa[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], scoresCol])
        }
        
        scores.Kappa['threshold', , tail(mod.name.split,1)[1], tail(mod.name.split,2)[1] ] <- as.numeric(Evaluation.results.Kappa[[mod.name]][names(ens.choice)[which(ens.choice == TRUE)], "Cutoff"])
      }
    }


    ## get bad models names based on wheight meth
    models.to.keep <- models.to.remove <- c()
    models.to.keep.thresh <- c()
    models.to.keep.scores <- models.to.remove.scores <- c()
    models.to.keep.thresh.Roc <- models.to.keep.thresh.TSS <- models.to.keep.thresh.Kappa <- c()
    
    if(qual.th < 1){
      cat("\n   > Removing models having a", weight.method, "<", qual.th)
      scores.tmp <- get(paste("scores.", weight.method, sep=""))
      
      for(d2 in dimnames(scores.tmp)[[2]]){
        for(d3 in dimnames(scores.tmp)[[3]]){
          for(d4 in dimnames(scores.tmp)[[4]]){
            if(!is.na(scores.tmp['score',d2,d3,d4])){
              if(scores.tmp['score',d2,d3,d4] > qual.th){
                models.to.keep <- c(models.to.keep, paste(d4,d3,d2, sep="_"))
                models.to.keep.thresh <- c(models.to.keep.thresh, scores.tmp['threshold',d2,d3,d4])
                models.to.keep.scores <- c(models.to.keep.scores, scores.tmp['score',d2,d3,d4])
                ## keep also model thresh for other metrics
                if(!is.null(scores.Roc)){
                  # define an artificial threshold for SRE
                  if(d2 == 'SRE'){
                    models.to.keep.thresh.Roc <- c(models.to.keep.thresh.Roc, 500)
                  } else{
                    models.to.keep.thresh.Roc <- c(models.to.keep.thresh.Roc, scores.Roc['threshold',d2,d3,d4])
                  }
                }
                if(!is.null(scores.TSS)){
                  models.to.keep.thresh.TSS <- c(models.to.keep.thresh.TSS, scores.TSS['threshold',d2,d3,d4])
                }
                if(!is.null(scores.Kappa)){
                  models.to.keep.thresh.Kappa <- c(models.to.keep.thresh.Kappa, scores.Kappa['threshold',d2,d3,d4])
                }
              
              } else{
                models.to.remove <- c(models.to.remove, paste(d4,d3,d2, sep="_"))
                models.to.remove.scores <- c(models.to.remove.scores, scores.tmp['score',d2,d3,d4])
              }     
            }
          }
        }
      }
      rm('scores.tmp')
      if(length(models.to.keep)>0){
#         cat("\n      ", toString(models.to.keep), "models will be kept to compute EF")
        cat("\n")
        cat("\n   Models kept :", toString(models.to.keep))
        cat("\n")
        cat("\n   Models removed :", toString(models.to.remove))
        cat("\n")
      } else{
        cat("\n   ! all models were excluded for this treshold.. EF not done for this specie !")
        break
      }
    }

    # 2. creating outputs -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    
    # 3. doing EF -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    # get proj names we will work with...
    proj.files <- grep(Biomod.material$species.names[i] ,list.files(path=paste(getwd(),"/proj.", Proj.name, sep=""), pattern=".raster" ), value=TRUE)

    # remove Bin and Filt
    if(sum(grepl("Bin", proj.files)) > 0){
      proj.files <- proj.files[-grep("Bin", proj.files)]
    }
    if(sum(grepl("Filt", proj.files)) > 0){
      proj.files <- proj.files[-grep("Filt", proj.files)]
    }
    
    # if projection are not stacked, we can keep only the selected ones else we will have to make the selection later
    if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
      proj.files.kept <- paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_", models.to.keep, ".raster",sep="")
      # remove the part "_allData" if exists
      proj.files.kept <- unlist( lapply(proj.files.kept, function(x){
        if(grepl("_allData", x)){
          return(sub("_allData", "", x))
        } else{
          return(x)
        }
      }) )
      
      # check validity of kept names
      if(sum(!(proj.files.kept %in% proj.files)) == 0){
        proj.files <- proj.files.kept
      } else{
        stop("Internal model selection issue")
      }
    }
    

    if(work.with.stack){
      # make a enormous stack containing all projections
      eval(parse(text= paste("BigStack <- raster::stack(",paste("get(load('", rep(getwd(), length(proj.files)),"/proj.", Proj.name, "/" , proj.files ,"'))", collapse = ", ", sep="") ,")", sep="") ))
      
      rm(list=proj.files)
      
      if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
        layer.names <- proj.files
        # remove all what is useless in layer.names
        layer.names <- gsub(paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_",sep=""), "", layer.names)
        layer.names <- gsub(".raster", "", layer.names)
        
        # add allData for run made with all data
        layer.names <- unlist(lapply(layer.names, function(x){
          if(length(unlist(strsplit(x,"_"))) < 3){
            x.split <- unlist(strsplit(x,"_"))
            return(paste(x.split[1],"_allData_",x.split[2], sep=""))
          } else { return(x) }
        }))
        
        layerNames(BigStack) <- layer.names
      } else{
        # remove all what is useless in layer.names
        layer.names <- layerNames(BigStack)
        layer.names <- unlist(lapply(layer.names, function(x){
          return(head(unlist(strsplit(x, ".", fixed=T)),1))
        }))
          
        # add allData for run made with all data
        layer.names <- unlist(lapply(layer.names, function(x){
          if(length(unlist(strsplit(x,"_"))) < 2){
            return(paste(x,"_allData", sep=""))
          } else { return(x) }
        }))
          
        mod <- unique(unlist(lapply(proj.files,function(xx){
          return(tail(unlist(strsplit(gsub(".raster","", xx), "_")),1))
        })))

        layer.names <- paste(layer.names, rep(mod,each = as.integer(length(layer.names) / length(mod))), sep="_")

        layerNames(BigStack) <- layer.names
                                     
        # remove bad models
        if(length(models.to.keep) > 0){
          BigStack <- raster::subset(BigStack, models.to.keep, drop=FALSE)
        }
      }

    }

    for(algo in ef.algo){
      cat("\n  -", algo)
      
      ## 1. Probabilities mean =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'prob.mean'){
        if(exists('BigStack')){
          ef.prob.mean <- round(mean(BigStack))
        } else{
          if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
            # initialisation of raster
            ef.prob.mean <- get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[1], sep="")))
            rm(list=proj.files[1])
            if(length(proj.files > 1)){
              for(pf in proj.files){
                ef.prob.mean <- ef.prob.mean+ get(load(paste(getwd(), "/proj.", Proj.name, "/" , pf, sep="")))
                rm(list=pf)
              }
              ef.prob.mean <- round( ef.prob.mean / length(proj.files) )             
            }
          } else{
            stop("!!! NOT SUPPORTED YET \nPlease redo your projections with arg stack.out=FALSE ")
          }
        }
          
        ef.prob.mean.thresh <- mean(models.to.keep.thresh, na.rm=T) 
        
      }
      
      ## 2. Probabilities median =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'median'){
        # works only if we can construct the BigStack
        if(exists('BigStack')){
          ef.median <- round(calc(BigStack, median))
          ef.median.thresh <- median(models.to.keep.thresh, na.rm=T)
        } else{
          cat("\n   !!! Not Computed !!!")
        }
      }

      ## 3. Probabilities weighted mean -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'prob.mean.weighted'){
        
        models.kept.tmp <- models.to.keep
        models.kept.scores.tmp <- models.to.keep.scores
        models.kept.thresh.tmp <- models.to.keep.thresh
        proj.files.tmp <- proj.files
        
        # remove SRE from models if Roc
        if(sum(grepl("SRE", models.to.keep) > 0)){
          cat("\n   ! SRE algo have been removed for this EM algo !")
          models.kept.tmp <- models.to.keep[- grep("SRE", models.to.keep)]
          models.kept.scores.tmp <- models.to.keep.scores[- grep("SRE", models.to.keep)]
          models.kept.thresh.tmp <- models.to.keep.thresh[- grep("SRE", models.to.keep)]
        }
           
        # weights are "decay" times decreased for each subsequent model in model quality order.                              
        models.kept.scores.tmp <- round(models.kept.scores.tmp, 10) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
        DecayCount <- sum(models.kept.scores.tmp>0)
    		WOrder <- order(models.kept.scores.tmp, decreasing=T)
    		Dweights <- models.kept.scores.tmp
    		for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * prob.mean.weigth.decay
    		#If 2 or more score are identical -> make a mean weight between the ones concerned
    		for(J in 1:length(models.kept.scores.tmp)){
    			if(sum(models.kept.scores.tmp[J]==models.kept.scores.tmp)>1) Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)] <- mean(Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)])
    		}
    		models.kept.scores.tmp <- Dweights
    		rm(Dweights,DecayCount,WOrder)
        
        ### Standardise model weights
        models.kept.scores.tmp <- models.kept.scores.tmp/sum(models.kept.scores.tmp)
        
      	### Compute ensemble forecast
        if(exists('BigStack')){
          if(length(models.kept.tmp) > 1){
            ef.prob.mean.weighted <- round(sum(raster::subset(BigStack, models.kept.tmp) * models.kept.scores.tmp))
          } else{
            ef.prob.mean.weighted <- round(raster::subset(BigStack, models.kept.tmp) * models.kept.scores.tmp)
          }
          
        } else{
          if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
            if(sum(grepl("SRE", proj.files.tmp) > 0)){
              proj.files.tmp <- proj.files.tmp[ - grep("SRE", proj.files.tmp)]
            }
            
            # initialisation of raster
            ef.prob.mean.weighted <- get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files.tmp[1], sep=""))) *  models.kept.scores.tmp[1]
            rm(list=proj.files.tmp[1])
            if(length(proj.files.tmp > 1)){
              for(pf in proj.files.tmp[-1]){
                ef.prob.mean.weighted <- ef.prob.mean.weighted + get(load(paste(getwd(), "/proj.", Proj.name, "/" , pf, sep=""))) * models.kept.scores.tmp[which(proj.files.tmp == pf)]
                rm(list=pf)
              }
#               ef.prob.mean.weighted <- round( ef.prob.mean.weighted / length(proj.files.tmp) )             
            }
          } else{
            stop("!!! NOT SUPPORTED YET \nPlease redo your projections with arg stack.out=FALSE ")
          }
        }
        ef.prob.mean.weighted.thresh <- sum( models.kept.thresh.tmp * models.kept.scores.tmp , na.rm=T)
 
      }
      
      ## 4. Roc comitee averaging -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'Roc.mean'){
        if(weight.method != 'Roc'){
          cat("\n      ! model selection was based on", weight.method, "score, that can lead to strange results !")
        }
        
        # check if we are working with bigStack
        if(exists('BigStack')){      
          ef.Roc.mean <- mean(BinaryTransformation(BigStack, models.to.keep.thresh.Roc))*1000
        } else{
          if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
            # initialisation of raster
            ef.Roc.mean <- BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[1], sep=""))),models.to.keep.thresh.Roc[1])
            rm(list=proj.files[1])
            if(length(proj.files) > 1){
              for(pfID in 2:length(proj.files)){
                ef.Roc.mean <- ef.Roc.mean + BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[pfID], sep=""))), models.to.keep.thresh.Roc[pfID])
                rm(list=proj.files[pfID])
              }
              ef.Roc.mean <- round( ef.Roc.mean * 1000 / length(proj.files) )             
            }
          } else{
            stop("!!! NOT SUPPORTED YET \nPlease redo your projections with arg stack.out=FALSE ")
          }
        }
          
        ef.Roc.mean.thresh <- 500
        
      }
        
      ## 5. TSS comitee averaging -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'TSS.mean'){
        if(weight.method != 'TSS'){
          cat("\n      ! model selection was based on", weight.method, "score, that can lead to strange results !")
        }
        
        # check if we are working with bigStack
        if(exists('BigStack')){
          ef.TSS.mean <- mean(BinaryTransformation(BigStack, models.to.keep.thresh.TSS))*1000
        } else{
          if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
            # initialisation of raster
            ef.TSS.mean <- BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[1], sep=""))),models.to.keep.thresh.TSS[1])
            rm(list=proj.files[1])
            if(length(proj.files) > 1){
              for(pfID in 2:length(proj.files)){
                ef.TSS.mean <- ef.TSS.mean + BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[pfID], sep=""))), models.to.keep.thresh.TSS[pfID])
                rm(list=proj.files[pfID])
              }
              ef.TSS.mean <- round( ef.TSS.mean * 1000 / length(proj.files) )             
            }
          } else{
            stop("!!! NOT SUPPORTED YET \nPlease redo your projections with arg stack.out=FALSE ")
          }
        }
          
        ef.TSS.mean.thresh <- 500
 
      }
        
      ## 6. Kappa comitee averaging -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(algo == 'Kappa.mean'){
        if(weight.method != 'Kappa'){
          cat("\n      ! model selection was based on", weight.method, "score, that can lead to strange results !")
        }
        
        # check if we are working with bigStack
        if(exists('BigStack')){
          ef.Kappa.mean <- mean(BinaryTransformation(BigStack, models.to.keep.thresh.Kappa))*1000
        } else{
          if(!Biomod.material[[paste("proj.", Proj.name , ".stack", sep="")]]){
            # initialisation of raster
            ef.Kappa.mean <- BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[1], sep=""))),models.to.keep.thresh.Kappa[1])
            rm(list=proj.files[1])
            if(length(proj.files) > 1){
              for(pfID in 2:length(proj.files)){
                ef.Kappa.mean <- ef.Kappa.mean + BinaryTransformation(get(load(paste(getwd(), "/proj.", Proj.name, "/" , proj.files[pfID], sep=""))), models.to.keep.thresh.Kappa[pfID])
                rm(list=proj.files[pfID])
              }
              ef.Kappa.mean <- round( ef.Kappa.mean * 1000 / length(proj.files) )             
            }
          } else{
            stop("!!! NOT SUPPORTED YET \nPlease redo your projections with arg stack.out=FALSE ")
          }
        }
          
        ef.Kappa.mean.thresh <- 500
      
      }
 
      
    }
    
    # 7. Assembling all computed models -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    ef.potential <- c('ef.prob.mean','ef.prob.mean.weighted','ef.median','ef.Roc.mean','ef.TSS.mean','ef.Kappa.mean')
    ef.computed <- ef.potential[unlist(lapply(ef.potential, exists, envir = environment()))]

    eval(parse(text = paste("ef.proj <- raster::stack(", toString(ef.computed), ")", sep="")))
    layerNames(ef.proj) <- ef.computed

    eval(parse(text = paste("ef.proj.thresh <- data.frame(", paste(ef.computed," = ", ef.computed, ".thresh" , sep="", collapse=","), ")", sep="")))

    rm(list=ef.computed)
    
    eval(parse(text=paste("EF.",Biomod.material$species.names[i], " <- list(projections = ef.proj, thresholds = ef.proj.thresh, models.kept = models.to.keep)", sep="")))

#     eval(parse(text = paste("list.out$", Biomod.material$species.names[i], " <- list(projections = ef.proj, thresholds = ef.proj.thresh)")))
        
    # 7. save specie EF
    eval(parse(text = paste("save(EF.",Biomod.material$species.names[i],", file = '",getwd(),"/proj.", Proj.name , "/EF.",Biomod.material$species.names[i] ,"')", sep="") ))
        
    if(return.stack){
      eval(parse(text = paste("list.out$", Biomod.material$species.names[i], " <- EF.",Biomod.material$species.names[i], sep="")))
    } else {
      eval(parse(text = paste("list.out$", Biomod.material$species.names[i], " <- list(projection ='",getwd(),"/proj.", Proj.name , "/EF.",Biomod.material$species.names[i],"', thresholds = ef.proj.thresh, models.kept = models.to.keep)", sep="")))
    }

    
      
    # 8. remove usless objects
    rm(list = ls()[- which(ls() %in% work.space.species.loop)])

    
  }
  return(list.out)
}
 
####################################################################################################
####################################################################################################
####################################################################################################

.check.Ensemble.Forecasting.raster <- function(ANN=TRUE,CTA=TRUE,GAM=TRUE,GBM=TRUE,GLM=TRUE,MARS=TRUE,FDA=TRUE,RF=TRUE,SRE=TRUE, Proj.name, weight.method, decay=1.6, PCA.median=TRUE, binary=TRUE, bin.method='Roc', Test=FALSE, repetition.models=TRUE, final.model.out=FALSE, qual.th=0, compress="xz", return.stack=TRUE){
  
  # 1. eval method params
  Th <- c('Roc', 'Kappa', 'TSS')
  if(!(bin.method %in% Th)) stop("\n bin.method should be one of 'Roc', 'Kappa' or 'TSS'  \n")
  if(!(weight.method %in% Th)) stop("\n weight.method should be one of 'Roc', 'Kappa', or 'TSS' \n") 
  if(!(weight.method %in% names(Biomod.material$evaluation.choice))) stop("\n the weight.method selected was not run in Models \n")
  if(is.null(Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]])) stop("unknown Projection name \n")

  Comp <- c(FALSE, 'gzip', 'xz')
  if(!(compress %in% Comp)) stop("\n compress should be one of FALSE, 'gzip' or 'xz'  \n")
  if(compress == 'xz' & .Platform$OS.type == 'windows'){
    cat("\n compress arg was switch to 'gzip' for OS compatibility")
    compress  <- 'gzip'
  }

  #store the models wanted and shut down the models that cannot run (not trained)
  ens.choice <- Biomod.material[[paste("proj.", Proj.name, ".choice", sep="")]]
  for(j in Biomod.material$algo) if(!eval(parse(text=j))) ens.choice[j] <- FALSE
  
  #turn repetition models off if were not used for projections
  if(Biomod.material[[paste("proj.", Proj.name, ".repetition.models", sep="")]]==FALSE){
      repetition.models <- FALSE
      ProjRunEval <- 1     
      cat("repetition models cannot be used for consensus : they have not been used to render projections \n")
  } else ProjRunEval <- Biomod.material$NbRunEval + 1 # ProjRunEval is an object only used for constituting the stack to do consensus on considering if rep.models were run or not in Proj().
  
  if(!repetition.models & final.model.out){
    cat("only PA models are available for total consensus, they cannot be taken out \n ")          
    #incident of rep.models in Proj() -> data storage format which obliges to build a stack with reps even if not wanted  
    final.model.out <- FALSE                                                                                                         
  } 
  
  ### Test if we can work with raster stack or if we need to take raster layer one by one
#   if(Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]]){
    # load a stack created
    layerToLoad <- grep(Biomod.material$species.names[1] ,list.files(path=paste(getwd(),"/proj.", Proj.name, sep=""), pattern=".raster" ), value=TRUE)
    
    if( sum(grepl('_Bin', layerToLoad)) > 0 ){
      layerToLoad <- layerToLoad[- grep('_Bin', layerToLoad)]
    }
    
    if( sum(grepl('_Filt', layerToLoad)) > 0 ){
      layerToLoad <- layerToLoad[- grep('_Filt', layerToLoad)]
    }
    
    layerToLoad <- layerToLoad[1]

#     layerToLoad <- paste("Proj_", Proj.name, "_", Biomod.material$species.names[1], sep="")
#     
# #     if(Biomod.material$NbRepPA > 0){
# #       layerToLoad <- paste(layerToLoad,"_PA1", sep="")
# #     } else{
# #       layerToLoad <- paste(layerToLoad,"_full", sep="")
# #     }
# #     if(Biomod.material$NbRunEval > 0){
# #       layerToLoad <- paste(layerToLoad,"_rep1", sep="")
# #     }
#     
#     layerToLoad = paste(layerToLoad,"_", names(ens.choice)[which(ens.choice == TRUE)[1]], ".raster", sep="")
    
    load(paste(getwd(), "/proj.", Proj.name, "/", layerToLoad, sep=""))
  if(Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]]){ 
    work.with.stack = canProcessInMemory(get(layerToLoad), sum(ens.choice) + 2 )
    if(return.stack){
      return.stack = canProcessInMemory(get(layerToLoad), (sum(ens.choice) + 2) + (Biomod.material$NbSpecies*6) )
      if(!return.stack){
        cat("\n No stack will be returned because to avoid memory issue")
      }
    }
    
  } else {
    work.with.stack = canProcessInMemory(get(layerToLoad), Biomod.material$NbRun + 6 )
    if(return.stack){
      return.stack = canProcessInMemory(get(layerToLoad), Biomod.material$NbRun + 6 + Biomod.material$NbSpecies)
      if(!return.stack){
        cat("\n No stack will be returned because to avoid memory issue")
      }
    }
  }
    
#     proj.files <- grep(Biomod.material$species.names[i] ,list.files(path=paste(getwd(),"/proj.", Proj.name, sep=""), pattern=".raster" ), value=TRUE)[1]
#     eval(parse(text=paste("load('", proj.files, "')", sep="")))
#     
# #     load(paste(getwd(), "/proj.", Proj.name, "/Proj_", Proj.name, "_", Biomod.material$species.names[1], "_", names(ens.choice)[which(ens.choice == TRUE)[1]], ".raster", sep=""))
#     
# #     work.with.stack = canProcessInMemory(get(paste("Proj_", Proj.name, "_", Biomod.material$species.names[1], "_", names(ens.choice)[which(ens.choice == TRUE)[1]], ".raster", sep="")), sum(ens.choice) + 2 )
#     work.with.stack = canProcessInMemory(get(proj.files), sum(ens.choice) + 2 )

#   } else {
#     # load a layer created
#     layerToLoad <- paste("Proj_", Proj.name, "_", Biomod.material$species.names[1], sep="")
#     
#     if(Biomod.material$NbRepPA > 0){
#       layerToLoad <- paste(layerToLoad,"_PA1", sep="")
#     } else{
#       layerToLoad <- paste(layerToLoad,"_full", sep="")
#     }
#     if(Biomod.material$NbRunEval > 0){
#       layerToLoad <- paste(layerToLoad,"_rep1", sep="")
#     }
#     
#     layerToLoad = paste(layerToLoad,"_", names(ens.choice)[which(ens.choice == TRUE)[1]], ".raster", sep="")
#     
#     load(paste(getwd(), "/proj.", Proj.name, "/", layerToLoad, sep=""))
#     
#     work.with.stack = canProcessInMemory(get(layerToLoad), Biomod.material$NbRun + 6 )
#   }
  
  # determine the algo that will be computed
  ef.algo <- c('prob.mean','prob.mean.weighted')
  
  if(work.with.stack){
    ef.algo <- c( ef.algo, 'median')
  } else {cat("\n   ! Median of probabilities will not be computed !")}
  
  if(!is.na(Evaluation.results.Roc[1])){
    ef.algo <- c( ef.algo, 'Roc.mean')
  } else {cat("\n   ! Roc comitee averaging will not be computed !")}
  
  if(!is.na(Evaluation.results.Kappa[1])){
    ef.algo <- c( ef.algo, 'Kappa.mean')
  } else {cat("\n   ! Kappa comitee averaging will not be computed !")}
  
  if(!is.na(Evaluation.results.TSS[1])){
    ef.algo <- c( ef.algo, 'TSS.mean')
  } else {cat("\n   ! TSS comitee averaging will not be computed !")}
  
  # determine the decay
  if(is.numeric(decay)){
    prob.mean.weigth.decay <- decay
  } else if(decay == 'proportional'){
    prob.mean.weigth.decay <- 1
  } else{
    stop("decay arg must be a numeric or proportional")
  }
  
  return(list(compress = compress,
              repetition.models = repetition.models,
              ProjRunEval = ProjRunEval,
              final.model.out = final.model.out,
              ens.choice = ens.choice,
              work.with.stack = work.with.stack,
              ef.algo = ef.algo,
              prob.mean.weigth.decay = prob.mean.weigth.decay,
              return.stack = return.stack))
}


