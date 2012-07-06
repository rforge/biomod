# require(abind,quietly=TRUE)
# function(Proj=NULL, Proj.name=NULL, GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, quant=0.025,
# FDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE,
# repetition.models=TRUE, compress="xz")

setGeneric( "Projection_v2", 
            def = function(models.name,
                           modeling.work.dir = getwd(),
                           new.env.data ,
                           ...){
                            standardGeneric( "Projection_v2" )
                            } )

setMethod( 'Projection_v2', signature(new.env.data = 'data.frame'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           compress="xz",
           rescaled.models=TRUE){
        
    # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    .Models.dependencies(silent=TRUE)
  
    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')
    
    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    cat('\nDoing Models Projections...')
    if(length(grep('EF.',models.name)) > 0 ){
      kept.models.name <- models.name[-grep('EF.',models.name)] 
      kept.algo.run <- algo.run[-grep('EF.',algo.run)]
    } else {
      kept.models.name <- models.name
      kept.algo.run <- algo.run
    }
    
    proj.array <- lapply(kept.models.name, .Projection.do.proj, env=new.env.data, xy=xy)
    proj.array <- as.data.frame(proj.array)
    names(proj.array) <- kept.models.name
    
    # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(binary.proj)>0){
      cat("\nBinary transformations...")
      lapply(binary.proj, function(bin.proj){
        
        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
          }))

        proj.bin.array <- BinaryTransformation(proj.array, cuts)
        proj.bin.array <- DF_to_ARRAY(proj.bin.array)

        eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_array <- proj.bin.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_bin_",bin.proj,"_array' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_array , proj.bin.array, cuts)", sep="" )))
      })
    }

    
    # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(filtred.proj)>0){
      cat("\nFiltred transformations...")
      lapply(filtred.proj, function(filt.proj){
        
        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
        }))
        
        proj.filt.array <- FilteringTransformation(proj.array, cuts)
        proj.filt.array <- DF_to_ARRAY(proj.filt.array)
        
        eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_array <- proj.filt.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_filt_",filt.proj,"_array' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_array , proj.filt.array, cuts)", sep="" )))
      })
    }
    
    proj.array <- DF_to_ARRAY(proj.array)
    
    # 7. Saving projection on hard disk
    eval(parse(text = paste(proj.name,"_",sp.name, " <- proj.array", sep="")))
    eval(parse(text = paste("save(",proj.name,"_",sp.name, ", file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                            proj.name,"_",sp.name,"' )",sep="")))
    eval(parse(text = paste("rm(",proj.name,"_",sp.name,")", sep="" )))
    gc(reset=TRUE)
    
    return(invisible(proj.array))   
  })


setMethod( 'Projection_v2', signature(new.env.data = 'RasterStack'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,           
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           stack = TRUE,
           compress="xz",
           rescaled.models=TRUE){
        
    # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    .Models.dependencies(silent=TRUE)
  
    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')
       
    
    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    cat('\nDoing Models Projections...')
    if(length(grep('EF.',models.name)) > 0 ){
      kept.models.name <- models.name[-grep('EF.',models.name)] 
      kept.algo.run <- algo.run[-grep('EF.',algo.run)]
    } else {
      kept.models.name <- models.name
      kept.algo.run <- algo.run
    }
      
    proj.ras <- lapply(kept.models.name, .Projection.do.proj, env=new.env.data)

    # transform list of rasterLayers into a rasterStack
    proj.stack <- stack(x = proj.ras)

    layerNames(proj.stack) <- kept.models.name #names(proj.ras.mod)
    rm(proj.ras)
    
    # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(binary.proj)>0){
      cat("\nBinary transformations...")
      lapply(binary.proj, function(bin.proj){
        
        cuts <- unlist(lapply(layerNames(proj.stack), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
          }))

        proj.bin.stack <- BinaryTransformation(proj.stack, cuts)
        layerNames(proj.bin.stack) <- paste(layerNames(proj.stack), ".bin", sep="")

        eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_RasterStack <- proj.bin.stack", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack , proj.bin.stack, cuts)", sep="" )))
      })
    }
    
    # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(filtred.proj)>0){
      cat("\nFiltered transformations...")
      lapply(filtred.proj, function(filt.proj){
        
        cuts <- unlist(lapply(layerNames(proj.stack), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
        }))
        
        proj.filt.stack <- FilteringTransformation(proj.stack, cuts)
        layerNames(proj.filt.stack) <- paste(layerNames(proj.stack), ".filt", sep="")
        
        eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_RasterStack <- proj.filt.stack", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack , proj.filt.stack, cuts)", sep="" )))
      })
    }
    
    # 7. Saving projection on hard disk
     eval(parse(text = paste(proj.name,"_",sp.name, " <- proj.stack", sep="")))
     eval(parse(text = paste("save(",proj.name,"_",sp.name, ", file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                             proj.name,"_",sp.name,"' )",sep="")))
     eval(parse(text = paste("rm(",proj.name,"_",sp.name,")", sep="" )))
    gc(reset=TRUE)
    
    
    return(invisible(proj.stack))   
  })


setGeneric( ".Projection.do.proj", 
            def = function(model.name, env, model.dir = NULL,...){
                    standardGeneric( ".Projection.do.proj" )
                    } )

setMethod('.Projection.do.proj', signature(env='data.frame'),
  function(model.name, env, xy = NULL, model.dir = NULL, rescaled.models=TRUE){
    cat('\n***', model.name)
    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',unlist(strsplit(model.name, split="_"))[1],'/models', sep="")
    }
    
    # loading model
    if(length(c(grep('SRE',model.name), grep('MAXENT',model.name), grep('EF.', model.name))) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }
  
    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT') )){
      if(!grep('EF.',model.type))
        stop('Unknow model type')
    }
    
    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      # proj automaticly rescaled
      return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "raw")),
                                              name = model.name ) * 1000)))
    }
    
    if(model.type == 'CTA'){
      proj <- as.integer(as.numeric(predict(model.sp, env,type="prob")[,2]) * 1000)
    
      if(rescaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame( proj = proj) )
    }
    
    if(model.type == 'FDA'){
      # proj automaticly rescaled
      return( data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "posterior")[, 2]),
                                              name = model.name) * 1000)))
    }
    
    if(model.type == 'GBM'){
      best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      proj <- as.integer(predict.gbm(model.sp, env, best.iter, type = "response") * 1000)

      if(rescaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return(data.frame( proj = proj))
    }
    
    if(model.type == 'GLM'){
      proj <- as.integer(.testnull(model.sp, Prev, env) * 1000)

      if(rescaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame(proj = proj) )
    }
    
    if(model.type == 'GAM'){
      proj <- as.integer(.testnull(model.sp, Prev, env) * 1000)

      if(rescaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
 
      return( data.frame( proj = proj ) )
    }
    
    if(model.type == 'MARS'){
      # proj automaticly rescaled
      return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env)), 
                                              name = model.name) * 1000)))
    }
    
    if(model.type == 'RF'){
      proj <- as.integer(as.numeric(predict(model.sp,env, type='prob')[,'1']) *1000)
      
      if(rescaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame( proj = proj ))
    }
    
    if(model.type == 'SRE'){
      # loading data of the correspunding run
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))                
      return(eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
                                   model.name,"$Explanatory, env, Data_",model.name,
                                   "$Quant)*1000", sep=""))))
    }
    
    if(model.type == 'MAXENT'){
      if(!is.null(xy)){
        .Prepare.Maxent.Proj.WorkDir(env, xy)
        system(command=paste("java -cp maxent.jar density.Project ", model.dir,"/",
                             model.name,"/",sub("_MAXENT","",model.name),
                             ".lambdas MaxentTmpData/Proj_swd.csv MaxentTmpData/projMaxent", sep=""))
        
        maxent.proj <- read.asciigrid("MaxentTmpData/projMaxent.asc")@data
        
        .Delete.Maxent.WorkDir()
        return(proj = as.integer(.Rescaler5(as.numeric(maxent.proj), 
                                              name = model.name) * 1000))
      } else {
        cat('\n MAXENT need coordinates to run! NA returned ')
        return(data.frame(rep(NA,nrow(env))))
      }
    }
          
  })

setMethod('.Projection.do.proj', signature(env='RasterStack'),
  function(model.name, env, model.dir = NULL, rescaled.models=TRUE){
    cat('\n***', model.name)
    
    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',unlist(strsplit(model.name, split="_"))[1],'/models', sep="")
    }
    
    # loading model
    if(length(c(grep('SRE',model.name), grep('MAXENT',model.name), grep('EF.', model.name))) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }
  
    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT') )){
      if(!grep('EF.',model.type))
        stop('Unknow model type')
    }
    
    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model.sp, type="raw")
      proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(proj.ras * 1000))
    }
    
    if(model.type == 'CTA'){
      set.seed(123) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(rescaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'FDA'){
      pred.ras <- predict(env, model.sp, type="post", index=2)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000))
    }
    
    if(model.type == 'GBM'){
      if(file.exists(paste(model.dir,'/',model.name,'_best.iter'))){
        load(paste(model.dir,'/',model.name,'_best.iter'))
      } else{
        best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      }
        
      proj.ras <- predict(env, model.sp, n.trees=best.iter, type='response')
      if(rescaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'GLM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(rescaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'GAM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(rescaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'MARS'){
      pred.ras <- predict(env, model.sp)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000) )
    }
    
    if(model.type == 'RF'){
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(rescaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )      
    }
    
    if(model.type == 'SRE'){
#       cat('\n SRE prediction not supported yet ! ')
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))
      data.sre <- get(paste('Data_',model.name, sep=""))
      rm(list=paste('Data_',model.name, sep=""))
#       sre.out <- eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
#                                    model.name,"$Explanatory, env, Data_",model.name,
#                                    "$Quant)*1000", sep="")))
      sre.out <- raster::subset(sre(data.sre$Response, data.sre$Explanatory, env, data.sre$Quant), 1, drop=TRUE) * 1000
      
      return(sre.out)
    }
    
    if(model.type == 'MAXENT'){
        .Prepare.Maxent.Proj.Raster.WorkDir(env)
        system(command=paste("java -cp maxent.jar density.Project ", model.dir,"/",
                             model.name,"/",sub("_MAXENT","",model.name),
                             ".lambdas MaxentTmpData/Proj MaxentTmpData/projMaxent", sep=""))
        
        proj.ras <- raster("MaxentTmpData/projMaxent.asc")
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
        return(round(proj.ras*1000))
    }
    
       
  })

