`LoadProj` <-
function(Sp=1, Bin=FALSE, Filt=FALSE, which.pred="pred", raster=FALSE){

    #check the projection names wanted
    check.names <- c()
    
    for(i in 1:length(which.pred)){ 
        if(which.pred[i] == 'pred'){ check.names <- c(check.names, 'pred')
        } else {
            if(is.null(Biomod.material[[paste("proj.", which.pred[i], ".length", sep="")]])){ cat(paste("Warning \n Projection '", which.pred[i], "' unknown \n", sep=""))
            } else check.names <- c(check.names, which.pred[i])
         }
    }
    which.pred <- check.names
    
    
    #Check dimension calls for species
    if(length(Sp) > (Biomod.material$NbSpecies) | Sp=='all') Sp <- 1:(Biomod.material$NbSpecies)



    if(!raster){ 
    #------- Looping the loading of all the pred selected -------   
    for(i in 1:length(which.pred)){
         
            if(which.pred[i]=='pred'){ pname <- paste("/pred/Pred_", sep="")
            } else pname <- paste("/proj.", which.pred[i], "/Proj_", which.pred[i], "_", sep="") 
         
            for(k in Sp){
            
                #load probabilities
                pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], sep="") 
                if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
 
                #load bin and filt
                if(Bin){ for(j in 1:3){ if(Biomod.material$evaluation.choice[j]){
                    pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], "_Bin", names(Biomod.material$evaluation.choice)[j], sep="") 
                    if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
                }}}
                if(Filt){ for(j in 1:3){ if(Biomod.material$evaluation.choice[j]){
                    pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], "_Filt", names(Biomod.material$evaluation.choice)[j], sep="") 
                    if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
                }}}
                           
        }
    }
    } else{           #raster projections
    
     for(i in 1:length(which.pred)){
         
          pname <- paste("/proj.", which.pred[i], "/Proj_", which.pred[i], "_", sep="") 
          MoD <- Biomod.material[[paste("proj.", which.pred, ".choice", sep="")]]
       
          for(m in Biomod.material$algo[MoD]){
              for(k in Sp){
              
                  #load probabilities
                  pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], "_", m, ".raster", sep="") 
                  if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
   
                  #load bin and filt
                  if(Bin){ for(j in 1:3){ if(Biomod.material$evaluation.choice[j]){
                      pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], "_", m, "_Bin", names(Biomod.material$evaluation.choice)[j], ".raster", sep="") 
                      if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
                  }}}
                  if(Filt){ for(j in 1:3){ if(Biomod.material$evaluation.choice[j]){
                      pred.name <- paste(getwd(), pname, Biomod.material$species.names[k], "_", m, "_Filt", names(Biomod.material$evaluation.choice)[j], ".raster", sep="") 
                      if(file.exists(pred.name)) load(pred.name, envir=.GlobalEnv)
                  }}}
                             
              }
          }
     }}
}