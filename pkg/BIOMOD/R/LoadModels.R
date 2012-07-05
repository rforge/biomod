`LoadModels` <-
function(Sp=1, PA='all', rep='all', models='all'){

    
    #check and info messages for the models 
    if(models!='all'){
        algo.c  <- c()  
        for(i in 1:9) if(sum(Biomod.material$algo[i]==models)>1) algo.c <- c(algo.c, Biomod.material$algo[i]==TRUE) else algo.c <- c(algo.c, Biomod.material$algo[i]==FALSE)  
        algo.c[names(which(!Biomod.material$algo.choice))] <- FALSE  #shut off the ones not available  
    
        w <- names(which(!Biomod.material$algo.choice[names(which(algo.c))]))
        ww <- ""
        for(i in 1:length(w)) ww <- paste(ww, w[i])
        if(length(w) > 0) cat(paste("\n\n The following models can not be used to render projections : ", ww,"\n they have not been trained in Models() \n\n", sep=""))     
        algo.c[names(which(!Biomod.material$algo.choice))] <- FALSE
        
    } else  algo.c <- Biomod.material$algo.choice 
    algo.c['SRE'] <- FALSE
    
    
    #Check dimension calls for species, rep, PA
    if(length(Sp) > (Biomod.material$NbSpecies) | Sp=='all') Sp <- 1:(Biomod.material$NbSpecies)
    if(length(rep) > (Biomod.material$NbRunEval+1) | rep=='all') rep <- 1:(Biomod.material$NbRunEval+1)
    if(Biomod.material$NbRepPA!=0){ if(length(PA) > (Biomod.material$NbRepPA) | PA=='all') PA <- 1:(Biomod.material$NbRepPA)
    } else PA <- 1

    
      
    #------- Looping the loading of all the models selected -------   
    for(i in Sp){
     
        for(j in PA){
            if(j==1){ if(Biomod.material$NbRepPA==0) pa<-"full" else pa<-"PA1"
            } else pa <- paste("PA", j, sep="")
        
            for(k in rep){
                if(k==1) re <- "" else re <- paste("_rep", k, sep="")
            
                for(kk in 1:9){ if(algo.c[kk]){   
     
     
                    model.name <- paste(Biomod.material$species.names[i], "_", names(algo.c[kk]), "_", pa, re, sep="") 
                    if(file.exists(paste(getwd(), "/models/", model.name, sep="")))
                    load(paste(getwd(), "/models/", model.name, sep=""), envir=.GlobalEnv)
     
                }}
            }   
        }
    }

}

