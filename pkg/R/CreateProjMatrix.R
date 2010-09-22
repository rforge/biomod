`CreateProjMatrix` <- function(ANN=TRUE, CTA=TRUE, GLM=TRUE, GBM=TRUE, GAM=TRUE, SRE=TRUE, FDA=TRUE, MARS=TRUE, RF=TRUE, Sp, Probs=TRUE, Bin=FALSE, Filt=FALSE, which.proj){

    #check for models that are wanted but not available
    algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, FDA=FDA, RF=RF, SRE=SRE)    
    algo.c[names(which(!Biomod.material$algo.choice))] <- F
    Format <- c('Probs','Bin','Filt')
    Th <- c('Roc','Kappa','TSS')


    for(spi in 1:length(Sp)){             
        for(PBF in Format){ if(eval(parse(text=PBF))){           #format loop
        
            Spname <- Sp[spi]
            DataFrame <- data.frame(NA)
            
            #load the first proj to know the dimensions
            projdata <- paste("Proj_", which.proj[1], "_", Spname, sep="")
            projj <- eval(parse(text=load(paste(getwd(), "/proj.", which.proj[1], "/", projdata, sep=""))))
            NbProjPerScen <- length(algo.c)*dim(projj)[4]*dim(projj)[3]
            
            #constitute the Groups matrix    
            Groups <- matrix(NA, 4, NbProjPerScen*length(which.proj))
            Groups[1,] <- c(rep(which.proj, each=NbProjPerScen))
            Groups[2,] <- c(rep(rep(dimnames(projj)[4][[1]], each=length(algo.c)*dim(projj)[3]), length(which.proj)))
            Groups[3,] <- c(rep(rep(dimnames(projj)[3][[1]], each=length(algo.c)), length(which.proj)*dim(projj)[4]))
            Groups[4,] <- c(rep(names(algo.c), dim(projj)[4]*dim(projj)[3]*length(which.proj)))
            
            if(PBF != "Probs"){
                Groups2 <- Groups
                for(Thi in 1:(sum(Biomod.material$evaluation.choice)-1)) Groups2 <- cbind(Groups2, Groups)
                Groups2 <- rbind(rep(names(Biomod.material$evaluation.choice[Biomod.material$evaluation.choice==T]), each=ncol(Groups)), Groups2)
            }
               
               
            if(PBF=='Probs'){
                
                for(ip in 1:length(which.proj)){  #projs
                     
                    #load the projection
                    projdata <- paste("Proj_", which.proj[ip], "_", Spname, sep="")
                    if(!exists(projdata)){                                                                            #not loading if already there in R 
                        projj <- eval(parse(text=load(paste(getwd(), "/proj.", which.proj[ip], "/", projdata, sep=""))))
                    } else  projj <- eval(parse(text=projdata))
                
                    for(i in 1:dim(projj)[4])        #PAs
                        for(j in 1:dim(projj)[3])        #reps        
                            DataFrame <- cbind(DataFrame, projj[, algo.c, j, i])                             
                }  #projs 
                assign(paste(Spname, "_ProjMat", sep=""), DataFrame[,-1], pos=1)
                assign(paste(Spname, "_ProjMat_Groups", sep=""), Groups, pos=1)    
            
            } else{         
            
                for(ip in 1:length(which.proj)){  #projs
                    
                    for(jj in 1:3){ if(Biomod.material$evaluation.choice[Th[jj]]){ 
                        #load the projection
                        projdata <- paste("Proj_", which.proj[ip], "_", Spname, "_", PBF, Th[jj], sep="")
                        if(!exists(projdata)){                                                                            #not loading if already there in R 
                            projj <- eval(parse(text=load(paste(getwd(), "/proj.", which.proj[ip], "/", projdata, sep=""))))
                        } else  projj <- eval(parse(text=projdata))
                    
                        for(i in 1:dim(projj)[4])        #PAs
                            for(j in 1:dim(projj)[3])        #reps        
                                DataFrame <- cbind(DataFrame, projj[, algo.c, j, i])                             
                
                    }}
                }  #projs               
                
                assign(paste(Spname, "_ProjMat_", PBF, sep=""), DataFrame[,-1], pos=1)
                assign(paste(Spname, "_ProjMat_", PBF, "_Groups", sep=""), Groups2, pos=1)    
            }               
              
        }} #format
    }#sp
}

