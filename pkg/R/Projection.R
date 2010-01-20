`Projection` <-
function(Proj=NULL, Proj.name, GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, Perc025=FALSE, Perc05=TRUE,
MDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE,
repetition.models=TRUE)
{
    require(nnet, quietly=T)
    require(rpart, quietly=T)
    require(Hmisc, quietly=T)
    require(Design, quietly=T)
    require(MASS, quietly=T)
    require(gbm, quietly=T)
    require(mda, quietly=T)
    require(randomForest, quietly=T)
    require(gam, quietly=T)
    
    dir.create(paste(getwd(), "/proj.", Proj.name, sep=""), showWarnings=F) #showWarnings=F -> permits overwritting of an already existing directory without signaling (dangerous?)
    
    if(BinRoc && !Biomod.material$evaluation.choice["Roc"] | FiltRoc && !Biomod.material$evaluation.choice["Roc"]) { BinRoc <- FiltRoc <- F ; cat("Roc cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinKappa && !Biomod.material$evaluation.choice["Kappa"] | FiltKappa && !Biomod.material$evaluation.choice["Kappa"]) { BinKappa <- FiltKappa <- F ; cat("Kappa cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinTSS && !Biomod.material$evaluation.choice["TSS"] | FiltTSS && !Biomod.material$evaluation.choice["TSS"]) { BinTSS <- FiltTSS <- F ; cat("TSS cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    
    
    #checking for the variable name compatibility with initial data
    nb <- 0
    for(i in 1:ncol(Proj)) if(sum(colnames(Proj)[i]==Biomod.material$VarNames) == 1) nb <- nb+1
    if(nb != Biomod.material$NbVar) stop("The variable names given do not correspond to the one used for calibrating the models \n Projections cannot proceed. \n") 
    
    #reorder the variables correctly 
    Proj <- Proj[match(Biomod.material$VarNames, colnames(Proj))]

    #check and error messages for the models that are wanted but not available
    algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    w <- names(which(!Biomod.material$algo.choice[names(which(algo.c))]))
    ww <- ""
    for(i in 1:length(w)) ww <- paste(ww, w[i])
    if(length(w) > 0) cat(paste("\n\n The following models can not be used to render projections : ", ww,"\n they have not been trained in Models() \n\n", sep=""))     
    algo.c[names(which(!Biomod.material$algo.choice))] <- F

    #save information on the projection
    Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]] <- nrow(Proj)
    Biomod.material[[paste("proj.", Proj.name, ".choice", sep="")]] <- algo.c
    Biomod.material[[paste("proj.", Proj.name, ".repetition.models", sep="")]] <- repetition.models
    assign("Biomod.material", Biomod.material, pos=1) 

    #the proj will be transformed for the models set to true in algo.cc. The point of this is for the GAM and GLM which will variably need or not transformations.
    algo.cc <- algo.c
    algo.cc['SRE'] <- F
    
     
    #------- projection loop per species -------   
    i <- 1
    while(i <= Biomod.material$NbSpecies){ 
        cat(paste(Biomod.material$species.names[i], " \n"))
        
        #------- defining the number of models to use for projecting -------
        NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)     #how many PA runs (even if NbRepPA was set to 0)
        if(repetition.models) Nbrep <- Biomod.material$NbRunEval +1 else Nbrep <- 1
        
        
        #------- create arrays to store projections -------
        PAs <- reps <- c()
        if(Biomod.material$NbRunEval != 0) for(j in 1:Biomod.material$NbRunEval) reps <- c(reps, paste("rep", j, sep="")) 
        if(Biomod.material$NbRepPA != 0) for(j in 1:NbPA) PAs <- c(PAs, paste("PA", j, sep="")) else PAs <- "no.PA"

        ARRAY <- array(NA, c(nrow(Proj), 9, Biomod.material$NbRunEval+1, NbPA), dimnames=list(1:nrow(Proj), Biomod.material$algo, c("total.data", reps), PAs))
        g  <- gg <- ggg <- gggg <- k <- kk <- kkk   <-   ARRAY    
        
       
        #------- looping for PAs, reps, and models -------    
        for(jj in 1:NbPA){
        
            if(Biomod.material$NbRepPA == 0) run.name <- "full"  else  run.name <- paste("PA", jj, sep="")
            
            for(Nrep in 1:Nbrep){ 
                for(a in Biomod.material$algo[algo.c]){ 
                    
                    run.name2 <- run.name
                    if(Nrep > 1) run.name2 <- paste(run.name2, "_rep", Nrep-1, sep="")
                    
                    if(exists("object")) rm(object)
                    
                    if(a != 'SRE'){
                        if(file.exists(paste(getwd(), "/models/", Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))){  # if model exists on hardisk
                        
                            if(!exists(paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))){    #no loading if already there in R 
                                object <- eval(parse(text=load(paste(getwd(), "/models/", Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))))
                            } else  object <- eval(parse(text=paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep="")))
                        
                        #remove the model loaded and kept in memory under its real name
                        ModelName <- paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep="")
                        rm(list=ModelName)
                        gc(reset=T)
                            
                        } else cat("WARNING: Could not find data for model", a, "evaluation repetition", Nrep, ". Probable cause : failure when running Models()", "\n")
                    } else object <- "SRE"
                    
                    if(exists("object")){
                        #------- making the projections with the model loaded -------# Special case of GLM and GAM which fill all arrays at once
                        
                        if(a == 'GLM') {
                            if(object$deviance == object$null.deviance) {  
                                algo.cc["GLM"] <- F    #in this case, the projections need no binary or filtered transformation 
                                if((sum(DataBIOMOD[,Biomod.material$NbVar+i])/nrow(DataBIOMOD)) < 0.5) g[,a,Nrep,jj] <- rep(0, nrow(Proj))
                                else  g[,a,Nrep,jj] <- gg[,a,Nrep,jj] <- ggg[,a,Nrep,jj] <- gggg[,a,Nrep,jj] <- k[,a,Nrep,jj] <- kk[,a,Nrep,jj] <- kkk[,a,Nrep,jj] <- rep(1000, nrow(Proj)) 
                            } else g[,a,Nrep,jj] <- as.integer(as.numeric(predict(object, Proj, type="response")) *1000)
                        }
            
                        if(a == 'GAM') {
                            if(object$deviance == object$null.deviance) {  
                                algo.cc["GAM"] <- F    #in this case, the projections need no binary or filtered transformation 
                                if((sum(DataBIOMOD[,Biomod.material$NbVar+i])/nrow(DataBIOMOD)) < 0.5) g[,a,Nrep,jj] <- rep(0, nrow(Proj))
                                else g[,a,Nrep,jj] <- gg[,a,Nrep,jj] <- ggg[,a,Nrep,jj] <- gggg[,a,Nrep,jj] <- k[,a,Nrep,jj] <- kk[,a,Nrep,jj] <- kkk[,a,Nrep,jj] <- rep(1000, nrow(Proj)) 
                            } else g[,a,Nrep,jj] <- as.integer(as.numeric(predict(object, Proj, type="response")) *1000)
                        }
                        
                        if(a == 'GBM') g[,a,Nrep,jj] <- as.integer(as.numeric(predict.gbm(object, Proj, Models.information[[i]]$GBM[[paste("PA", jj, sep="")]][[1]]$best.iter[[run.name2]], type='response')) *1000)
                        if(a == 'CTA') g[,a,Nrep,jj] <- as.integer(as.numeric(predict(object, Proj, type="vector")) *1000)
                        if(a == 'ANN') g[,a,Nrep,jj] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="raw")), type="range", OriMinMax=Models.information[[i]]$ANN[[paste("PA", jj, sep="")]][[1]]$RawPred[[run.name2]]) *1000) 
                        if(a == 'SRE') g[,a,Nrep,jj] <- as.integer(as.numeric(sre(DataBIOMOD[,Biomod.material$NbVar+i], DataBIOMOD[, 1:Biomod.material$NbVar], Proj, Perc025, Perc05)) *1000)
                        if(a == 'MDA') g[,a,Nrep,jj] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="post")[,2]), type="range", OriMinMax=Models.information[[i]]$MDA[[paste("PA", jj, sep="")]][[1]]$RawPred[[run.name2]]) *1000) 
                        if(a == 'MARS') g[,a,Nrep,jj] <- as.integer(Rescaler2(as.numeric(predict(object, Proj)), type="range", OriMinMax=Models.information[[i]]$MARS[[paste("PA", jj, sep="")]][[1]]$RawPred[[run.name2]]) *1000) 
                        if(a == 'RF') g[,a,Nrep,jj] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="prob")[,2]), type="range", OriMinMax=Models.information[[i]]$RF[[paste("PA", jj, sep="")]][[1]]$RawPred[[run.name2]]) *1000) 
                        
                        
                        #------- making the binary and filtered transformations if wanted -------#
                        if(algo.cc[a]){
                            if(BinRoc) gg[,a,Nrep,jj] <- as.numeric(BinaryTransformation(g[,a,Nrep,jj], as.numeric(Evaluation.results.Roc[[i]][a,4])))
                            if(FiltRoc) ggg[,a,Nrep,jj] <- as.numeric(FilteringTransformation(g[,a,Nrep,jj], as.numeric(Evaluation.results.Roc[[i]][a,4])))
                            if(BinKappa) gggg[,a,Nrep,jj] <- as.numeric(BinaryTransformation(g[,a,Nrep,jj], Evaluation.results.Kappa[[i]][a,4]))
                            if(FiltKappa) k[,a,Nrep,jj] <- as.numeric(FilteringTransformation(g[,a,Nrep,jj], Evaluation.results.Kappa[[i]][a,4]))
                            if(BinTSS) kk[,a,Nrep,jj] <- as.numeric(BinaryTransformation(g[,a,Nrep,jj], Evaluation.results.TSS[[i]][a,4]))
                            if(FiltTSS) kkk[,a,Nrep,jj] <- as.numeric(FilteringTransformation(g[,a,Nrep,jj], Evaluation.results.TSS[[i]][a,4]))
                        }
                         ggg[,'SRE',Nrep,jj] <- k[,'SRE',Nrep,jj] <- kkk[,'SRE',Nrep,jj] <-  g[,'SRE',Nrep,jj]   #filtered values for SRE = projections by SRE
                         gg[,'SRE',Nrep,jj] <- gggg[,'SRE',Nrep,jj] <- kk[,'SRE',Nrep,jj] <-  g[,'SRE',Nrep,jj]/1000   #binary values for SRE = projections by SRE /1000
                    }
                } #models       
            } #Nbrep -> coresponds to if repetition models were selected (==1 or ==NbRunEval+1)     
        } #NbPA 



        #------- exportation of the objects created in the working directory -------#           
        
        #list storing the names of the projections arrays produced to delete them afterwards
        ProjNameInList <- c()
        
        assign(paste("Proj",Proj.name,Biomod.material$species.names[i], sep="_"), g)
        eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],", file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"')", sep="")))
        #write.table(g, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],".txt", sep=""), row.names=F)
        ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i], sep=""))
        
        if(BinRoc){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"BinRoc", sep="_"), gg)
                   eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinRoc, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinRoc')", sep="")))
                   ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinRoc", sep=""))
        }  
                    
        if(FiltRoc){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"FiltRoc", sep="_"), ggg)
                    eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltRoc, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltRoc')", sep="")))
                    ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltRoc", sep=""))
        }   
                         
        if(BinKappa){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"BinKappa", sep="_"), gggg)
                     eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinKappa, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinKappa')", sep="")))
                     ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinKappa", sep=""))
        }
                     
        if(FiltKappa){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"FiltKappa", sep="_"), k)
                      eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltKappa, file='", getwd(),"/proj.", Proj.name, "/Proj_", Proj.name,"_",Biomod.material$species.names[i],"_FiltKappa')", sep="")))
                      ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltKappa", sep=""))
        }
                      
        if(BinTSS){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"BinTSS", sep="_"), kk)
                   eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinTSS, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinTSS')", sep="")))
                   ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_BinTSS", sep=""))
        }
    
        if(FiltTSS){assign(paste("Proj",Proj.name,Biomod.material$species.names[i],"FiltTSS", sep="_"), kkk)
                    eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltTSS, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltTSS')", sep="")))
                    ProjNameInList <- c(ProjNameInList, paste("Proj_",Proj.name,"_",Biomod.material$species.names[i],"_FiltTSS", sep=""))
        }
                             
        rm(list=c("g", "gg", "ggg", "gggg", "k" ,"kk", "kkk", "ARRAY", "object", ProjNameInList))
        gc(reset=T)
               
        i <- i+1
    }  #while
}

