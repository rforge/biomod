`Projection.raster` <-
function(RasterProj=NULL, Proj.name, GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, quant=0.025,
FDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE,
repetition.models=TRUE, stack.out=TRUE)
{
    require(nnet, quietly=TRUE)
    require(rpart, quietly=TRUE)
    require(Hmisc, quietly=TRUE)
    require(Design, quietly=TRUE)
    require(MASS, quietly=TRUE)
    require(gbm, quietly=TRUE)
    require(mda, quietly=TRUE)
    require(randomForest, quietly=TRUE)
    require(gam, quietly=TRUE)
    
    require(foreign, quietly=TRUE)
    require(sp, quietly=TRUE)
    require(rgdal, quietly=TRUE)
    require(raster, quietly=TRUE)
    require(maptools, quietly=TRUE)
    
    
    
    #Check wether data input is a raster stack and corresponds to the variables used for calibration
    if(class(RasterProj) != "RasterStack") stop("\n Entry in 'RasterProj' should be an object of the class 'RasterStack' ")
    #checking for the variable name compatibility with initial data
    nb <- 0
    for(i in layerNames(RasterProj)) if(sum(i==Biomod.material$VarNames) == 1) nb <- nb+1
    if(nb != Biomod.material$NbVar) stop("The variable names in RasterProj do not correspond to the one used for calibrating the models \n Projections cannot proceed. \n") 

    #reorder the variables correctly 
    #RasterProj <- RasterProj[match(Biomod.material$VarNames, colnames(RasterProj))]
    
    if(BinRoc && !Biomod.material$evaluation.choice["Roc"] | FiltRoc && !Biomod.material$evaluation.choice["Roc"]) { BinRoc <- FiltRoc <- FALSE ; cat("Roc cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinKappa && !Biomod.material$evaluation.choice["Kappa"] | FiltKappa && !Biomod.material$evaluation.choice["Kappa"]) { BinKappa <- FiltKappa <- FALSE ; cat("Kappa cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinTSS && !Biomod.material$evaluation.choice["TSS"] | FiltTSS && !Biomod.material$evaluation.choice["TSS"]) { BinTSS <- FiltTSS <- FALSE ; cat("TSS cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
 
    #check and error messages for the models that are wanted but not available
    algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, FDA=FDA, RF=RF, SRE=SRE)
    w <- names(which(!Biomod.material$algo.choice[names(which(algo.c))]))
    ww <- ""
    for(i in 1:length(w)) ww <- paste(ww, w[i])
    if(length(w) > 0) cat(paste("\n\n The following models can not be used to render projections : ", ww,"\n they have not been trained in Models() \n\n", sep=""))     
    algo.c[names(which(!Biomod.material$algo.choice))] <- FALSE
    
  
   
   
    dir.create(paste(getwd(), "/proj.", Proj.name, sep=""), showWarnings=FALSE) #showWarnings=FALSE -> permits overwritting of an already existing directory without signaling (dangerous?)
    
    #save information on the projection 
    Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]] <- dim(RasterProj)[1] * dim(RasterProj)[2]
    Biomod.material[[paste("proj.", Proj.name, ".choice", sep="")]] <- algo.c
    Biomod.material[[paste("proj.", Proj.name, ".repetition.models", sep="")]] <- repetition.models
    Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]] <- stack.out
    assign("Biomod.material", Biomod.material, pos=1) 

    
    
    #Create blank raster
    BLANK.ras <- RasterProj@layers[[1]]
    BLANK.ras[!is.na(BLANK.ras)] <- 0
     
      
    #------- projection loop per species -------   
    i <- 1
    while(i <= Biomod.material$NbSpecies){ 
        cat(paste(Biomod.material$species.names[i], " \n"))
        
        #------- defining the number of models to use for projecting -------
        NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)     #how many PA runs (even if NbRepPA was set to 0)
        if(repetition.models) Nbrep <- Biomod.material$NbRunEval +1 else Nbrep <- 1
        
        #Defining the final run
        if(repetition.models) FinalRun <- Biomod.material[["NbRun"]][i] else FinalRun <- Biomod.material[["NbRun"]][i] / (Biomod.material$NbRunEval+1)

            
        #------- looping for PAs, reps, and models -------    
        for(a in Biomod.material$algo[algo.c]){  
            #for storing the stacks if needed
            pile.proj <- stack()
            pile.names <- pile.thT <- pile.thK <- pile.thR <- c()
        
            for(jj in 1:NbPA){
                if(Biomod.material$NbRepPA == 0) run.name <- "full"  else  run.name <- paste("PA", jj, sep="")

                for(Nrep in 1:Nbrep){
                    run.name2 <- run.name
                    if(Nrep > 1) run.name2 <- paste(run.name2, "_rep", Nrep-1, sep="")
                    
                    
                    
                    #Loading the model for this run
                    if(exists("object")) rm(object) 
                    if(a != 'SRE'){
                        if(file.exists(paste(getwd(), "/models/", Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))){  # if model exists on hardisk
                        
                            if(!exists(paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))){   #no loading if already there in R 
                                object <- eval(parse(text=load(paste(getwd(), "/models/", Biomod.material$species.names[i], "_", a, "_", run.name2, sep=""))))
                            } else  object <- eval(parse(text=paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep="")))
                            
                        } else cat("WARNING: Could not find data for model", a, "evaluation repetition", Nrep, ". Probable cause : failure when running Models()", "\n")
                    } else object <- "SRE"
                    
                    
                    
                    #------- making the projections with the model loaded -------# 
                    if(exists("object")){   
                        
                        if(a == 'GLM')  g <- predict(RasterProj, model=object, type='response')
                        if(a == 'GAM')  g <- predict(RasterProj, model=object, type='response')
                        if(a == 'GBM')  g <- predict(RasterProj, object, n.trees=GBM.perf[[i]][[run.name2]], type='response') 
                        if(a == 'ANN')  {
                        	set.seed(555)
                        	g <- predict(RasterProj, object, type="raw")
                        }	
                        if(a == 'FDA')  g <- predict(RasterProj, object, type="post", index=2)
                        if(a == 'MARS') g <- predict(RasterProj, object)  
                        if(a == 'RF')   g <- predict(RasterProj, model=object, type='prob', index=2)
                        if(a == 'CTA')  g <- predict(RasterProj, model=object, type='prob', index=2)
                        if(a == 'SRE')  g <- sre(DataBIOMOD[,Biomod.material$NbVar+i], DataBIOMOD[, 1:Biomod.material$NbVar], RasterProj, quant)
                                             
                        #g <- as.numeric(g)  
                        #Rescale prediction for the models that need to
                        if(any(c("ANN", "FDA", "MARS")==a)) g <- .Rescaler4(g, run=paste(Biomod.material$species.names[i], "_", a, "_", run.name2, sep="")) 
                        #Transform into integer values
                        g <- round(g*1000) 
                                                                                                                                           
                     } else g <- BLANK.ras                                                                                                       #assign blank raster -> no model = no prediction
                        
                    
                    #------ making the binary and filtered transformations if wanted ------#
                    #------ exportation of the objects created in the working directory -----#                      
                    evals <- rep(c('Roc', 'Kappa', 'TSS'), 2)
                    trans <- c('BinRoc','BinKappa','BinTSS','FiltRoc','FiltKappa','FiltTSS')
                    
                    if(stack.out!=TRUE){ 
                    
                        eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material$species.names[i],"_", run.name2, "_", a, ".raster <- g", sep="")))
                        eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"_", run.name2, "_", a, ".raster, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_", run.name2, "_", a, ".raster', compress='xz')", sep="")))
 
                        if(a != 'SRE'){
                            for(ET in 1:6){ if(eval(parse(text=trans[ET]))){
                                Thresh <- as.numeric(eval(parse(text=paste("Evaluation.results.", evals[ET], "[[i]][a,4]", sep=""))))   #get the threshold for that proj (defined by 'trans')
                              #  gg <- g                                                                                                  #order between 'evals' and 'trans' is important
                              
						        gg <- (g >= Thresh)
                                if(ET>3) gg <- (gg * g)
                               
                             
                               # gg[g<Thresh] <- 0
                                #if(ET<4) gg[g>=Thresh] <- 1                                                                              #transforming values over threshold to 1 if we want binary data only
                                        
                                nam <- paste("Proj_", Proj.name,"_", Biomod.material$species.names[i],"_", trans[ET], "_", run.name2, "_", a, ".raster", sep="")
                                #assign(nam, gg)                                                                                         #assign the data to the name wanted
                                eval(parse(text=paste(nam, "<-gg", sep="")))
                                eval(parse(text=paste("save(" ,nam, ", file='", getwd(),"/proj.", Proj.name, "/", nam, "', compress='xz')", sep="")))   #and save it on disk
                            }}   
                        } 
                    } else{
                    	pile.proj <- stack(pile.proj, run.name2=g)                                                                      #store the projection (g) for each value of the loop to its name (run.name2)
                        pile.names <- c(pile.names, run.name2)                                                                          #store the names of the projection to assign later to the stack
                        if(Biomod.material$evaluation.choice[[1]]) pile.thR <- c(pile.thR, as.numeric(Evaluation.results.Roc[[i]][a,4]))      #store the thresholds for each eval technic 
                        if(Biomod.material$evaluation.choice[[2]]) pile.thK <- c(pile.thK, as.numeric(Evaluation.results.Kappa[[i]][a,4]))
                        if(Biomod.material$evaluation.choice[[3]]) pile.thT <- c(pile.thT, as.numeric(Evaluation.results.TSS[[i]][a,4]))
                    
                        
                        if(jj*Nrep==FinalRun){                                                                                          #only run this if final run of the loop accroos reps -> all proj have been stacked in pile.proj
                            assign("pile.proj", pile.proj, pos=1)                                                                          
                            assign("pile.names", pile.names, pos=1) 
                            layerNames(pile.proj) <- pile.names                                                                         #assign names of proj to layers of the stack
                           # assign(paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", a, ".raster", sep=""), pile.proj)
                            eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", a, ".raster <- pile.proj", sep="")))
                            eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i], "_", a, ".raster, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_", a, ".raster', compress='xz')", sep="")))
                            rm(list=paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", a, ".raster", sep=""))        #delete the proj created with correct name just for storage
                            
                            
                            if(a != 'SRE'){ 
                                for(ET in 1:6){ if(eval(parse(text=trans[ET]))){  #proceed to transform in bin and filt
                                    
                                    gg <- stack()  ### create a stack to store projection for each repetition
                                    if(ET==1 | ET==4) Thresh <- pile.thR                                                                #set thresh to Roc, Kappa or TSS considering the run in loop
                                    if(ET==2 | ET==5) Thresh <- pile.thK
                                    if(ET==3 | ET==6) Thresh <- pile.thT
                                    
                                    for(NB in 1:(jj*Nrep)){    
                                       temp <- (pile.proj@layers[[NB]] >= Thresh[NB]) #transforming values over threshold to 1, and 0 else
                                        if(ET>3) temp <- (pile.proj@layers[[NB]] * temp) # for filtering => actual values above the threshold. 
                                    	gg <- addLayer(gg, temp)
                                    }
                                    layerNames(gg) <- pile.names        
                                    nam <- paste("Proj_", Proj.name,"_", Biomod.material$species.names[i],"_", trans[ET], "_", a, ".raster", sep="")
                                   eval(parse(text=paste(nam, "<-  gg", sep="")))
                                    eval(parse(text=paste("save(" ,nam, ", file='", getwd(),"/proj.", Proj.name, "/", nam, "', compress='xz')", sep=""))) 
                                    rm(gg, list=nam)
                                }}   
                            } 
                            rm(pile.proj, g)
                        } #last run per model                        
                    } #stack.out
                    
                    
                    
                    
                    
                    #SRE
                    #else{
                    #    assign(paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_Bin.raster", sep=""), g/1000)
                    #    eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material$species.names[i],"Bin.raster, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material$species.names[i],"_Bin.raster')", sep="")))                            
                    #}                             
             
             
             
                 
                } #Nbrep -> coresponds to if repetition models were selected (==1 or ==NbRunEval+1)        
            } #NbPA  
        } #models 'a' loop                      
        
        i <- i+1
    }  #while species 'i' loop
    
    #save the history and workspace
    if(Biomod.material[["NbSpecies"]]==1) filename <- paste(Biomod.material[["species.names"]], "_run", sep="") else filename <- 'Biomod_run' 
    save.image(paste(filename, ".RData", sep=""))

}
                       