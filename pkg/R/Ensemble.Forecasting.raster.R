`Ensemble.Forecasting.raster` <-
function(ANN=TRUE,CTA=TRUE,GAM=TRUE,GBM=TRUE,GLM=TRUE,MARS=TRUE,FDA=TRUE,RF=TRUE,SRE=TRUE, Proj.name, weight.method, decay=1.6, PCA.median=TRUE, binary=TRUE, bin.method='Roc', Test=FALSE, repetition.models=TRUE, final.model.out=FALSE, qual.th=0)
{
    #Verify if user input is correct
    Th <- c('Roc', 'Kappa', 'TSS')
    if(!any(Th == bin.method)) stop("\n bin.method should be one of 'Roc', 'Kappa' or 'TSS'  \n")
    if(!any(Th == weight.method)) stop("\n weight.method should be one of 'Roc', 'Kappa', or 'TSS' \n") 
    if(sum(weight.method==names(Biomod.material$evaluation.choice))!=1) stop("\n the weight.method selected was not run in Models \n")
    if(is.null(Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]])) stop("unknown Projection name \n")

    #store the models wanted and shut down the models that cannot run (not trained)
    ens.choice <- p.choice <- Biomod.material[[paste("proj.", Proj.name, ".choice", sep="")]]
    for(j in Biomod.material$algo) if(!eval(parse(text=j))) ens.choice[j] <- F
    
    #create storing list of outputs
    list.out <- vector('list', Biomod.material$NbSpecies)
    names(list.out) <- Biomod.material$species.names
    
    #turn repetition models off if were not used for projections
    if(Biomod.material[[paste("proj.", Proj.name, ".repetition.models", sep="")]]==F){
        repetition.models <- F
        ProjRunEval <- 1     
        cat("repetition models cannot be used for consensus : they have not been used to render projections \n")
    } else ProjRunEval <- Biomod.material$NbRunEval + 1                                                                             #ProjRunEval is an object only used for constituting the stack to do
                                                                                                                                     #consensus on considering if rep.models were run or not in Proj().
    if(!repetition.models & final.model.out){                                                                                        #incident of rep.models in Proj() -> data storage format which 
        cat("only PA models are available for total consensus, they cannot be taken out \n ")                                        #obliges to build a stack with reps even if not wanted  
        final.model.out <- F                                                                                                         
    }                                                                                                                                

                                                                                                                                     
    
    #--------- start species loop ---------       
    for(i in 1:Biomod.material$NbSpecies){
        cat(paste(Biomod.material$species.names[i], " \n"))
        
        if(Test) nbrun <- 2 else nbrun <- 1
        nbrun <- 1        # no test can work yet
        for(ii in 1:nbrun){

            #considering the number of PA runs and reps that were actually done 
            NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)                                
            nbrep <- 1
            if(repetition.models) nbrep <- nbrep + Biomod.material$NbRunEval
            #store the thresholds produced by the ensemble forecasts
            ths <- vector('list', 6)
            #create a list "out" that will be used for storing the information on weights awarded, evaluation results, pca.median model selected, etc.
            out <- list()
            out[["weights"]]    <- matrix(NA, nr=nbrep*NbPA, nc=9, dimnames=list(rep("rep", nbrep*NbPA), Biomod.material$algo))
            out[["PCA.median"]] <- matrix(NA, nr=nbrep*NbPA, nc=1, dimnames=list(rep("rep", nbrep*NbPA), "model.selected"))
            out[["thresholds"]] <- matrix(NA, nr=6, nc=nbrep*NbPA, dimnames=list(c('prob.mean','prob.mean.weighted','median','Roc.mean','Kappa.mean','TSS.mean'), rep("rep", nbrep*NbPA)))
          
            #define which models correspond to each run -> will be used to drop layers from full stack of raster projections
            if(Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]]==T){ RunXModels <- rep(F, ProjRunEval*NbPA)                                     #difference if stack or not -> later : STK not containing same info, hence not same indices to select from
            } else RunXModels <- rep(F, nbrep*NbPA)    
            RunXModels[1] <- T
            RunXModels <- rep(RunXModels, sum(ens.choice))                             
            
            #goal : 'gnames' containing the names for loading the projections for the case when stack.out=F when running projection.raster()
            #and assigning the layer names for the storing stacks per species 
            if(Biomod.material$NbRepPA == 0) { PAs <- rep("full", nbrep)
            } else{ 
                PAs <- c()
                for(j in 1:NbPA) PAs <- c(PAs, paste("PA", j, sep=""))
                PAs <- rep(PAs, each=nbrep)
            }
            reps <- "" 
            if(nbrep != 1) for(j in 1:(nbrep-1)) reps <- c(reps, paste("rep", j, sep=""))
            if(Biomod.material$NbRepPA != 0) reps <- rep(reps, NbPA)
            gnames <- rep(NA, length(reps))
            for(j in 1:length(reps)){ if(reps[j]!="") gnames[j] <- paste(PAs[j], reps[j], sep="_") else gnames[j] <- PAs[j] }

            #Storing rasters for the 6 methods ( mean, median, etc.)
            STACK.m <- STACK.med <- STACK.w <- STACK.R <- STACK.K <- STACK.T <- stack()


         
            #load the projection data (generated by the Projection.raster() function)
            if(ii==1){
                STK <- stack()                                                                                                                                #Load Rasters
                for(m in Biomod.material$algo[ens.choice]){
            
                    if(Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]]==T){                                                                    #if rasters were saved in stacks
                        load(paste(getwd(), "/proj.", Proj.name, "/Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", m, ".raster", sep=""))
                        STK <- stack(STK, eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", m, ".raster", sep=""))))
                    } else{
                        for(Rep in 1:(nbrep*NbPA)){                                                                                                           #load rasters saved individualy
                            load(paste(getwd(), "/proj.", Proj.name, "/Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", gnames[Rep], "_", m, ".raster", sep=""))
                            STK <- stack(STK, eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", gnames[Rep], "_", m, ".raster", sep=""))))
                            #remove the raster piled into the stack
                            eval(parse(text = paste("rm(Proj_", Proj.name, "_", Biomod.material$species.names[i], "_", gnames[Rep], "_", m, ".raster)", sep = "")))
                        }
                    }
                }
            } else {                                                                                                                                          #Load Test dataset
                load(paste(getwd(), "/pred/Pred_", Biomod.material$species.names[i], sep=""))
                sp.data <- eval(parse(text=paste("Pred_", Biomod.material$species.names[i], sep="")))
            }              
            

           
            
            #--------- PAs and reps loops ---------  
            for(j in 1:NbPA){ #Loop through PAs replicates
                for(k in 1:nbrep){ #Loop through Eval replicates (nbrep set to 1 if no reps wanted -> 1 is PA 100% calib model)
                 
                 
                    #writing the name to use for getting the right info in Evaluation.results lists
                    if(Biomod.material$NbRepPA==0) nam <- "full" else nam <- paste("PA", j, sep="")
                    if(k!=1) nam <- paste(nam, "_rep", k-1, sep="")
                    
                    #Check if there was model fails in Models() -> turn this model off for that run
                    RUNens.choice <- ens.choice
                    for(a in Biomod.material$algo) if(ens.choice[a]) 
                        if(!file.exists(paste(getwd(), "/models/", Biomod.material$species.names[i], "_", a, "_", nam, sep=""))) RUNens.choice[a] <- F

                    #adding species name after check of model fails for convenience
                    nam <- paste(Biomod.material$species.names[i], nam, sep="_")

                    #set models to false if under the quality threshold
                    if(Biomod.material$NbRunEval!=0) whichEval <- 1 else whichEval <- 3
                    for(a in Biomod.material$algo) if(RUNens.choice[a]) if(as.numeric(eval(parse(text=paste("Evaluation.results.", weight.method, sep='')))[[nam]][a,whichEval]) < qual.th) RUNens.choice[a] <- F  #Weights are based on the Cross-validated evaluation values.                       

                    
                    
                    
                    #If all models are set to false, skip to next rep
                    if(sum(RUNens.choice)!=0){  

                        #drop layers of the original full stack and constitute one with just this runs' (PA/rep) projections
                        if(Biomod.material[[paste("proj.", Proj.name, ".stack", sep="")]]==T){
                            if(!repetition.models){ LayersToDrop <- c(1:length(RunXModels))[-c(which(RunXModels) + ((j-1)*(ProjRunEval)))]                              #layers to drop to create the stack to use for consensus for this run
                            } else LayersToDrop <- c(1:length(RunXModels))[-c(which(RunXModels)+((j-1)*(ProjRunEval)+k-1))]                                             #what is inside the -c() is in fact the layers we want to get             
                        } else                                                                                            #no stack.out in Proj() -> no extra layers than wanted in STACK
                            LayersToDrop <- c(1:length(RunXModels))[-c(which(RunXModels)+((j-1)*(nbrep)+k-1))]                     
                        
                        STKcons <- dropLayer(STK, LayersToDrop)                                                           #STKcons = data containing all models to do consensus for that run (PA, run)
                        
                        #then restrict data to models set to True
                        STKcons <- dropLayer(STKcons, as.numeric(which((ens.choice!=RUNens.choice)[ens.choice])))
    

    
      
    
                        if(sum(RUNens.choice)>1){                                                                            #if more than 1 model is wanted for ensemble forecating
                                                                       
                            #Mean and Median ensemble forecasting
                            #Take out models that failed (still in STKcons)
                            STACK.m <- stack(STACK.m, round(mean(STKcons)))
                            STACK.med <- stack(STACK.med, round(mean(STKcons)))                                           #median not functioning
    

                            #-------- binary results means ensemble forecating ---------#                                 #mean of the binary projections accross all selected techniques
                            for(jj in 1:3){ if(Biomod.material$evaluation.choice[Th[jj]]){
                            
                                #create a raster to accumulate the binary prediction for each model successively 
                                BLANK.ras <- STK@layers[[1]]
                                BLANK.ras[!is.na(BLANK.ras)] <- 0
                                ras.bin <- BLANK.ras
                                
                                for(kk in Biomod.material$algo[RUNens.choice]){ if(kk!='SRE'){
                                    BIN <- STKcons@layers[[which(Biomod.material$algo[RUNens.choice]==kk)]] >= as.numeric(eval(parse(text=paste("Evaluation.results.", Th[jj], sep="")))[[nam]][kk,4])
                                    BIN["TRUE"] <- 1
                                    ras.bin <- ras.bin + BIN
                                }}
                                if(RUNens.choice['SRE']) ras.bin <- ras.bin + STKcons@layers[[sum(RUNens.choice)]]/1000            #SRE already in binary (sum(ens) -> gives position of SRE (because last))
                                
                                if(Th[jj]=='Roc') STACK.R <- stack(STACK.R, round(ras.bin/sum(RUNens.choice)*1000))  
                                if(Th[jj]=='Kappa') STACK.K <- stack(STACK.K, round(ras.bin/sum(RUNens.choice)*1000)) 
                                if(Th[jj]=='TSS') STACK.T <- stack(STACK.T, round(ras.bin/sum(RUNens.choice)*1000))  
                            }}
    

                            #---------- Weighted Average Ensemble Forecasting --------# 
                            #This is like a mean accross all selected methods but with a weight associated to each technique depending on its score during evaluation
                            #Recover the weights (depending on the chosen "weight.method" from the "Evaluation.results.weightMethod" object)
                            wk <- p.choice
                            if(Biomod.material$NbRunEval!=0) whichEval <- 1 else whichEval <- 3                   # Weights are based on the Cross-validated evaluation values.
                            for(a in Biomod.material$algo) if(RUNens.choice[a]) wk[a] <- as.numeric(eval(parse(text=paste("Evaluation.results.", weight.method, sep='')))[[nam]][a,whichEval]) else wk[a] <- NA  
                            if(weight.method=='Roc') wk['SRE'] <- 0
                            #deal with cases where there are no scores >0, or rep model failed
                            if(sum(wk!=0, na.rm=T)==0) wk[wk==0] <- 0.1                                           # 0.1 = arbitrary value >0  -> those models will be used and not set to NA by next line
                            wk[wk==0] <- NA
    
                            # Calculate and attribute Weights to each modelling techniques
                            if(decay=="proportional"){                                                            # proportional: the weights are proportional to the chosen evaluation value (also relatively to the other scores)
                                wk[is.na(wk)] <- 0
                                if(weight.method=='Roc') wk[wk!=0] <- (wk[wk!=0]-0.5)*2
                                W <- wk/sum(wk)
                            }
                            if(is.numeric(decay)){                                                                # weights are "decay" times decreased for each subsequent model in model quality order.                            
                                if(sum(is.na(wk))<8){
                                    z <-rep(1,sum(!is.na(wk)))
                                    for(wj in 2:sum(!is.na(wk))) z[wj] <- z[wj-1]*decay
                                    z <- c(rep(0,(length(wk)-sum(!is.na(wk)))), z/sum(z))
                                    
                                    #determine which weight for which model
                                    wk[is.na(wk)] <- 0
                                    W <- rep(0,9)
                                    for(m in 1:9) {
                                      	if(sum(wk[m]==wk)!=1){ #if 2 or more score are identical -> make a mean weight between the ones concerned
                                      		  if(!is.na(wk[m])){
                                      		      for(nbm in 1:sum(wk[m]==wk)) W[m] <- W[m] + z[sum(wk[m]>wk)+nbm]
                                      		    	W[m] <- W[m] / sum(wk[m]==wk) 
                                      		  }
                                      	} else W[m] <- z[sum(wk[m]>wk)+1]
                                    }
                                } else if(sum(is.na(wk))==8) { wk <- is.na(wk) ; wk[T] <- 1 }
                            }                        
                            

                            #applying weights to projections
                            ras.mean.w <- BLANK.ras 
                            for(m in 1:9) if(RUNens.choice[m]) ras.mean.w <- ras.mean.w + STKcons@layers[[which(Biomod.material$algo[RUNens.choice]==Biomod.material$algo[m])]] * W[m]
                            STACK.w <- stack(STACK.w, round(ras.mean.w))                  
                            
                            #calculating the weighted threshold to convert the weighted probabilities to binary and/or filtered values     
                            thmi <- thpondi <- c()
                            for(a in Biomod.material$algo[RUNens.choice]) {
                                thmi    <- c(thmi,    eval(parse(text=paste("Evaluation.results.", bin.method, sep="")))[[nam]][a,4])
                                thpondi <- c(thpondi, eval(parse(text=paste("Evaluation.results.", weight.method, sep="")))[[nam]][a,4])
                            }
                            ths[[1]] <- c(ths[[1]], mean(as.numeric(thmi), na.rm=T))
                            thpondi[is.na(thpondi)] <- 0
                            ths[[2]] <- c(ths[[2]], sum(as.numeric(thpondi)*W[RUNens.choice]))   
                            #-----------------------end weights------------------------#
                                          
                                
                            #determine the model selected by the PCA consensus approach
                            #DISABLED  
                        
                            #store the information for each run
                            out[["thresholds"]][,(j-1)*nbrep+k] <-  c(mean(as.numeric(thmi), na.rm=T), sum(as.numeric(thpondi)*W[RUNens.choice]) ,NA,500,500,500)
                            out[["weights"]][(j-1)*nbrep+k, ] <- round(W,digits=4)              
                          
                          

                        
                        } else {                                                                                #only one model is available/wanted  
                                                                                                                                                                                  
                            #STKcons contains only the proj wanted
                            STACK.m <- stack(STACK.m, STKcons)
                            STACK.w <- stack(STACK.w, STKcons)
                            STACK.med <- stack(STACK.med, STKcons)
                          
                            if(Biomod.material$algo[RUNens.choice] != 'SRE'){  
                                #binary values
                                for(jj in 1:3){ if(Biomod.material$evaluation.choice[Th[jj]]){
                                    thresh <- as.numeric(eval(parse(text=paste("Evaluation.results.", Th[jj], sep="")))[[nam]][Biomod.material$algo[RUNens.choice],4])
                                    assign("thresh", thresh, pos=1)
                                    if(Th[jj]=='Roc'){   ST <- STKcons@layers[[1]] >= thresh ;  STACK.R <- stack(STACK.R, ST*1000) }   
                                    if(Th[jj]=='Kappa'){ ST <- STKcons@layers[[1]] >= thresh ;  STACK.K <- stack(STACK.K, ST*1000) } 
                                    if(Th[jj]=='TSS'){   ST <- STKcons@layers[[1]] >= thresh ;  STACK.T <- stack(STACK.T, ST*1000) }
                                
                                    #store thresholds
                                    if(Th[jj] == weight.method){
                                        ths[[1]] <- c(ths[[1]], thresh)
                                        ths[[2]] <- c(ths[[2]], thresh)
                                    }
                                }}
                            } else {
                                STACK.R <- stack(STACK.R, STKcons)
                                STACK.K <- stack(STACK.K, STKcons)
                                STACK.T <- stack(STACK.T, STKcons)                        
                            }   
                            out[["thresholds"]][,(j-1)*nbrep+k] <- c(ths[[1]][length(ths[[1]])], ths[[2]][length(ths[[2]])],NA,500,500,500)
                            out[["weights"]][(j-1)*nbrep+k, ] <- rep(0,9) ; out[["weights"]][(j-1)*nbrep+k, RUNens.choice] <- 1
                        }                                  
                    } else{ #if RUNens.choice!=0
                     ths[[1]] <- c(ths[[1]], NA)
                     ths[[2]] <- c(ths[[2]], NA)
                    }
                    
                    
                } #Evaluation replicates k loop
            } #PAs replicates j loop
            
                        
            ths[[4]] <- ths[[5]] <- ths[[6]] <- rep(500, nbrep*NbPA) 
            ths[[3]] <- ths[[1]]                                                                          #attribute the mean thresholds to median data (median th not relevant)
            rownames(out[["weights"]]) <- colnames(out[["thresholds"]]) <- gnames
            list.out[[i]] <- out 
            
            
            

            

            #----------------------If normal ensFor run or Test run-------------------------#
            if(ii==1){ #normal ensFor run --> saving on disk, transform in binary, total consensus
                
                #drop PA layers if not wanted for total consensus with final.model.out option
                #NB : take into consideration if all runs had at least one model for not selecting the wrong layers
                STACK.all <- stack(STACK.m, STACK.w, STACK.med)                                                             #the full stack takes all data even if PAs not wanted in total consensus
                thsC <- ths                                                                                                 #only used for line 315 : take PAs out if necessary
                if(final.model.out){
                    PaLayers <- rep(c(T, rep(F, nbrep-1)), NbPA) 
                    STACK.m <- dropLayer(STACK.m, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    STACK.w <- dropLayer(STACK.w, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    STACK.med <- dropLayer(STACK.med, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    for(j in 1:3) thsC[[j]] <- ths[[j]][PaLayers]
                }

                #save results on hard disk per species, stack all rasters, get the final consensus, transform to binary 
                #Final consensus : on all data available (average across all Evaluation Replicates and All PAs replicates (if wanted))   
                if(length(layerNames(STACK.m)) != 1) STACK.consTot <- stack(round(mean(STACK.m)), round(mean(STACK.w)), round(mean(STACK.med)))
                LNa <- c(paste("mean_", gnames, sep=""), paste("weigthed.mean_", gnames, sep=""), paste("median_", gnames, sep=""))                   #storing the layer names to assign to the final stacks
                LNaTot <- c("mean", "weighted.mean", "median")
                
                #convert the total consensus for the first 3 methods into binary
                if((nbrep*NbPA) != 1) STACK.consTot.bin <- stack()
                if((nbrep*NbPA) != 1) for(j in 1:3) STACK.consTot.bin <- stack(STACK.consTot.bin, STACK.consTot@layers[[j]] >= mean(thsC[[j]], na.rm=T))
                
                
                
                 
                #for each eval method available : stack the meanMethods, calculate the mean accross reps and convert it to binary    
                if(Biomod.material$evaluation.choice['Roc']){                
                    STACK.all <- stack(STACK.all, STACK.R)
                    if(final.model.out) STACK.R <- dropLayer(STACK.R, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    if((nbrep*NbPA) != 1) STACK.consTot <- stack(STACK.consTot, round(mean(STACK.R)))
                    if((nbrep*NbPA) != 1) STACK.consTot.bin <- stack(STACK.consTot.bin, mean(STACK.R) >= 500)
                    LNa <- c(LNa, paste("meanRoc_", gnames, sep="")) 
                    LNaTot <- c(LNaTot, "meanRoc")               
                }
                if(Biomod.material$evaluation.choice['Kappa']){
                    STACK.all <- stack(STACK.all, STACK.K)
                    if(final.model.out) STACK.K <- dropLayer(STACK.K, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    if((nbrep*NbPA) != 1) STACK.consTot <- stack(STACK.consTot, round(mean(STACK.K)))
                    if((nbrep*NbPA) != 1) STACK.consTot.bin <- stack(STACK.consTot.bin, mean(STACK.K) >= 500)
                    LNa <- c(LNa, paste("meanKappa_", gnames, sep=""))
                    LNaTot <- c(LNaTot, "meanKappa")
                }
                if(Biomod.material$evaluation.choice['TSS']){
                    STACK.all <- stack(STACK.all, STACK.T) 
                    if(final.model.out) STACK.T <- dropLayer(STACK.T, c(which(PaLayers[!is.na(out[["weights"]][,1])])))
                    if((nbrep*NbPA) != 1) STACK.consTot <- stack(STACK.consTot, round(mean(STACK.T))) 
                    if((nbrep*NbPA) != 1) STACK.consTot.bin <- stack(STACK.consTot.bin, mean(STACK.T) >= 500)        
                    LNa <- c(LNa, paste("meanTSS_", gnames, sep=""))
                    LNaTot <- c(LNaTot, "meanTSS")
                }                                                             
               
                assign("LNa1", LNa, pos=1)
                #Adapt layerNames to data : consider runs with no models selected
                LNa <- LNa[rep(!is.na(out[["weights"]][,1]), (length(LNa)/length(gnames)))]                 
                 
                                 
                                  
                #if binary=T, we transform the ensemble forecasting values into binary ones using the thresholds stored in ths (already done for final consensus)
                #first step (ths.lin) -> constitute a vector with all thresholds associated to the stack in the same order (for looping the binary conversion)
                if(binary){ 
                    ths.lin <- c()
                    for(j in 1:6){
                        if(j>3){ if(Biomod.material$evaluation.choice[Th[j-3]]) ths.lin <- c(ths.lin, ths[[j]][!is.na(ths[[1]])])                  
                        } else ths.lin <- c(ths.lin, ths[[j]][!is.na(ths[[1]])]) 
                    }
                    STACK.all.bin <- stack()
                    for(j in 1:length(layerNames(STACK.all))){                                                                                        #loop of binary conversion
                        J.bin <- STACK.all@layers[[j]] >= ths.lin[j]
                        STACK.all.bin <- stack(STACK.all.bin, J.bin)                        
                    }      
                }
               

                #if only one run was done there is no further calculation possible
                if((nbrep*NbPA) == 1){ STACK.consTot <- STACK.all ; STACK.consTot.bin <- STACK.all.bin ; layerNames(STACK.consTot) <- layerNames(STACK.consTot.bin) <- LNa }
                #save objects on disk (prior : assign objects to the name they should)
                layerNames(STACK.all) <- layerNames(STACK.all.bin) <- LNa
                layerNames(STACK.consTot) <- layerNames(STACK.consTot.bin) <- LNaTot 
                    
                #STACK.all
                assign(paste("consensus_", Biomod.material$species.names[i], "_", Proj.name, ".raster", sep=""), STACK.all) 
                eval(parse(text=paste("save(consensus_", Biomod.material$species.names[i], "_", Proj.name, ".raster, file='", getwd(),"/proj.", Proj.name, "/consensus_", Biomod.material$species.names[i], "_", Proj.name, ".raster')", sep="")))
                #STACK.all.bin
                assign(paste("consensus_", Biomod.material$species.names[i], "_", Proj.name, "_Bin.raster", sep=""), STACK.all.bin) 
                eval(parse(text=paste("save(consensus_", Biomod.material$species.names[i], "_", Proj.name, "_Bin.raster, file='", getwd(),"/proj.", Proj.name, "/consensus_", Biomod.material$species.names[i], "_", Proj.name,"_Bin.raster')", sep="")))
                #STACK.consTot
                assign(paste("Total_consensus_", Biomod.material$species.names[i], "_", Proj.name, ".raster", sep=""), STACK.consTot) 
                eval(parse(text=paste("save(Total_consensus_", Biomod.material$species.names[i], "_", Proj.name, ".raster, file='", getwd(),"/proj.", Proj.name, "/Total_consensus_", Biomod.material$species.names[i], "_", Proj.name,".raster')", sep="")))
                #STACK.consTot.bin
                assign(paste("Total_consensus_", Biomod.material$species.names[i], "_", Proj.name, "_Bin.raster", sep=""), STACK.consTot.bin) 
                eval(parse(text=paste("save(Total_consensus_", Biomod.material$species.names[i], "_", Proj.name, "_Bin.raster, file='", getwd(),"/proj.", Proj.name, "/Total_consensus_", Biomod.material$species.names[i], "_", Proj.name,"_Bin.raster')", sep="")))
                
            }#if ii==1 
            
            
            
            
            
                           
            
            
            #### Not working -> format problem (idea: call Ensemble.Foecasting ?) 
            
            if(ii==2){ #Test run : consensus methods done on current data (the test on the methods is done with AUC)
                test <- matrix(nc=nbrep*NbPA, nr=6, dimnames=list(c('prob.mean','prob.mean.weighted','median','Roc.mean','Kappa.mean','TSS.mean'),dimnames(ARRAY)[[2]]))

                for(j in 1:NbPA){
                    if(Biomod.material$NbRepPA != 0) lin <- Biomod.PA.sample[[Biomod.material$species.names[i]]][[j]]
                    else lin <- 1:dim(ARRAY.tot)[1]
                    
                    for(k in 1:nbrep){                                                                    # Loop through Model Evaluation replicates
                        for(m in 1:6){                                                                    # Loop through ensemble forecasting methods
                            if(m<4) {test[m,(j-1)*nbrep+k] <- somers2(ARRAY[,(j-1)*nbrep+k,m], DataBIOMOD[lin, Biomod.material$NbVar+i])["C"]   #to check if method was chosen
                            } else if(Biomod.material$evaluation.choice[Th[m-3]]) test[m,(j-1)*nbrep+k] <- somers2(ARRAY[,(j-1)*nbrep+k,m], DataBIOMOD[lin, Biomod.material$NbVar+i])["C"]    
                        }
                    }
                }
                list.out[[i]][["test.results"]] <- test       
            } 
            
          
            
             
               
        } #end of ii loop (test run (ensemble forecast calculated on calibration data) or normal ensFor run)       
    } #end of species i loop 
    
    #save list of info
    assign(paste("consensus_", Proj.name,"_results", sep=""), list.out, pos=1)
    eval(parse(text=paste("save(consensus_", Proj.name,"_results, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name,"_results')", sep="")))
    
    
    cat(paste("\n consensus_", Proj.name,"_results \n", sep=""))
    return(list.out)
}
