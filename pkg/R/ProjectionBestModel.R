`ProjectionBestModel` <-
function(Proj.name, Bin.trans=TRUE, Filt.trans=TRUE, method='all')
{
    Th <- c('Kappa','TSS','Roc','all')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc' 'all'")
    if(is.null(Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]])) stop("unknown Projection name \n") 
    
    
    #apply the function to the 3 possible transformation methods if all are selected
    if(method == 'all') for(k in 1:3) ProjectionBestModel(Proj.name, Bin.trans, Filt.trans, method=Th[k])
    
    #run the function for one method at a time 
    else { if(Biomod.material$evaluation.choice[method]) {
    
        leng <- Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]]
        ARRAY.bin <- ARRAY.filt <- ARRAY <- array(NA, c(leng, max(Biomod.material$NbRun), Biomod.material$NbSpecies), dimnames=list(1:leng, rep(NA,max(Biomod.material$NbRun)), Biomod.material$species.names))

        #load PredictionBestModel resultsfor that method
        eval(parse(text=paste("load('pred/BestModelBy", method, "')", sep="")))

        i <- 1
        while(i <= Biomod.material$NbSpecies){     
            load(paste(getwd(), "/proj.", Proj.name, "/Proj_",Proj.name, "_", Biomod.material$species.names[i],sep=''))
            
            
            NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)   #considering the number of PA runs that were done   
            nbrep <- Biomod.material$NbRunEval + 1
    
            for(j in 1:NbPA){
                for(k in 1:nbrep){ 
            
                    #writing the name to use for getting the right info in Evaluation.results lists
                    if(Biomod.material$NbRepPA == 0) nam <- "full" else nam <- paste("PA", j, sep="")
                    if(k!=1) nam <- paste(nam, "_rep", k-1, sep="")
                    #assign name to column of matrix                    
                    dimnames(ARRAY)[[2]][(j-1)*nbrep+k] <- dimnames(ARRAY.bin)[[2]][(j-1)*nbrep+k] <- dimnames(ARRAY.filt)[[2]][(j-1)*nbrep+k] <- nam  
                    nam <- paste(Biomod.material$species.names[i], nam, sep="_")
            
            
                    best <- as.character(eval(parse(text=paste("BestModelBy",method,"[[Biomod.material$species.names[i]]][(j-1)*nbrep+k, 1]",sep=""))))   
                    
                    ARRAY[, (j-1)*nbrep+k, i] <- eval(parse(text=paste("Proj", Proj.name, Biomod.material$species.names[i], sep='_')))[, best,k,j]
                    if(Bin.trans)    eval(parse(text=paste("ARRAY.bin[, (j-1)*nbrep+k,i] <- BinaryTransformation(ARRAY[, (j-1)*nbrep+k,i], as.numeric(Evaluation.results.", method, "[[nam]][best,4]))",sep='')))
                    if(Filt.trans)   eval(parse(text=paste("ARRAY.filt[, (j-1)*nbrep+k,i] <- FilteringTransformation(ARRAY[, (j-1)*nbrep+k,i], as.numeric(Evaluation.results.", method, "[[nam]][best,4]))",sep='')))
                                
                }
            }
            i <- i + 1
        }
        
        
        #objects assignations and storage     
        assign(paste("Proj_",Proj.name,"_BestModelBy",method,sep=""), ARRAY)
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,", file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"')", sep="")))
            
        if(Bin.trans) {assign(paste("Proj_",Proj.name,"_BestModelBy",method,"_Bin",sep=""), ARRAY.bin)
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,"_Bin, file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Bin')", sep="")))}
        
        if(Filt.trans) {assign(paste("Proj_",Proj.name,"_BestModelBy",method,"_Filt",sep=""), ARRAY.filt)
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,"_Filt, file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Filt')", sep="")))}
        
    }}
}

