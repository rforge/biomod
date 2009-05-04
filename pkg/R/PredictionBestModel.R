`PredictionBestModel` <-
function(ANN=TRUE, CTA=TRUE, GAM=TRUE, GBM=TRUE,GLM=TRUE, MARS=TRUE, MDA=TRUE, RF=TRUE, SRE=TRUE, Bin.trans=TRUE, Filt.trans=TRUE, method='all')
{
    Th <- c('Kappa','TSS','Roc', 'all')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc'")
    
    
    #apply the function to the 3 possible transformation methods if all are selected
    if(method == 'all') for(k in 1:3) PredictionBestModel(ANN, CTA, GAM, GBM, GLM, MARS, MDA, RF, SRE, Bin.trans, Filt.trans, method=Th[k])   #runs the function alternatively for each method
                 
    #run the function for one method at a time                  
    else { if(Biomod.material$evaluation.choice[method]){
    
    
        NbSp <- Biomod.material$NbSpecies
        SpNames <- Biomod.material$species.names
        algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
        algo.c[names(which(!Biomod.material$algo.choice))] <- F  #switch off the models that are wanted but have not been trained    
    
        gg <- list()
    
        i <- 1
        while(i <= NbSp) {
            
            jj <- 3 #if(exists("DataEvalBIOMOD"))   jj <- 2   else    jj <- 3
   
            #load prob data from the pred directory
            eval(parse(text=paste("load('", getwd(), "/pred/Pred_", SpNames[i],"')", sep="")))
            sp.data <- eval(parse(text=paste("Pred_", SpNames[i], sep="")))
            
            
            #define an order to select models if there are equals : GAM > GLM > GBM > RF > CTA > MDA > MARS > ANN > SRE
            G <- matrix(0, ncol=9, nrow=Biomod.material$NbRun[i], dimnames=list(1:Biomod.material$NbRun[i], c("GAM", "GLM", "GBM", "RF", "CTA", "MDA", "MARS", "ANN", "SRE")))            
            #storing info on best model for each species
            g <- as.data.frame(matrix(0, nrow=Biomod.material$NbRun[i], ncol=7, dimnames=list(1:Biomod.material$NbRun[i], c("Best.Model", 'Cross.validation','indepdt.data','total.score','Cutoff','Sensitivity','Specificity'))))
            #storing matrices of outputs 
            MAT.bin <- MAT.filt <- MAT <- matrix(NA, nr=dim(sp.data)[1], nc=Biomod.material$NbRun[i], dimnames=list(1:dim(sp.data)[1], rep(NA,Biomod.material$NbRun[i])))      
    
    
            NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)   #considering the number of PA runs that were done   
            nbrep <- Biomod.material$NbRunEval + 1
    
       
            for(j in 1:NbPA){
                for(k in 1:nbrep){ 
                 
                    #writing the name to use for getting the right info in Evaluation.results lists
                    if(Biomod.material$NbRepPA == 0) nam <- "full" else nam <- paste("PA", j, sep="")
                    if(k!=1) nam <- paste(nam, "_rep", k-1, sep="")
                    #assign name to column of matrix                    
                    dimnames(MAT)[[2]][(j-1)*nbrep+k] <- nam 
                    rownames(g)[(j-1)*nbrep+k] <- nam   
                    nam <- paste(Biomod.material$species.names[i], nam, sep="_")
                    
                    #determine the best model
                    for(a in Biomod.material$algo[algo.c]) if(a != 'SRE') eval(parse(text=paste("G[(j-1)*nbrep+k, a] <- as.numeric(Evaluation.results.", method, "[[nam]][a,jj])", sep="")))
                    temp <- factor(which.max(G[(j-1)*nbrep+k,]), levels=seq(along=colnames(G)), labels=colnames(G))   
                                                                                                      
                    for(a in Biomod.material$algo[algo.c]){
                         if(a == temp){
                              g[(j-1)*nbrep+k, 2:7] <- eval(parse(text=paste("Evaluation.results.", method, sep="")))[[nam]][a,1:6]
                              MAT[,(j-1)*nbrep+k] <- sp.data[,a,k,j]
                              if(Bin.trans)  MAT.bin[,(j-1)*nbrep+k] <- BinaryTransformation(MAT[,(j-1)*nbrep+k], g[(j-1)*nbrep+k,5])
                              if(Filt.trans) MAT.filt[,(j-1)*nbrep+k] <- FilteringTransformation(MAT[,(j-1)*nbrep+k], g[(j-1)*nbrep+k,5])
                         }
                    }
                    
                } #nbrep k loop
            } #NbPA j loop
          
            colnames(MAT.bin) <- colnames(MAT.filt) <- colnames(MAT)     
            g[,1] <- factor(max.col(G), levels=seq(along=colnames(G)), labels=colnames(G))     
            gg[[SpNames[i]]] <- g

            #objects assignations and storage on hard disk           
            assign(paste("PredBestModelBy",method, "_", SpNames[i],sep=""), as.data.frame(MAT))
            eval(parse(text=paste("save(PredBestModelBy",method, "_", SpNames[i],", file='", getwd(), "/pred/PredBestModelBy",method, "_", SpNames[i],"')", sep="")))
            write.table(MAT, file=paste(getwd(),"/pred/PredBestModelBy",method, "_", SpNames[i],".txt", sep=""), row.names=F) 
            
            if(Bin.trans) {
            assign(paste("PredBestModelBy",method, "_", SpNames[i],"_Bin",sep=""), as.data.frame(MAT.bin))
            eval(parse(text=paste("save(PredBestModelBy",method, "_", SpNames[i],"_Bin, file='", getwd(), "/pred/PredBestModelBy",method, "_", SpNames[i],"_Bin')", sep="")))
            write.table(MAT.bin, file=paste(getwd(),"/pred/PredBestModelBy",method, "_", SpNames[i],"_Bin.txt", sep=""), row.names=F)
            }
            
            if(Filt.trans) {
            assign(paste("PredBestModelBy",method, "_", SpNames[i],"_Filt",sep=""), as.data.frame(MAT.filt))
            eval(parse(text=paste("save(PredBestModelBy",method, "_", SpNames[i],"_Filt, file='", getwd(), "/pred/PredBestModelBy",method, "_", SpNames[i],"_Filt')", sep="")))
            write.table(MAT.filt, file=paste(getwd(),"/pred/PredBestModelBy",method, "_", SpNames[i],"_Filt.txt", sep=""), row.names=F)
            }                            
                    
            i <- i + 1
        } #species i loop
    
        #saving the list of best models per method
        assign(paste("BestModelBy", method, sep=""), gg)
        eval(parse(text=paste("save(BestModelBy",method,", file='", getwd(), "/pred/BestModelBy",method,"')", sep="")))
    
    
    }}
}

