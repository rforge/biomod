`PredictionBestModel` <-
function(ANN=TRUE, CTA=TRUE, GAM=TRUE, GBM=TRUE,GLM=TRUE, MARS=TRUE, FDA=TRUE, RF=TRUE, SRE=TRUE, Bin.trans=TRUE, Filt.trans=TRUE, method='all')
{
    Th <- c('Kappa','TSS','Roc', 'all')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc'")
    
    
    #apply the function to the 3 possible transformation methods if all are selected
    if(method == 'all') for(k in 1:3) PredictionBestModel(ANN, CTA, GAM, GBM, GLM, MARS, FDA, RF, SRE, Bin.trans, Filt.trans, method=Th[k])   #runs the function alternatively for each method
                 
    #run the function for one method at a time                  
    else { if(Biomod.material$evaluation.choice[method]){
    
        NbSp <- Biomod.material$NbSpecies
        SpNames <- Biomod.material$species.names
        algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, FDA=FDA, RF=RF, SRE=SRE)
        algo.c[names(which(!Biomod.material$algo.choice))] <- FALSE  #switch off the models that are wanted but have not been trained    
    
        gg <- list()
 
        #storing arrays of outputs
        ARRAY.bin <- ARRAY.filt <- ARRAY <- array(NA, c(nrow(DataBIOMOD), max(Biomod.material$NbRun), Biomod.material$NbSpecies), dimnames=list(1:nrow(DataBIOMOD), rep(NA,max(Biomod.material$NbRun)), Biomod.material$species.names))
    
    
        i <- 1
        while(i <= NbSp) {
            
            jj <- 1 #selecting according to the cross validation values     #if(exists("DataEvalBIOMOD"))   jj <- 2   else    jj <- 1
   
            #load prob data from the pred directory
            eval(parse(text=paste("load('", getwd(), "/pred/Pred_", SpNames[i],"')", sep="")))
            sp.data <- eval(parse(text=paste("Pred_", SpNames[i], sep="")))
            
            
            #define an order to select models if there are equals : GAM > GLM > GBM > RF > CTA > FDA > MARS > ANN > SRE
            G <- matrix(0, ncol=9, nrow=Biomod.material$NbRun[i], dimnames=list(1:Biomod.material$NbRun[i], c("GAM", "GLM", "GBM", "RF", "CTA", "FDA", "MARS", "ANN", "SRE")))            
            #storing info on best model for each species
            g <- as.data.frame(matrix(0, nrow=Biomod.material$NbRun[i], ncol=7, dimnames=list(1:Biomod.material$NbRun[i], c("Best.Model", 'Cross.validation','indepdt.data','total.score','Cutoff','Sensitivity','Specificity'))))
                 
    
            NbPA <- Biomod.material$NbRun[i] / (Biomod.material$NbRunEval+1)   #considering the number of PA runs that were done   
            nbrep <- Biomod.material$NbRunEval + 1
       
            for(j in 1:NbPA){
                for(k in 1:nbrep){ 
                 
                    #writing the name to use for getting the right info in Evaluation.results lists
                    if(Biomod.material$NbRepPA == 0) nam <- "full" else nam <- paste("PA", j, sep="")
                    if(k!=1) nam <- paste(nam, "_rep", k-1, sep="")
                    #assign name to column of matrix                    
                    dimnames(ARRAY)[[2]][(j-1)*nbrep+k] <- nam 
                    rownames(g)[(j-1)*nbrep+k] <- nam   
                    nam <- paste(Biomod.material$species.names[i], nam, sep="_")
                    
                    #determine the best model
                    for(a in Biomod.material$algo[algo.c]) if(a != 'SRE') eval(parse(text=paste("G[(j-1)*nbrep+k, a] <- as.numeric(Evaluation.results.", method, "[[nam]][a,jj])", sep="")))
                    temp <- factor(which.max(G[(j-1)*nbrep+k,]), levels=seq(along=colnames(G)), labels=colnames(G))   
                                                                                                      
                    for(a in Biomod.material$algo[algo.c]){
                         if(a == temp){
                              g[(j-1)*nbrep+k, 2:7] <- eval(parse(text=paste("Evaluation.results.", method, sep="")))[[nam]][a,1:6]

                              #get the best model and fill with NAs to have predictions of the same length across species
                              sdata <- sp.data[,a,k,j]
                              sdataNA <- c(sdata, rep(NA, nrow(DataBIOMOD) - length(sdata)))
                              sdata.bin <- BinaryTransformation(sdata, g[(j-1)*nbrep+k,5])
                              sdata.bin <- c(sdata.bin, rep(NA, nrow(DataBIOMOD) - length(sdata)))
                              sdata.filt <- FilteringTransformation(sdata, g[(j-1)*nbrep+k,5])
                              sdata.filt <- c(sdata.filt, rep(NA, nrow(DataBIOMOD) - length(sdata)))

                              ARRAY[ ,(j-1)*nbrep+k, i] <- sdataNA
                              if(Bin.trans)  ARRAY.bin[ ,(j-1)*nbrep+k, i] <- sdata.bin
                              if(Filt.trans) ARRAY.filt[ ,(j-1)*nbrep+k, i] <- sdata.filt

                         }
                    }
                    
                } #nbrep k loop
            } #NbPA j loop
          
            dimnames(ARRAY.bin)[[2]] <- dimnames(ARRAY.filt)[[2]] <- dimnames(ARRAY)[[2]]     
            g[,1] <- factor(max.col(G), levels=seq(along=colnames(G)), labels=colnames(G))     
            gg[[SpNames[i]]] <- g                         
                    
            i <- i + 1
        } #species i loop
    

        #objects assignations and storage on hard disk           
        assign(paste("PredBestModelBy",method,sep=""), ARRAY)
        eval(parse(text=paste("save(PredBestModelBy",method,", file='", getwd(), "/pred/PredBestModelBy",method, "', compress='xz')", sep="")))
        write.table(ARRAY, file=paste(getwd(),"/pred/PredBestModelBy",method, ".txt", sep=""), row.names=FALSE) 
        
        if(Bin.trans) {
        assign(paste("PredBestModelBy",method, "_Bin",sep=""), ARRAY.bin)
        eval(parse(text=paste("save(PredBestModelBy",method, "_Bin, file='", getwd(), "/pred/PredBestModelBy",method, "_Bin', compress='xz')", sep="")))
        write.table(ARRAY.bin, file=paste(getwd(),"/pred/PredBestModelBy",method, "_Bin.txt", sep=""), row.names=FALSE)
        }
        
        if(Filt.trans) {
        assign(paste("PredBestModelBy",method, "_Filt",sep=""), ARRAY.filt)
        eval(parse(text=paste("save(PredBestModelBy",method, "_Filt, file='", getwd(), "/pred/PredBestModelBy",method, "_Filt', compress='xz')", sep="")))
        write.table(ARRAY.filt, file=paste(getwd(),"/pred/PredBestModelBy",method, "_Filt.txt", sep=""), row.names=FALSE)
        }
    
        #saving the list of best models per method
        assign(paste("BestModelBy", method, sep=""), gg)
        eval(parse(text=paste("save(BestModelBy",method,", file='", getwd(), "/pred/BestModelBy",method,"', compress='xz')", sep="")))
    
    }}
}

