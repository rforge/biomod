`Models` <-
function(GLM=FALSE, TypeGLM=simple, Test=AIC, GBM=FALSE, No.trees= 2000, GAM=FALSE, Spline=3, CTA=FALSE, CV.tree=50,
ANN=FALSE, CV.ann=5, SRE=FALSE, Perc025=FALSE, Perc05=FALSE, MDA=FALSE, MARS=FALSE, RF=FALSE, NbRunEval=1, DataSplit=100, Yweights=NULL, Roc=FALSE,
Optimized.Threshold.Roc=FALSE, Kappa=FALSE, TSS=FALSE, KeepPredIndependent=FALSE, VarImport=0)
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
    
    if(!exists("DataBIOMOD")) stop("Initial.State should be run first") 
    if(!Roc && !TSS && !Kappa) stop("at least one evaluation technique (Roc, TSS or Kappa) must be selected \n") 
    if(Roc != T && Optimized.Threshold.Roc != F) stop("Roc must be TRUE to derive optimized threshold value")
    if(DataSplit < 50) cat("Warning : You choose to allocate more data to evaluation than to calibration of your model \n Make sure you really wanted to do that. \n")
      
    
    # Check that the weight matrix (if any was specified) was entered correctly:
    if(!is.null(Yweights)){
        if(is.null(dim(Yweights))) Yweights <- as.data.frame(Yweights)
        if(ncol(Yweights) != Biomod.material[["NbSpecies"]]) stop("The number of 'Weight' columns does not match the number of species. Simulation cannot proceed.")
        if(nrow(Yweights) != nrow(DataBIOMOD)) stop("The number of 'Weight' rows does not match with the input calibration data. Simulation cannot proceed.")
    }
    
    dir.create(paste(getwd(), "/models", sep=""), showWarnings=F)
    dir.create(paste(getwd(), "/pred", sep=""), showWarnings=F)
  
    #create usefull vectors for condensing code   
    Biomod.material[["algo"]] <- c("ANN","CTA","GAM","GBM","GLM","MARS","MDA","RF","SRE")
    Biomod.material[["algo.choice"]] <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    Biomod.material[["evaluation.choice"]] <- c(Roc=Roc, TSS=TSS, Kappa=Kappa)
    assign("Biomod.material", Biomod.material, pos=1)
 
    if(NbRunEval==0){
        DataSplit <- 100
        if(!exists("DataEvalBIOMOD")) cat("\n\n Warning : The models will be evaluated on the calibration data only (NbRunEval=0) \n\t it could lead to over-optimistic predictive performances. \n\n")
    }
    if(DataSplit==100) NbRunEval <- 1
    Ids <- data.frame(matrix(0, nrow=ceiling(nrow(DataBIOMOD)*(DataSplit/100)), ncol=NbRunEval))
       
    #Create list to store results from each model  
    if(GBM | ANN | RF | MARS | MDA) assign('Models.information', list(), pos=1)
    
    if(VarImport != 0){
        VI <- vector('list', Biomod.material[["NbSpecies"]])
        names(VI) <- Biomod.material[["species.names"]]
        for(i in 1:Biomod.material[["NbSpecies"]]) VI[[i]] <- data.frame(matrix(NA, nrow=sum(Biomod.material[["algo.choice"]]), ncol=Biomod.material[["NbVar"]], dimnames=list(Biomod.material[["algo"]][Biomod.material[["algo.choice"]]], names(DataBIOMOD)[1:Biomod.material[["NbVar"]]])))
        assign("VarImportance", VI, pos=1)
    } else assign("VarImportance", NA, pos=1)
    
    g <- vector('list',Biomod.material[["NbSpecies"]])
    names(g) <- Biomod.material[["species.names"]]
    for(i in 1:Biomod.material[["NbSpecies"]]) g[[i]] <- data.frame(matrix(NA, nrow=sum(Biomod.material[["algo.choice"]]), ncol=6, dimnames=list(Biomod.material[["algo"]][Biomod.material[["algo.choice"]]], c('Cross.validation','indepdt.data','Final.model','Cutoff','Sensitivity','Specificity'))))
    if(Roc) assign("Evaluation.results.Roc", g, pos=1) else assign("Evaluation.results.Roc", NA, pos=1)
    if(Kappa) assign("Evaluation.results.Kappa", g, pos=1) else assign("Evaluation.results.Kappa", NA, pos=1)
    if(TSS) assign("Evaluation.results.TSS", g, pos=1) else assign("Evaluation.results.TSS", NA, pos=1)


    i <- 1
    while(i <= Biomod.material[["NbSpecies"]]) {
        cat("#####\t\t\t", Biomod.material[["species.names"]][i], "\t\t\t#####\n")
        assign("i", i, pos= 1)
        mat <- matrix(NA, nr=nrow(DataBIOMOD), nc=9, dimnames=list(1:nrow(DataBIOMOD), Biomod.material[["algo"]]))
        if(exists("DataEvalBIOMOD") && KeepPredIndependent) mat.ind <- matrix(NA, nr=nrow(DataEvalBIOMOD), nc=9, dimnames=list(1:nrow(DataEvalBIOMOD), Biomod.material[["algo"]]))
       
        if(exists("Models.information")){
           Models.information[[Biomod.material[["species.names"]][i]]] <- list()
           assign('Models.information', Models.information, pos=1) 
        }
        
        if(NbRunEval==1 && DataSplit==100) { Ids <- data.frame(matrix(seq(1:nrow(DataBIOMOD)), ncol=1))
        } else for(j in 1:NbRunEval) Ids[,j] <- SampleMat2(DataBIOMOD[,Biomod.material[["NbVar"]]+i], DataSplit/100)$calibration
        
        for(a in Biomod.material[["algo"]][Biomod.material[["algo.choice"]]]){
            
              g <- list(Biomod.Models(a, Ids, TypeGLM, Test, No.trees, CV.tree, CV.ann, Perc025, Perc05, NbRunEval, Spline, DataSplit,
                        Yweights, Roc, Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, VarImport))
          
              mat[,a] <- g.pred
              if(exists("DataEvalBIOMOD") && KeepPredIndependent) mat.ind[,a] <- predtest
          
              if(a=='ANN' | a=='GBM' | a=='RF' | a=='MARS' | a=='MDA'){
                  Models.information[[Biomod.material[["species.names"]][i]]][a] <- g 
                  assign('Models.information', Models.information, pos=1)
              }
        }  
       
        assign(paste("Pred_",Biomod.material[["species.names"]][i], sep=""), mat)
        write.table(mat, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],".txt", sep=""), row.names=F)
        eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],", file='",getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i], "')", sep="")))
  
        if(exists("DataEvalBIOMOD") && KeepPredIndependent) {
            assign(paste("Pred_",Biomod.material[["species.names"]][i], "_indpdt", sep=""), mat.ind)
            write.table(mat, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i], "_indpdt.txt", sep=""), row.names=F)
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i], "_indpdt, file='",getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i], "_indpdt')", sep="")))
        }
  
        i <- i + 1
    }
    rm(sp, i, g.pred, pos=1)
    if(exists("Models.information")){
      save(Evaluation.results.Roc, Evaluation.results.TSS, Evaluation.results.Kappa, VarImportance, Biomod.material, Models.information,
        file=if(Biomod.material[["NbSpecies"]]==1) paste(Biomod.material[["species.names"]], "_run.RData", sep="") else file='Biomod_run.RData')
    }
    else {
      save(Evaluation.results.Roc, Evaluation.results.TSS, Evaluation.results.Kappa, VarImportance, Biomod.material,
        file=if(Biomod.material[["NbSpecies"]]==1) paste(Biomod.material[["species.names"]], "_run.RData", sep="") else file='Biomod_run.RData')
    }
 
  
}

