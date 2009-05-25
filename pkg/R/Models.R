`Models` <-
function(GLM=FALSE, TypeGLM="simple", Test="AIC", GBM=FALSE, No.trees= 2000, GAM=FALSE, Spline=3, CTA=FALSE, CV.tree=50, ANN=FALSE, CV.ann=5, SRE=FALSE, Perc025=FALSE, Perc05=FALSE, MDA=FALSE, MARS=FALSE, RF=FALSE,
NbRunEval=1, DataSplit=100, NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=NULL, Yweights=NULL, VarImport=0, 
Roc=FALSE, Optimized.Threshold.Roc=FALSE, Kappa=FALSE, TSS=FALSE, KeepPredIndependent=FALSE)
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
    
    #checking possible mistakes in the argument selections
    if(!exists("DataBIOMOD")) stop("Initial.State should be run first") 
    if(!Roc && !TSS && !Kappa) stop("at least one evaluation technique (Roc, TSS or Kappa) must be selected \n") 
    if(Roc != T && Optimized.Threshold.Roc != F) stop("Roc must be TRUE to derive optimized threshold value")
    if(DataSplit < 50) cat("Warning : You choose to allocate more data to evaluation than to calibration of your model \n Make sure you really wanted to do that. \n")
      
      
    #check that the weight matrix (if any was specified) was entered correctly:
    if(!is.null(Yweights)){
        if(is.null(dim(Yweights))) Yweights <- as.data.frame(Yweights)
        if(ncol(Yweights) != Biomod.material$NbSpecies) stop("The number of 'Weight' columns does not match the number of species. Simulation cannot proceed.")
        if(nrow(Yweights) != nrow(DataBIOMOD)) stop("The number of 'Weight' rows does not match with the input calibration data. Simulation cannot proceed.")
    }
    
    #create the directories in which various objects will be stored (models, predictions and projection). The projection directories are created in the Projection() function.
    dir.create(paste(getwd(), "/models", sep=""), showWarnings=F)
    dir.create(paste(getwd(), "/pred", sep=""), showWarnings=F)
  
    #create usefull vectors for condensing code   
    Biomod.material[["algo"]] <- c("ANN","CTA","GAM","GBM","GLM","MARS","MDA","RF","SRE")
    Biomod.material[["algo.choice"]] <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    Biomod.material[["evaluation.choice"]] <- c(Roc=Roc, Kappa=Kappa, TSS=TSS)
    Biomod.material[["NbRunEval"]] <- NbRunEval
    Biomod.material[["NbRepPA"]] <- NbRepPA
    assign("Biomod.material", Biomod.material, pos=1)
 
    #run the pseudo.absence function once for each species, keeping all the absences possible each time
    #only the row names are stored in PA.data
    if(NbRepPA > 0){
        PA.data <- vector('list', Biomod.material$NbSpecies)
        names(PA.data) <- Biomod.material$species.names  
        Biomod.PA.sample <- PA.data
                 
        for(i in 1:Biomod.material$NbSpecies){
            PA.data[[i]] <- pseudo.abs(coor=coor, status=DataBIOMOD[,Biomod.material$NbVar+i], env=DataBIOMOD[,1:Biomod.material$NbVar], 
              strategy=strategy, distance=distance, nb.points=NULL, add.pres=TRUE, species.name= 'SpNoName', create.dataset=FALSE)                                        
        }
        assign("Biomod.PA.data",PA.data, pos=1)
    }
 
    #defining evaluation runs
    if(NbRunEval==0){
        DataSplit <- 100
        if(!exists("DataEvalBIOMOD")) cat("\n\n Warning : The models will be evaluated on the calibration data only (NbRunEval=0 and no independent data) \n\t it could lead to over-optimistic predictive performances. \n\n")
    }
    if(DataSplit==100) NbRunEval <- 0   
            
    #create list to store results from some specific models   
    if(GBM | ANN | RF | MARS | MDA) assign('Models.information', list(), pos=1)
    
    #create matrix to store variable importance results
    if(VarImport != 0){
        VI <- vector('list', Biomod.material$NbSpecies)
        names(VI) <- Biomod.material$species.names
        for(i in 1:Biomod.material$NbSpecies) VI[[i]] <- data.frame(matrix(NA, nrow=sum(Biomod.material[["algo.choice"]]), ncol=Biomod.material$NbVar, dimnames=list(Biomod.material$algo[Biomod.material$algo.choice], names(DataBIOMOD)[1:Biomod.material$NbVar])))
        assign("VarImportance", VI, pos=1)
    } else assign("VarImportance", NA, pos=1)  
    
    
    #create output objects to store the results of the evaluation procedures 
    if(NbRepPA == 0) { PAs <- rep("full", (NbRunEval+1)*Biomod.material$NbSpecies)
    } else{ PAs <- c() ; for(j in 1:NbRepPA) PAs <- c(PAs, paste("PA", j, sep="")) 
           PAs <- rep(rep(PAs, each=NbRunEval+1), Biomod.material$NbSpecies)
    }
    reps <- ""
    if(NbRunEval != 0) for(j in 1:NbRunEval) reps <- c(reps, paste("rep", j, sep=""))
    sp <- rep(Biomod.material$species.names, each=length(reps))
    if(NbRepPA != 0) sp <- rep(sp, each=NbRepPA)
    reps <- rep(reps, Biomod.material$NbSpecies)
    if(NbRepPA != 0) reps <- rep(reps, NbRepPA)
    
    g <- vector('list',length(reps))
    gnames <- rep(NA, length(g))
    for(i in 1:length(g)) if(reps[i]!="") gnames[i] <- paste(sp[i],"_", PAs[i], "_", reps[i], sep="") else gnames[i] <- paste(sp[i],"_", PAs[i], sep="")
                      
    for(i in 1:length(g)) g[[i]] <- data.frame(matrix(NA, nrow=sum(Biomod.material$algo.choice), ncol=6, dimnames=list(Biomod.material$algo[Biomod.material$algo.choice], c('Cross.validation','indepdt.data','total.score','Cutoff','Sensitivity','Specificity'))))
    names(g) <- gnames
    if(Roc) assign("Evaluation.results.Roc", g, pos=1) else assign("Evaluation.results.Roc", NA, pos=1)
    if(Kappa) assign("Evaluation.results.Kappa", g, pos=1) else assign("Evaluation.results.Kappa", NA, pos=1)
    if(TSS) assign("Evaluation.results.TSS", g, pos=1) else assign("Evaluation.results.TSS", NA, pos=1)


    #start species loop
    i <- 1
    while(i <= Biomod.material$NbSpecies) {
        cat("#####\t\t\t", Biomod.material$species.names[i], "\t\t\t#####\n")
        assign("i", i, pos= 1)
               
        if(exists("Models.information")){
           Models.information[[Biomod.material$species.names[i]]] <- list()
           assign('Models.information', Models.information, pos=1) 
        }  
        
        #check the number of absences wanted compared to those available
        #for the case where the number of absences wanted is higher than those available -> set NbRepPA.pos to 1 and do a single PA run with all the absences available
        NbRepPA.pos <- 1
        if(is.null(nb.absences)) nb.absences <- nrow(DataBIOMOD)
        nb.absences.pos <- nb.absences
            
        if(NbRepPA != 0){
            NbRepPA.pos <- NbRepPA
            nbpres <- sum(DataBIOMOD[,Biomod.material$NbVar+i])
            if(nb.absences > (length(Biomod.PA.data[[i]]) - nbpres)) {  #if not enough absences -> only one PA run with nb.absences set to max
                NbRepPA.pos <- 1                                           
                nb.absences.pos <- length(Biomod.PA.data[[i]]) - nbpres
            }  
        }          
        
            
        #constructing the storing array for the species considering NbRepPA.pos, NbRunEval, nb.absences.pos
        reps <- c()
        if(NbRunEval != 0) for(j in 1:NbRunEval) reps <- c(reps, paste("rep", j, sep="")) 
        PAs <- c()
        if(NbRepPA.pos == 1 && NbRepPA==0) PAs <- "no.PA" else for(j in 1:NbRepPA.pos) PAs <- c(PAs, paste("PA", j, sep=""))
    
        if(NbRepPA == 0) ndata <- nrow(DataBIOMOD) 
        if(NbRepPA != 0) ndata <- nb.absences.pos + nbpres
        ARRAY <- array(NA, c(ndata, 9, NbRunEval+1, NbRepPA.pos), dimnames=list(1:ndata, Biomod.material$algo, c("total.data", reps), PAs))
        if(exists("DataEvalBIOMOD") && KeepPredIndependent) ARRAY.ind <- array(NA, c(nrow(DataEvalBIOMOD), 9, NbRunEval+1, NbRepPA.pos), dimnames=list(1:nrow(DataEvalBIOMOD), Biomod.material[["algo"]], c("total.data", reps), PAs))

        assign("Array", ARRAY, pos=1)
        Ids <- data.frame(matrix(0, nrow=ceiling(ndata*(DataSplit/100)), ncol=NbRunEval))
        

        for(pa in 1:NbRepPA.pos){

            assign("pa", pa, pos=1)
            if(NbRepPA != 0) cat("#####\t\t   pseudo-absence run", pa, "       \t\t#####\n")           
                
            #defining the data (as lines to take from DataBIOMOD) to be used for calibration
            if(NbRepPA == 0) PA.samp <- 1:nrow(DataBIOMOD)
            else {
                absamp <- sort(sample((nbpres+1):length(Biomod.PA.data[[i]]), nb.absences.pos))
                PA.samp <- sort(Biomod.PA.data[[i]][c(1:nbpres,absamp)])
                
                Biomod.PA.sample[[i]][[paste("PA", pa, sep="")]] <- PA.samp  #storing the lines selected for each PA run 
            }        
                
            #defining the Ids to be selected for the evaluation runs
            if(NbRunEval != 0) for(j in 1:NbRunEval) Ids[,j] <- sort(SampleMat2(DataBIOMOD[PA.samp,(Biomod.material$NbVar+i)], DataSplit/100)$calibration)
                   
            for(a in Biomod.material$algo[Biomod.material$algo.choice]){
                
                  g <- list(Biomod.Models(a, Ids, PA.samp, TypeGLM, Test, No.trees, CV.tree, CV.ann, Perc025, Perc05, NbRunEval, Spline, DataSplit,
                            Yweights, Roc, Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, VarImport))
                            
                  if(exists("DataEvalBIOMOD") && KeepPredIndependent) ARRAY.ind[,a,1,pa] <- predind
                  
                  if(a=='ANN' | a=='GBM' | a=='RF' | a=='MARS' | a=='MDA'){
                      Models.information[[Biomod.material$species.names[i]]][[a]][[paste("PA", pa, sep="")]] <- g         # Models.information[[Biomod.material$species.names[i]]][[paste(a, "_PA", pa, sep="")]]
                      assign('Models.information', Models.information, pos=1)
                  } 
            }    
        }
        
        #save the prediction for that species
        assign(paste("Pred_", Biomod.material$species.names[i], sep=""), Array)
        eval(parse(text=paste("save(Pred_",Biomod.material$species.names[i],", file='", getwd(),"/pred/Pred_",Biomod.material$species.names[i], "')", sep="")))     
  
        if(exists("DataEvalBIOMOD") && KeepPredIndependent) {
            assign(paste("Pred_",Biomod.material$species.names[i], "_indpdt", sep=""), ARRAY.ind)
            eval(parse(text=paste("save(Pred_",Biomod.material$species.names[i], "_indpdt, file='",getwd(),"/pred/Pred_",Biomod.material$species.names[i], "_indpdt')", sep="")))
        }
        
        #store the total number of models produced for that species. It will be used by Projection()
        if(NbRepPA!=0) Biomod.material[["NbRun"]][i] <- dim(ARRAY)[3]*dim(ARRAY)[4] else Biomod.material[["NbRun"]][i] <- dim(ARRAY)[3]
  
        i <- i + 1
    }
    
    rm(i, pa, calib.lines, Array, pos=1)
    assign("Biomod.material", Biomod.material, pos=1)
    if(NbRepPA != 0) assign('Biomod.PA.sample', Biomod.PA.sample, pos=1)
    
    if(exists("Models.information")){
      save(Evaluation.results.Roc, Evaluation.results.TSS, Evaluation.results.Kappa, VarImportance, Biomod.material, Models.information,
        file=if(Biomod.material[["NbSpecies"]]==1) paste(Biomod.material[["species.names"]], "_run.RData", sep="") else file='Biomod_run.RData')
    }
    else {
      save(Evaluation.results.Roc, Evaluation.results.TSS, Evaluation.results.Kappa, VarImportance, Biomod.material, 
        file=if(Biomod.material[["NbSpecies"]]==1) paste(Biomod.material[["species.names"]], "_run.RData", sep="") else file='Biomod_run.RData')
    }
  
}

