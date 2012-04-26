`Models` <-
function(GLM=FALSE, TypeGLM="simple", Test="AIC", GBM=FALSE, No.trees= 5000, GAM=FALSE, Spline=3, CTA=FALSE, CV.tree=50, ANN=FALSE, CV.ann=5, SRE=FALSE, quant=0.025, FDA=FALSE, MARS=FALSE, RF=FALSE,
NbRunEval=1, DataSplit=100, NbRepPA=0, strategy="sre", coor=NULL, distance=0, nb.absences=NULL, Yweights=NULL, VarImport=0, 
Roc=FALSE, Optimized.Threshold.Roc=FALSE, Kappa=FALSE, TSS=FALSE, KeepPredIndependent=FALSE)
{
    require(nnet, quietly=TRUE)
    require(rpart, quietly=TRUE)
#     require(Hmisc, quietly=TRUE)
#     require(Design, quietly=TRUE)
    require(MASS, quietly=TRUE)
    require(gbm, quietly=TRUE)
    require(mda, quietly=TRUE)
    require(randomForest, quietly=TRUE)
    require(gam, quietly=TRUE)	
    

    
    #checking possible mistakes in the argument selections
    if(!exists("DataBIOMOD")) stop("Initial.State should be run first in order to procede")
    if(!any(GAM,GBM,GLM,RF,FDA,MARS,SRE,ANN,CTA)) stop("No models were selected \n") 
    if(!any(Roc,Kappa,TSS)) stop("At least one evaluation technique (Roc, TSS or Kappa) must be selected \n") 
    if(Roc != TRUE && Optimized.Threshold.Roc != FALSE) stop("Roc must be TRUE to derive optimized threshold value")
    if(DataSplit < 50) cat("Warning : You choose to allocate more data to evaluation than to calibration of your model (DataSplit<50) \n Make sure you really wanted to do that. \n") 
    if(quant>=0.5 | quant<0) stop("\n settings in 'quant' should be a value between 0 and 0.5 ")  
      
      
    #Check that the weight matrix was entered correctly with the pseudo.abs options    
    if(!is.null(Yweights)){
       if(is.null(dim(Yweights))) Yweights <- as.data.frame(Yweights)
       if(ncol(Yweights) != Biomod.material$NbSpecies) cat("\n Warning : Yweights and input data differ in length (nb columns). Check if this is wanted.")
       if(nrow(Yweights) != nrow(DataBIOMOD)) stop("The number of 'Weight' rows does not match with the input calibration data. Simulation cannot proceed.")
    }
    assign("isnullYweights", is.null(Yweights), pos=1)                                                                  #To keep track of Yweights state at origin (user's input)
    
    #checking var importance args
    if(VarImport > 0 && Biomod.material$NbVar < 2){
      cat("\nVar Importance calculation was automaticly switched off because of only one explanatory variable given\n")
      VarImport <- 0
    }
    
    if(NbRepPA!=0 && is.null(Yweights)) Yweights <- matrix(NA, nc=Biomod.material$NbSpecies, nr=nrow(DataBIOMOD))
   
    
    #create the directories in which various objects will be stored (models, predictions and projection). The projection directories are created in the Projection() function.
    dir.create(paste(getwd(), "/models", sep=""), showWarnings=FALSE)                   
    dir.create(paste(getwd(), "/pred", sep=""), showWarnings=FALSE)
    if(any(MARS, FDA, ANN)) dir.create(paste(getwd(), "/models/rescaling_models", sep=""), showWarnings=FALSE) 
  
  
    #switch SRE and MARS off if one of the variables is a non numeric
    Nbcat <- rep(0, Biomod.material$NbVar)
    for(i in 1:Biomod.material$NbVar) Nbcat[i] <- is.factor(DataBIOMOD[,i]) 
       
    if(sum(Nbcat) > 0 && MARS){ MARS <- FALSE ; cat(paste("MARS model was shut down, it cannot run on factorial variables : ", paste(Biomod.material$VarNames[Nbcat==TRUE], collapse=" "), "\n", sep="")) }
    if(sum(Nbcat) > 0 && SRE){ SRE <- FALSE ; cat(paste("SRE model was shut down, it cannot run on factorial variables : ", paste(Biomod.material$VarNames[Nbcat==TRUE], collapse=" "), "\n", sep="")) }
     
    #defining evaluation runs
    if(NbRunEval==0){
        DataSplit <- 100
        if(!exists("DataEvalBIOMOD")) cat("\n\n Warning : The models will be evaluated on the calibration data only (NbRunEval=0 and no independent data) \n\t it could lead to over-optimistic predictive performances. \n\n")
    }
    if(DataSplit==100) NbRunEval <- 0 
            
  
    #create usefull vectors for condensing code   
    Biomod.material[["algo"]] <- c("ANN","CTA","GAM","GBM","GLM","MARS","FDA","RF","SRE")
    Biomod.material[["algo.choice"]] <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, FDA=FDA, RF=RF, SRE=SRE)
    Biomod.material[["evaluation.choice"]] <- c(Roc=Roc, Kappa=Kappa, TSS=TSS)
    Biomod.material[["NbRunEval"]] <- NbRunEval
    Biomod.material[["NbRepPA"]] <- NbRepPA
    Biomod.material[["calibration.failures"]] <- NA
    assign("BM", Biomod.material, pos=1)
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
 
 
      #create list to store best.iter for GBM (needed for projection -> n.trees argument)  
    if(GBM){ GBM.perf <- list() ; GBMP <- c() }
    
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
    


    #modelling summary to be writen in console
    if(NbRepPA==0) NbRepPA.pos <- 1 else NbRepPA.pos=NbRepPA  #need to be given for calculation of total number of runs (bottom)
    cat(paste(
    "\n----------------------------------- \n",
    "Modelling summary \n",
    "----------------------------------- \n",
    "Number of species modelled : \t\t" , Biomod.material$NbSpecies, "\n",
    paste(Biomod.material$species.names, collapse=", "), "\n\n",
    
    "numerical variables : \t\t\t" , paste(Biomod.material$VarNames[Nbcat==FALSE], collapse=", "), "\n",  
    if(sum(Nbcat) > 0) "factorial variables : \t\t\t" , paste(Biomod.material$VarNames[Nbcat==TRUE], collapse=", "), "\n\n",
    
    "number of evaluation repetitions : \t" , NbRunEval, "\n",
    "number of pseudo-absences runs : \t" , NbRepPA, "\n",  
    "models selected : \t\t\t" ,  paste(Biomod.material$algo[Biomod.material$algo.choice], collapse=", "), "\n",
    "total number of model runs :  \t\t" , (NbRepPA.pos) * (NbRunEval+1) * Biomod.material$NbSpecies * sum(Biomod.material$algo.choice), "\n",
    "----------------------------------- \n\n\n",
    sep=""))




    #-----------------start species loop-------------------#
    i <- 1
    while(i <= Biomod.material$NbSpecies) {
        cat("#####\t\t\t", Biomod.material$species.names[i], "\t\t\t#####\n")
        assign("i", i, pos= 1)
               
        if(GBM) GBM.perf[[Biomod.material$species.names[i]]] <- list()                #Store the best.iter per species
        
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
              

        #defining the weights to be awarded for that species, considering the number of absences available.
        #we keep the same format as DataBIOMOD to make it easier to call the corresponding lines in Biomod.models. For that
        #reason, we award the same weight to all the absences even though it doesn't sum up to a 0.5 prevalence (but lower)
        #it will give 0.5 when considering the number of absences, hence the lines taken for calibration. 
        if(NbRepPA!=0 && isnullYweights){   
            Yweights[which(DataBIOMOD[,Biomod.material$NbVar+i]==1),i] <- 1
            Yweights[which(DataBIOMOD[,Biomod.material$NbVar+i]==0),i] <- nbpres/nb.absences.pos    
        }    
                
            
        #constructing the storing array for the species considering NbRepPA.pos(dim4), NbRunEval(dim3), nb.absences.pos(dim2)
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
        if(GBM) GBM.perf[[Biomod.material$species.names[i]]] <- list()
        
        

        for(pa in 1:NbRepPA.pos){
            assign("pa", pa, pos=1)
            if(NbRepPA != 0) cat("#####\t\t   pseudo-absence run", pa, "       \t\t#####\n")               
        
                
            #defining the data (as lines to take from DataBIOMOD) to be used for calibration
            #to constitute the data for that PA run
            if(NbRepPA != 0){ 
                absamp <- sort(sample((nbpres+1):length(Biomod.PA.data[[i]]), nb.absences.pos))
                PA.samp <- sort(Biomod.PA.data[[i]][c(1:nbpres,absamp)])
                Biomod.PA.sample[[i]][[paste("PA", pa, sep="")]] <- PA.samp  #storing the lines selected for each PA run   
            } else PA.samp <- 1:nrow(DataBIOMOD) 
                
            #defining the Ids to be selected for the evaluation runs
            if(NbRunEval != 0) for(j in 1:NbRunEval) Ids[,j] <- sort(.SampleMat2(DataBIOMOD[PA.samp,(Biomod.material$NbVar+i)], DataSplit/100)$calibration)      
                   
                   
            #Run Biomod.models       
            for(a in Biomod.material$algo[Biomod.material$algo.choice]){
                  .Biomod.Models(a, Ids, PA.samp, TypeGLM, Test, No.trees, CV.tree, CV.ann, quant, NbRunEval, Spline, DataSplit,
                            Yweights, Roc, Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, VarImport)
                  if(exists("DataEvalBIOMOD") && KeepPredIndependent) ARRAY.ind[,a,1,pa] <- predind
                   
            }
            if(GBM) GBMP <- c(GBMP, GBM.list)  #add info from each PA to preexisting list (total=one species)
        }
    
        
        #save the prediction for that species
        assign(paste("Pred_", Biomod.material$species.names[i], sep=""), Array)
        eval(parse(text=paste("save(Pred_",Biomod.material$species.names[i],", file='", getwd(),"/pred/Pred_",Biomod.material$species.names[i], "', compress='xz')", sep="")))     
  
        if(exists("DataEvalBIOMOD") && KeepPredIndependent) {
            assign(paste("Pred_",Biomod.material$species.names[i], "_indpdt", sep=""), ARRAY.ind)
            eval(parse(text=paste("save(Pred_",Biomod.material$species.names[i], "_indpdt, file='",getwd(),"/pred/Pred_",Biomod.material$species.names[i], "_indpdt', compress='xz')", sep="")))
        }
        
        #store the total number of models produced for that species. It will be used by Projection()
        if(NbRepPA!=0) BM[["NbRun"]][i] <- dim(ARRAY)[3]*dim(ARRAY)[4] else BM[["NbRun"]][i] <- dim(ARRAY)[3]
  
        if(GBM) GBM.perf[[BM$species.names[i]]] <- GBMP  
  
        i <- i + 1
    }
    
    #suppress useless objects
    rm(i, pa, Array, calib.lines, pos=1)
    if(exists("DataEvalBIOMOD") && KeepPredIndependent) rm(predind, pos=1)

    if(length(BM[["calibration.failures"]]) > 1) BM[["calibration.failures"]] <- BM[["calibration.failures"]][-1] 
    assign("Biomod.material", BM, pos=1)
    if(GBM) assign("GBM.perf", GBM.perf, pos=1)
    if(NbRepPA != 0) assign('Biomod.PA.sample', Biomod.PA.sample, pos=1)
    
    
    #save the history and workspace
    if(Biomod.material[["NbSpecies"]]==1) filename <- paste(Biomod.material[["species.names"]], "_run", sep="") else filename <- 'Biomod_run' 
    save.image(paste(filename, ".RData", sep=""))
    #savehistory(paste(filename, ".Rhistory", sep=""))
    
    #Final notice, runs are finished
    cat("\n\n--------- \n completed \n\n\n")  
    
    # If one or more of the evaluation runs failed, we report failure to user. This message is important.
    if(!is.na(BM$calibration.failures[1])) cat(paste("WARNING: the following repetition models failed : \n", BM$calibration.failures, "\n This might indicate serious problem in your data. \n\n", sep=""))
    rm(BM, pos=1)
        
} 
 
 
 