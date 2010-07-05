`Biomod.Models` <-
function(Model, Ids, PA.samp, TypeGLM, Test, No.trees, CV.tree, CV.ann, quant, NbRunEval, Spline, 
         DataSplit, Yweights, Roc, Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, VarImport)    
{

    #models messages in R's console
    if(Model == 'GLM') {
        if(TypeGLM == "poly"){Type <- "GLMp" ; cat("Model=GLM polynomial + quadratic")}
        if(TypeGLM == "quad"){Type <- "GLMq" ; cat("Model=GLM quadratic \n")}
        if(TypeGLM == "simple"){Type <- "GLMs" ; cat("Model=GLM linear \n")}
        if(Test == "AIC"){criteria <- 2 ; cat("\t Stepwise procedure using AIC criteria \n")}
        if(Test == "BIC"){criteria <- log(length(PA.samp)) ; cat("\t Stepwise procedure using BIC criteria \n")}
    }
     
    if(Model == 'GBM') { cat("Model=Generalised Boosting Regression \n") ; cat("\t", No.trees, "maximum different trees and lambda Fold Cross-Validation \n") }    
    if(Model == 'GAM') { cat("Model=GAM spline \n") ; cat("\t", Spline, " Degrees of smoothing \n") } 
    if(Model == 'CTA') { cat("Model=Classification tree \n") ; cat("\t", CV.tree, "Fold Cross-Validation \n") }    
    if(Model == 'ANN') { cat("Model=Artificial Neural Network \n") ; cat("\t", CV.ann, "Fold Cross Validation + 3 Repetitions \n") ; cat("Calibration and evaluation phase: Nb of cross-validations: ", ncol(Ids), "\n") }
    if(Model == 'SRE') cat("Model=Surface Range Envelop \n")  
    if(Model == 'FDA') cat("Model=Flexible Discriminant Analysis \n")
    if(Model == 'MARS') cat("Model=Multiple Adaptive Regression Splines \n")
    if(Model == 'RF') cat("Model=Breiman and Cutler's random forests for classification and regression \n")
    
    #setting some specific parameters
    if(Model == 'GLM' | Model == 'GAM') Prev <- sum(DataBIOMOD[,i])/nrow(DataBIOMOD)
    if(Model == 'CTA') set.seed(123)
    if(Model == 'RF') set.seed(71)                    
    
    #objects for storing results of the evaluation runs
    if(Roc & Model != 'SRE') AUC.train <- 0
    if(Kappa) Kappa.train <- 0
    if(TSS) TSS.train <- 0

	  SpNames <- Biomod.material$species.names
	  NbVar <- Biomod.material$NbVar
	  ErrorCounter <- 0 # Variable counting the number of models that fail in the "repetitions loop"
    if(Model == 'GBM') GBM.list <- list()	  
    RunWeights <- Yweights[,i]                                                   #for weights recalculation setting prevalence to 0.5


    #################################
	  #model loop for repetitions + final model
	  
    for(k in 1:(ncol(Ids)+1)){
        
        if(exists("model.sp")) rm(model.sp)
        if(exists("g.pred")) rm(g.pred)
        
        if(k == (ncol(Ids)+1)) { #if final model 
              if(Biomod.material$NbRepPA==0) nam <- "full"  else nam <- paste("PA", pa, sep="")
              calib.lines <- pred.lines <- PA.samp 
         } else {
              if(Biomod.material$NbRepPA==0) nam <- paste("full_rep", k, sep="")
              else  nam <- paste("PA", pa, "_rep", k, sep="")
              calib.lines <- PA.samp[Ids[,k]]
              pred.lines <- PA.samp[-Ids[,k]]
        }
        
    
        #recalculating weights to set prevalence of 0.5 ONLY if Yweights was given by user -> multiply absences weights by ratio
        #Has to be done after Ids is set because it defines which data is used for calibration
        #if Yweights NULL  ->  Yweights[whatever.lines, i] = NULL  too
        if(isnullYweights & Biomod.material$NbRepPA==0){} else{
            PW <- sum(RunWeights[calib.lines[DataBIOMOD[calib.lines, NbVar+i]==1]])  ;  assign("PW", PW, pos=1)
            AW <- sum(RunWeights[calib.lines[DataBIOMOD[calib.lines, NbVar+i]==0]])  ;  assign("AW", AW, pos=1)
            RunWeights[calib.lines[DataBIOMOD[calib.lines, NbVar+i]==0]] <- RunWeights[calib.lines[DataBIOMOD[calib.lines, NbVar+i]==0]] * (PW/AW)     #PW/AW = ratio between presWeights and absWeights
        }
        
      
        
        #necessary for the GAM, or else it doesn't read DataBIOMOD[calib.lines,] (not in the same environments)
        assign("calib.lines", calib.lines, pos=1)
        
        #building each model and making the full prediction
        #Using New Cvnnet
        if(Model == 'ANN'){
            set.seed(555)
            CV_nnet = CV.nnet(Input= DataBIOMOD[calib.lines, 1:NbVar], Target= DataBIOMOD[calib.lines, NbVar+i], nbCV=CV.ann, W=Yweights[calib.lines,i])             
            if(k==(ncol(Ids)+1)) model.sp <- nnet(DataBIOMOD[calib.lines, 1:NbVar], DataBIOMOD[calib.lines, NbVar+i], size=CV_nnet[1,1], rang=0.1, decay=CV_nnet[1,2], maxit=200, trace=F)
                        else try(model.sp <- nnet(DataBIOMOD[calib.lines, 1:NbVar], DataBIOMOD[calib.lines, NbVar+i], size=CV_nnet[1,1], rang=0.1, decay=CV_nnet[1,2], maxit=200, trace=F), silent=T) 
  
            if(exists("model.sp")) TempArray <- predict(model.sp, DataBIOMOD[PA.samp,], type="raw")
        }####     
        if(Model == 'CTA'){
            temp <- rpart.control(xval=CV.tree, minbucket=5, minsplit=5,cp=0.001, maxdepth=25)
            if(k==(ncol(Ids)+1)) model.sp <- rpart(eval(parse(text=paste("as.factor(", SpNames[i], ")" ,paste(scopeExpSyst(DataBIOMOD[1:10,1:NbVar], "CTA"),collapse="")))),DataBIOMOD[calib.lines,], weights=Yweights[calib.lines,i], control=temp)
                        else try(model.sp <- rpart(eval(parse(text=paste("as.factor(", SpNames[i], ")", paste(scopeExpSyst(DataBIOMOD[1:10,1:NbVar], "CTA"),collapse="")))),DataBIOMOD[calib.lines,], weights=Yweights[calib.lines,i], control=temp), silent=T)
            if(exists("model.sp")){
                tr <- as.data.frame(model.sp$cptable)
                tr[,6] <- tr[,4] + tr[,5]
                tr <- tr[tr[,2]!=0,]
                Cp <- tr[tr[,6] == min(tr[,6]), 1]
                if(length(Cp) ==1) model.sp <- prune(model.sp, cp=Cp)
                else model.sp <- prune(model.sp, cp=Cp[2])
                g.pred <- data.frame(as.integer(as.numeric(predict(model.sp, DataBIOMOD[PA.samp,1:NbVar], type="prob")[,2]) *1000))
            }
        }####
        if(Model == 'GAM'){
            gamStart <- gam(eval(parse(text=paste(paste(SpNames[i]),"~1",collapse=""))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
            if(k==(ncol(Ids)+1)) model.sp <- step.gam(gamStart, scope(DataBIOMOD[1:10, 1:NbVar],"s", Spline), keep=functionkeep, direction="both", trace=F, control=gam.control(maxit=50, bf.maxit=50))
                        else try(model.sp <- step.gam(gamStart, scope(DataBIOMOD[1:10, 1:NbVar],"s", Spline), keep=functionkeep, direction="both", trace=F, control=gam.control(maxit=50, bf.maxit=50)), silent=T)
            
            #New model.sp <- gam(eval(parse(text=paste(SpNames[i],paste(scopeGAM(DataBIOMOD[1:10,1:NbVar]),collapse="")))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
            if(exists("model.sp")) g.pred <- data.frame(as.integer(as.numeric(testnull(model.sp, Prev, DataBIOMOD[PA.samp,])) *1000))
        }####   
        if(Model == 'GBM'){
            if(k==(ncol(Ids)+1)) model.sp <- gbm(eval(parse(text=paste(SpNames[i],paste(scopeExpSyst(DataBIOMOD[1:10,1:NbVar], "GBM"),collapse="")))), data=DataBIOMOD[calib.lines,], distribution="bernoulli", var.monotone=rep(0, length=NbVar), w=Yweights[calib.lines,i], interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, n.trees=No.trees, verbose=F, cv.folds=5)        
                        else try(model.sp <- gbm(eval(parse(text=paste(SpNames[i],paste(scopeExpSyst(DataBIOMOD[1:10,1:NbVar], "GBM"),collapse="")))), data=DataBIOMOD[calib.lines,], distribution="bernoulli", var.monotone=rep(0, length=NbVar), w=Yweights[calib.lines,i], interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, n.trees=No.trees, verbose=F, cv.folds=5), silent=T)
            
            if(exists("model.sp")){
                best.iter <- gbm.perf(model.sp, method="cv", plot.it=F)
                g.pred <- data.frame(as.integer(predict.gbm(model.sp, DataBIOMOD[PA.samp,], best.iter, type="response")*1000))
            }
        }####        
        if(Model == 'GLM'){
            glmStart <- glm(eval(parse(text=paste(paste(SpNames[i]), "~1", collapse=""))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
            if(k==(ncol(Ids)+1)) model.sp <- stepAIC(glmStart, scopeExpSyst(DataBIOMOD[1:10, 1:NbVar], Type), direction="both", trace=F, control=glm.control(maxit=100), k=criteria)
                        else try(model.sp <- stepAIC(glmStart, scopeExpSyst(DataBIOMOD[1:10, 1:NbVar], Type), direction="both", trace=F, control=glm.control(maxit=100), k=criteria), silent=T)
            if(exists("model.sp")) g.pred <- data.frame(as.integer(testnull(model.sp, Prev, DataBIOMOD[PA.samp,]) *1000)) 
        }####                     
        if(Model == 'MARS'){
            if(k==(ncol(Ids)+1)){
                if(is.null(Yweights)) model.sp <- mars(x=DataBIOMOD[calib.lines, 1:NbVar], y=DataBIOMOD[calib.lines, NbVar+i], degree=2)
                                 else model.sp <- mars(x=DataBIOMOD[calib.lines, 1:NbVar], y=DataBIOMOD[calib.lines, NbVar+i], degree=2, w=Yweights[calib.lines,i])
            } else{
                 if(is.null(Yweights)) try(model.sp <- mars(x=DataBIOMOD[calib.lines, 1:NbVar], y=DataBIOMOD[calib.lines, NbVar+i], degree=2), silent=T)
                                 else try(model.sp <- mars(x=DataBIOMOD[calib.lines, 1:NbVar], y=DataBIOMOD[calib.lines, NbVar+i], degree=2, w=Yweights[calib.lines,i]), silent=T)
            }
            if(exists("model.sp")) TempArray <- predict(model.sp, DataBIOMOD[PA.samp,1:NbVar])
        }#### 
        if(Model == 'FDA') {
            if(k==(ncol(Ids)+1)) model.sp <- fda(eval(parse(text=paste(SpNames[i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:NbVar], "FDA"), collapse="")))), data=DataBIOMOD[calib.lines,], method=mars)
                        else try(model.sp <- fda(eval(parse(text=paste(SpNames[i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:NbVar], "FDA"), collapse="")))), data=DataBIOMOD[calib.lines,], method=mars), silent=T)
            if(exists("model.sp")) TempArray <- predict(model.sp, DataBIOMOD[PA.samp,1:NbVar], type="post")[,2]
        }#### 
        if(Model == 'RF') {
            if(k==(ncol(Ids)+1)) model.sp <- randomForest(x=DataBIOMOD[calib.lines, 1:NbVar], y=as.factor(DataBIOMOD[calib.lines, NbVar+i]), ntree=750, mtry=NbVar/2, importance=TRUE)
                        else try(model.sp <- randomForest(x=DataBIOMOD[calib.lines, 1:NbVar], y=as.factor(DataBIOMOD[calib.lines, NbVar+i]), ntree=750, mtry=NbVar/2, importance=TRUE), silent=T)
            if(exists("model.sp")) g.pred <- data.frame(as.integer(predict(model.sp, DataBIOMOD[PA.samp,1:NbVar], type="prob")[,2] *1000))
        }
        if(Model == 'SRE') g.pred <- data.frame(as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD[calib.lines,]$",paste(SpNames[i]), collapse=""))), DataBIOMOD[calib.lines,1:NbVar],DataBIOMOD[PA.samp,], quant)) *1000))

        
        
        #Building a Rescaling glm
        if(any(c("ANN", "FDA", "MARS")== Model)) g.pred <- data.frame(as.integer(Rescaler4(as.numeric(TempArray), ref=DataBIOMOD[PA.samp,NbVar+i], run=paste(SpNames[i], "_", Model ,"_", nam, sep=""), original=T) *1000))
        
        #compute the evaluation stats if a repetition run
        if(!exists("g.pred")){
            ErrorCounter <- ErrorCounter+1 # If the current evaluation run failed, we record this failure by increasing the counter by 1
            BM[["calibration.failures"]] <- c(BM[["calibration.failures"]], paste(SpNames[i], "_", Model, "_", nam, sep=""))
        }
        if(exists("g.pred")){  # The evaluation can only run if the current evaluation run did not fail 
        
            if(k != (ncol(Ids)+1)){
               if(Roc) auc.stat <-  somers2(g.pred[-Ids[,k],], DataBIOMOD[pred.lines,NbVar+i])["C"]
               if(Kappa) kappa.stat <- KappaRepet(DataBIOMOD[pred.lines,NbVar+i], g.pred[-Ids[,k],])$Kappa
               if(TSS) tss.stat <- KappaRepet(DataBIOMOD[pred.lines,NbVar+i], g.pred[-Ids[,k],], TSS=T)$TSS
            }
            
            #running the evaluation procedures for the evaluation runs
            #no extra prediction to be made -> using the g.pred (on full PA.samp) and getting the right lines
            if(k != (ncol(Ids)+1)){ #if repetition run           
                if(Model == 'SRE'){ 
                    if(Kappa) Kappa.train <- Kappa.train + KappaSRE(DataBIOMOD[pred.lines, NbVar+i], g.pred[-Ids[,k],])$Kappa 
                    if(TSS) TSS.train <- TSS.train + KappaSRE(DataBIOMOD[pred.lines, NbVar+i], g.pred[-Ids[,k],], TSS=T)$TSS 
                } else{
                    if(Roc) AUC.train <- AUC.train + auc.stat
                    if(Kappa) Kappa.train <- Kappa.train + kappa.stat
                    if(TSS) TSS.train <- TSS.train + tss.stat
                }
            } else {
                  if(Model == 'SRE'){ 
                    if(Kappa) Kappa.final <- KappaSRE(DataBIOMOD[pred.lines, NbVar+i], g.pred[,])$Kappa 
                    if(TSS) TSS.final <- KappaSRE(DataBIOMOD[pred.lines, NbVar+i], g.pred[,], TSS=T)$TSS 
                } else{
                    if(Roc) AUC.final <- somers2(g.pred[,], DataBIOMOD[pred.lines, NbVar+i])["C"]
                    if(Kappa) Kappa.final <- KappaRepet(DataBIOMOD[pred.lines, NbVar+i], g.pred[,])$Kappa
                    if(TSS) TSS.final <- KappaRepet(DataBIOMOD[pred.lines, NbVar+i], g.pred[,], TSS=T)$TSS
                }
            }
            
            
            #saving the evaluation stats for the repetition models  
            if(k != (ncol(Ids)+1)){
                if(Roc && Model != 'SRE') Evaluation.results.Roc[[paste(SpNames[i], "_", nam, sep="")]][Model,] <- c(round(auc.stat, digits=3), 'none', round(somers2(g.pred[,], DataBIOMOD[PA.samp, NbVar+i])["C"], digits=3), round(CutOff.Optimised(DataBIOMOD[PA.samp, NbVar+i], g.pred[,]), digits=3))  
                if(Kappa) Evaluation.results.Kappa[[paste(SpNames[i], "_", nam, sep="")]][Model,] <- c(round(kappa.stat,digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,NbVar+i], g.pred[,])[c(1,2,4,6)]) 
                if(TSS) Evaluation.results.TSS[[paste(SpNames[i], "_", nam, sep="")]][Model,] <- c(round(tss.stat,digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,NbVar+i], g.pred[,], TSS=T)[c(1,2,4,6)])   
            }
            
            #saving best.iter for GBM
            if(Model == 'GBM')  eval(parse(text=paste("GBM.list$", nam," <- best.iter", sep="")))
            
            #save the predictions in the array
            if(k == (ncol(Ids)+1)) Array[,Model,1,pa] <- g.pred[,]  else  Array[,Model,(k+1),pa] <- g.pred[,] 
            
            #saving the model on the hard disk
            if(Model != 'SRE') eval(parse(text=paste("assign('",SpNames[i],"_", Model, "_", nam,"', model.sp)", sep="")))
            if(Model != 'SRE') eval(parse(text=paste("save(",SpNames[i],"_", Model, "_", nam, ",file='", getwd(), "/models/", SpNames[i],"_", Model, "_", nam,"')", sep="")))   
        
        } #if model did not fail
    } #nbruns k loop
        
    
    
    assign("Array", Array, pos=1)
    assign("BM", BM, pos=1)
    if(Model == 'GBM') assign("GBM.list", GBM.list, pos=1)
    
    #mean evaluations from the calibration
    if(Roc & Model != "SRE") AUC.train <- AUC.train/(ncol(Ids)-ErrorCounter)
    if(Kappa) Kappa.train <- Kappa.train/(ncol(Ids)-ErrorCounter)
    if(TSS) TSS.train <- TSS.train/(ncol(Ids)-ErrorCounter)
    
    
    #Evaluation of Predictor Importance in the model:
    if(VarImport > 0){
        cat("Evaluating Predictor Contributions in " , Model, "...", "\n")
        TempVarImp <- as.data.frame(matrix(data=0, nrow=1, ncol=Biomod.material$NbVar))
        names(TempVarImp) <- names(DataBIOMOD)[1:Biomod.material$NbVar]
        
        for(J in 1:Biomod.material$NbVar){
            for(K in 1:VarImport){
                TempDS <- DataBIOMOD[PA.samp,1:Biomod.material$NbVar]
                TempDS[,J] <- sample(TempDS[,J])
                 
                if(Model == 'ANN') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler4(as.numeric(predict(model.sp, TempDS, type="raw")), ref=DataBIOMOD[PA.samp,NbVar+i], run=paste(SpNames[i], "_ANN_", nam, sep="")) *1000))
                if(Model == 'CTA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(predict(model.sp, TempDS, type="prob")[,2]) *1000))
                if(Model == 'GLM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GAM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GBM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(predict.gbm(model.sp, TempDS, best.iter, type='response')) *1000))
                if(Model == 'MARS')TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler4(predict(model.sp, TempDS[,1:NbVar]),        ref=DataBIOMOD[PA.samp,NbVar+i], run=paste(SpNames[i], "_MARS_", nam, sep="")) *1000))
                if(Model == 'FDA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler4(predict(model.sp, TempDS, type="post")[,2], ref=DataBIOMOD[PA.samp,NbVar+i], run=paste(SpNames[i], "_FDA_", nam, sep="")) *1000))
                if(Model == 'RF')  TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(predict(model.sp, TempDS, type="prob")[,2]) *1000))
                if(Model == 'SRE') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD[PA.samp,]$",paste(SpNames[i]), collapse=""))), TempDS,DataBIOMOD[PA.samp,], quant)) *1000))
           
           }
        }
        VarImportance[[i]][Model,] <-  round(1 - (TempVarImp/VarImport), digits=3)
    }
	  assign("VarImportance", VarImportance, pos=1)
	  
	 
    #Predictions on the independent data if any.
    if(exists("DataEvalBIOMOD")){
        if(Model == 'GLM') predind <- testnull(model.sp, Prev, DataEvalBIOMOD)
        if(Model == 'GAM') predind <- testnull(model.sp, Prev, DataEvalBIOMOD)
        if(Model == 'ANN') predind <- predict(model.sp, DataEvalBIOMOD, type="raw")
        if(Model == 'CTA') predind <- predict(model.sp, DataEvalBIOMOD[,1:NbVar], type="prob")[,2]
        if(Model == 'GBM') predind <- predict.gbm(model.sp, DataEvalBIOMOD, best.iter, type='response')        
        if(Model == 'MARS') predind <- predict(model.sp, DataEvalBIOMOD[,1:NbVar])
        if(Model == 'FDA') predind <- predict(model.sp, DataEvalBIOMOD[,1:NbVar], type="post")[,2]
        if(Model == 'RF') predind <- predict(model.sp, DataEvalBIOMOD[,1:NbVar], type="prob")[,2]
        if(Model == 'SRE') predind <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$", paste(SpNames[i]), collapse=""))), DataBIOMOD[1:NbVar], DataEvalBIOMOD, quant)) *1000)
    
        if(any(c("ANN", "FDA", "MARS")== Model)) predind <- Rescaler4(predind, ref=DataBIOMOD[PA.samp,NbVar+i], run=paste(SpNames[i], "_",Model ,"_", nam, sep=""))
        predind <- as.integer(as.numeric(predind) *1000)
    }



    #running the evaluation procedures
    if(Roc && Model != 'SRE'){
        Evaluation.results.Roc[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,] <- c(round(AUC.train, digits=3), 'none', round(somers2(g.pred[,], DataBIOMOD[PA.samp,NbVar+i])["C"], digits=3), round(CutOff.Optimised(DataBIOMOD[PA.samp,NbVar+i], g.pred[,]), digits=3))  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.Roc[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,2] <- round(somers2(predind,DataEvalBIOMOD[,NbVar +i])["C"], digits=3) 
        assign("Evaluation.results.Roc", Evaluation.results.Roc, pos=1)
	  }    
    if(Kappa){
	      Evaluation.results.Kappa[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,] <- c(round(Kappa.train, digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,NbVar+i], g.pred[,])[c(1,2,4,6)]) 
        if(exists("DataEvalBIOMOD")) Evaluation.results.Kappa[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,NbVar+i], predind)$Kappa, digits=3) 
        assign("Evaluation.results.Kappa", Evaluation.results.Kappa, pos=1)
    }	
    if(TSS){
	      Evaluation.results.TSS[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,] <- c(round(TSS.train, digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,NbVar+i], g.pred[,], TSS=T)[c(1,2,4,6)])  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.TSS[[paste(Biomod.material$species.names[i], "_", nam, sep="")]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,NbVar+i], predind, TSS=T)$TSS, digits=3) 
        assign("Evaluation.results.TSS", Evaluation.results.TSS, pos=1)
    }
    	
  	if(exists("DataEvalBIOMOD") && KeepPredIndependent) assign("predind", predind, pos=1)
  	assign("Biomod.material", Biomod.material, pos=1)
}
