`Biomod.Models` <-
function(Model, Ids, PA.samp, TypeGLM, Test, No.trees, CV.tree, CV.ann, Perc025, Perc05, NbRunEval, Spline, 
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
    if(Model == 'GAM') { cat("Model=GAM spline \n") ; cat("\t lambda Degrees of smoothing \n") } 
    if(Model == 'CTA') { cat("Model=Classification tree \n") ; cat("\t", CV.tree, "Fold Cross-Validation \n") }    
    if(Model == 'ANN') { cat("Model=Artificial Neural Network \n") ; cat("\t", CV.ann, "Fold Cross Validation + 3 Repetitions \n") ; cat("Calibration and evaluation phase: Nb of cross-validations: ", ncol(Ids), "\n") }
    if(Model == 'SRE') cat("Model=Surface Range Envelop \n")  
    if(Model == 'MDA') { cat("Model=Mixture Discriminant Analysis \n") ; library(reshape) }
    if(Model == 'MARS') cat("Model=Multiple Adaptive Regression Splines \n")
    if(Model == 'RF') cat("Model=Breiman and Cutler's random forests for classification and regression \n")

    #setting some specific parameters
    if(Model == 'GLM' | Model == 'GAM') Prev <- sum(DataBIOMOD[,i])/nrow(DataBIOMOD)
    if(Model == 'CTA') set.seed(123)
    if(Model == 'RF') set.seed(71)
    if(Model == 'ANN'){ decay.final <- 0 ; size.final <- 0 }                     
    
    #objects for storing results of the evaluation runs
    if(Model != 'SRE') AUC.train <- 0
    Kappa.train <- 0
    TSS.train <- 0
    #list for storing some model outputs (ex: best.iter for gbm)
    g.list <- list()
        
	
	
	  #################################
	  #################################
	
	  #model loop for repetitions + final model
	  
    for(k in 1:(ncol(Ids)+1)){
  	    
  	    if(k == (ncol(Ids)+1)) { 
              if(Biomod.material$NbRepPA==0) m.name <- "full"
              else m.name <- paste("PA", pa, sep="")
              calib.lines <- pred.lines <- PA.samp 
         } else {
              if(Biomod.material$NbRepPA==0) m.name <- paste("full_rep", k, sep="")
              else  m.name <- paste("PA", pa, "_rep", k, sep="")  
              calib.lines <- PA.samp[Ids[,k]]
              pred.lines <- PA.samp[-Ids[,k]]
        } 
  	
  	    assign("calib.lines", calib.lines, pos=1)
  	
  	    #building each model and making the full prediction
        if(Model == 'ANN'){
            
            if(k==(ncol(Ids)+1) && ncol(Ids)!=0){ #final model with prior evaluation runs (NbRunEval!=0)
            
                decay.final <- decay.final/ncol(Ids)
                size.final <- size.final/ncol(Ids) 
            		model.sp <- nnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"),collapse="")))), 
            		             data=DataBIOMOD[calib.lines,], weights=Yweights[calib.lines,i], size=size.final, linout=F, entropy=T, skip=F, decay=decay.final, maxit=100, trace=F)

            } else { #evaluation runs, or final model if NbRunEval==0
                truth <- DataBIOMOD[calib.lines,Biomod.material[["NbVar"]]+i]
                set.seed(200)
                tr <- CVnnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"), collapse="")))), 
                              data=DataBIOMOD[calib.lines,], truth=truth, linout= F, entropy=T, skip=F, maxit=100, nifold=CV.ann)
                decay <- tr[tr[,3] == max(tr[,3]),2]
                size <- tr[tr[,3] == max(tr[,3]),1]
                model.sp <- nnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"), collapse="")))), 
                             data=DataBIOMOD[calib.lines,], weights=Yweights[calib.lines,i], size=size, linout=F, entropy=T, skip=F, decay=decay, maxit=100,trace=F) 
                             
               decay.final <- decay.final + decay
		           size.final <- size.final + size
            } 
           
            TempArray <- predict(model.sp, DataBIOMOD[PA.samp,], type="raw")
      	    g.pred <- data.frame(as.integer(Rescaler2(as.numeric(TempArray), type="range") *1000)) 		
  	    }
        if(Model == 'CTA'){
            temp <- rpart.control(xval=CV.tree, minbucket=5, minsplit=5,cp=0.001, maxdepth=25)
            model.sp <- rpart(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeExpSyst(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]], "CTA"),collapse="")))),DataBIOMOD[calib.lines,], weights=Yweights[calib.lines,i], control=temp)
            tr <- as.data.frame(model.sp$cptable)
            tr[,6] <- tr[,4] + tr[,5]
            tr <- tr[tr[,2]!=0,]
            Cp <- tr[tr[,6] == min(tr[,6]), 1]
            if(length(Cp) ==1) model.sp <- prune(model.sp, cp=Cp)
            else model.sp <- prune(model.sp, cp=Cp[2])
        		g.pred <- data.frame(as.integer(as.numeric(predict(model.sp, DataBIOMOD[PA.samp,], type="vector")) *1000))
        }
        if(Model == 'GAM'){
        		gamStart <- gam(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]),"~1",collapse=""))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
          	model.sp <- step.gam(gamStart, scope(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]],"s", Spline), keep=functionkeep, direction="both", trace=F, control=gam.control(maxit=50, bf.maxit=50))	
      #New  model.sp <- gam(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeGAM(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]]),collapse="")))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
            g.pred <- data.frame(as.integer(as.numeric(testnull(model.sp, Prev, DataBIOMOD[PA.samp,])) *1000))
        }   
        if(Model == 'GBM'){
            model.sp <- gbm.fit(x=DataBIOMOD[calib.lines,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[calib.lines,Biomod.material[["NbVar"]]+i], distribution="bernoulli", w=Yweights[calib.lines,i], var.monotone=rep(0, length=Biomod.material[["NbVar"]]), n.trees=No.trees, 
                         interaction.depth=3, shrinkage=0.01, bag.fraction=0.5, train.fraction=1, verbose=F,  var.names=colnames(DataBIOMOD[1:Biomod.material[["NbVar"]]]),
                         response.name = Biomod.material[["species.names"]][i])
            best.iter <- gbm.perf(model.sp,method="OOB",plot.it=F)
        		g.pred <- data.frame(as.integer(predict.gbm(model.sp, DataBIOMOD[PA.samp,], best.iter, type="response")*1000))
        }
        if(Model == 'GLM'){
            glmStart <- glm(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]), "~1", collapse=""))), data=DataBIOMOD[calib.lines,], family=binomial, weights=Yweights[calib.lines,i])
            model.sp <- stepAIC(glmStart, scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], Type), direction="both", trace=F, control=glm.control(maxit=100), k=criteria)
        		g.pred <- data.frame(as.integer(testnull(model.sp, Prev, DataBIOMOD[PA.samp,]) *1000)) 
        }
        if(Model == 'MARS'){
            if(is.null(Yweights)) model.sp <- mars(x=DataBIOMOD[calib.lines,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[calib.lines,Biomod.material[["NbVar"]]+i], degree=2)
                             else model.sp <- mars(x=DataBIOMOD[calib.lines,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[calib.lines,Biomod.material[["NbVar"]]+i], degree=2, w=Yweights[calib.lines,i])
                             
        		TempArray <- predict(model.sp, DataBIOMOD[PA.samp,1:Biomod.material[["NbVar"]]])
        		g.pred <- data.frame(as.integer(Rescaler2(TempArray, type="range") *1000))
        }
        if(Model == 'MDA') {
            model.sp <- mda(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "MDA"), collapse="")))), data=DataBIOMOD[calib.lines,], method=mars)
            TempArray <- predict(model.sp, DataBIOMOD[PA.samp,1:Biomod.material[["NbVar"]]], type="post")[,2]
        		g.pred <- data.frame(as.integer(Rescaler2(TempArray, type="range") *1000))
        }

        if(Model == 'RF') {
            model.sp <- randomForest(x=DataBIOMOD[calib.lines, 1:Biomod.material[["NbVar"]]], y=as.factor(DataBIOMOD[calib.lines, Biomod.material[["NbVar"]]+i]), ntree=750, mtry=Biomod.material[["NbVar"]]/2, importance=TRUE)
            TempArray <- predict(model.sp, DataBIOMOD[PA.samp,1:Biomod.material[["NbVar"]]], type="prob")[,2]
        		g.pred <- data.frame(as.integer(Rescaler2(TempArray, type="range") *1000))
        }
        if(Model == 'SRE') g.pred <- data.frame(as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD[calib.lines,]$",paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[calib.lines,1:Biomod.material[["NbVar"]]],DataBIOMOD[PA.samp,], Perc025, Perc05)) *1000))
        
        #running the evaluation procedures for the evaluation runs
        #no extra prediction to be made -> using the g.pred (on full PA.samp) and getting the right lines
    		if(k!=(ncol(Ids)+1)){            
            if(Model == 'SRE'){ 
                Kappa.train <- Kappa.train + KappaSRE(DataBIOMOD[pred.lines,Biomod.material[["NbVar"]]+i], g.pred[-Ids[,k],])$Kappa 
                TSS.train <- TSS.train + KappaSRE(DataBIOMOD[pred.lines,Biomod.material[["NbVar"]]+i], g.pred[-Ids[,k],], TSS=T)$TSS 
            } else{
                AUC.train <- AUC.train + somers2(g.pred[-Ids[,k],], DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i])["C"]
          			Kappa.train <- Kappa.train + KappaRepet(DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i], g.pred[-Ids[,k],])$Kappa
          			TSS.train <- TSS.train + KappaRepet(DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i], g.pred[-Ids[,k],], TSS=T)$TSS
      		  }
        } else {
              if(Model == 'SRE'){ 
                Kappa.final <- KappaSRE(DataBIOMOD[pred.lines,Biomod.material[["NbVar"]]+i], g.pred[,])$Kappa 
                TSS.final <- KappaSRE(DataBIOMOD[pred.lines,Biomod.material[["NbVar"]]+i], g.pred[,], TSS=T)$TSS 
            } else{
                AUC.final <- somers2(g.pred[,], DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i])["C"]
          			Kappa.final <- KappaRepet(DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i], g.pred[,])$Kappa
          			TSS.final <- KappaRepet(DataBIOMOD[pred.lines, Biomod.material[["NbVar"]]+i], g.pred[,], TSS=T)$TSS
      		  }
        }
        
        #save the prediction in array
        if(Biomod.material[["NbRepPA"]] == 0){
            if(k == (ncol(Ids)+1)) ARRAY[,Model,1] <- g.pred[,]  else  ARRAY[,Model,(k+1)] <- g.pred[,]
        } else { 
            if(k == (ncol(Ids)+1)) ARRAY[,Model,1,pa] <- g.pred[,]  else  ARRAY[,Model,(k+1),pa] <- g.pred[,]
        }
        
        #saving the model on the hard disk
        if(Model != 'SRE') eval(parse(text=paste("assign('",Biomod.material[["species.names"]][i],"_", Model, "_", m.name,"', model.sp)", sep="")))
  	    if(Model != 'SRE') eval(parse(text=paste("save(",Biomod.material[["species.names"]][i],"_", Model, "_", m.name, ",file='", getwd(), "/models/", Biomod.material[["species.names"]][i],"_", Model, "_", m.name,"')", sep="")))
    
        # saving best.iter for GBM ; Save Minumum and Maxium of the calibration prediction range for rescaling steps:
        if(Model == 'GBM')  eval(parse(text=paste("g.list$best.iter$", m.name," <- best.iter", sep="")))
        if(Model == 'ANN' | Model == 'MDA' | Model == 'MARS' | Model == 'RF') eval(parse(text=paste("g.list$RawPred$", m.name," <- range(TempArray)", sep="")))       
    
    }
	  
    assign("ARRAY", ARRAY, pos=1)         
	  
    #mean evaluations from the calibration
    if(Model != 'SRE') AUC.train <- AUC.train/ncol(Ids)
    Kappa.train <- Kappa.train/ncol(Ids)
    TSS.train <- TSS.train/ncol(Ids)
      
    
    #Evaluation of Predictor Importance in the model:
    if(VarImport > 0){
        cat("Evaluating Predictor Contributions in " , Model, "...", "\n")
        TempVarImp <- as.data.frame(matrix(data=0, nrow=1, ncol=Biomod.material[["NbVar"]]))
        names(TempVarImp) <- names(DataBIOMOD)[1:Biomod.material[["NbVar"]]]
        
        for(J in 1:Biomod.material[["NbVar"]]){
            for(K in 1:VarImport){
                TempDS <- DataBIOMOD[PA.samp,1:Biomod.material[["NbVar"]]]
                TempDS[,J] <- sample(TempDS[,J])
                 
                if(Model == 'ANN') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler2(as.numeric(predict(model.sp, TempDS, type="raw")), type="range", OriMinMax=range(TempArray)) *1000))
                if(Model == 'CTA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(predict(model.sp, TempDS, type="vector")) *1000))
                if(Model == 'GLM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GAM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GBM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(predict.gbm(model.sp, TempDS, best.iter, type='response')) *1000))
                if(Model == 'MARS') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler2(predict(model.sp, TempDS[,1:Biomod.material[["NbVar"]]]), type="range", OriMinMax=range(TempArray)) *1000))
                if(Model == 'MDA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler2(predict(model.sp, TempDS, type="post")[,2], type="range", OriMinMax=range(TempArray)) *1000))
                if(Model == 'RF') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(Rescaler2(predict(model.sp, TempDS, type="prob")[,2], type="range", OriMinMax=range(TempArray)) *1000))
                if(Model == 'SRE') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred[,], as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD[PA.samp,]$",paste(Biomod.material[["species.names"]][i]), collapse=""))), TempDS,DataBIOMOD[PA.samp,], Perc025, Perc05)) *1000))    
           }
        }
        VarImportance[[i]][Model,] <-  round(1 - (TempVarImp/VarImport), digits=3)
    }
	 
	 
    #Predictions on the independent data if any.
    if(exists("DataEvalBIOMOD")){
        if(Model == 'GLM') predind <- as.integer(as.numeric(testnull(model.sp, Prev, DataEvalBIOMOD)) *1000)
        if(Model == 'GAM') predind <- as.integer(as.numeric(testnull(model.sp, Prev, DataEvalBIOMOD)) *1000)
        if(Model == 'ANN') predind <- as.integer(Rescaler2(as.numeric(predict(model.sp, DataEvalBIOMOD, type="raw")), type="range", OriMinMax=range(TempArray)) *1000)
        if(Model == 'CTA') predind <- as.integer(as.numeric(predict(model.sp, DataEvalBIOMOD, type="vector")) *1000)
        if(Model == 'GBM') predind <- as.integer(predict.gbm(model.sp, DataEvalBIOMOD, best.iter, type='response') *1000)        
        if(Model == 'MARS') predind <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]]), type="range", OriMinMax=range(TempArray)) *1000)
        if(Model == 'MDA') predind <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]], type="post")[,2], type="range", OriMinMax=range(TempArray)) *1000)
        if(Model == 'RF') predind <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]], type="prob")[,2], type="range", OriMinMax=range(TempArray)) *1000)
        if(Model == 'SRE') predind <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$", paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[1:Biomod.material[["NbVar"]]], DataEvalBIOMOD, Perc025, Perc05)) *1000)
    }

    #running the evaluation procedures
    if(Roc && Model != 'SRE'){
        Evaluation.results.Roc[[i]][Model,] <- c(round(AUC.train, digits=3), 'none', round(somers2(g.pred[,], DataBIOMOD[PA.samp,Biomod.material[["NbVar"]]+i])["C"], digits=3), round(CutOff.Optimised(DataBIOMOD[PA.samp,Biomod.material[["NbVar"]]+i], g.pred[,]), digits=3))  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.Roc[[i]][Model,2] <- round(somers2(predind,DataEvalBIOMOD[,Biomod.material[["NbVar"]] +i])["C"], digits=3) 
        assign("Evaluation.results.Roc", Evaluation.results.Roc, pos=1)
	  }    
    if(Kappa){
	      Evaluation.results.Kappa[[i]][Model,] <- c(round(Kappa.train, digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,Biomod.material[["NbVar"]]+i], g.pred[,], TSS=T)[c(1,2,4,6)]) 
        if(exists("DataEvalBIOMOD")) Evaluation.results.Kappa[[i]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,Biomod.material[["NbVar"]]+i], predind)$Kappa, digits=3) 
        assign("Evaluation.results.Kappa", Evaluation.results.Kappa, pos=1)
    }	
    if(TSS){
	      Evaluation.results.TSS[[i]][Model,] <- c(round(TSS.train, digits=3), 'none', KappaRepet(DataBIOMOD[PA.samp,Biomod.material[["NbVar"]]+i], g.pred[,], TSS=T)[c(1,2,4,6)])  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.TSS[[i]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,Biomod.material[["NbVar"]]+i], predind, TSS=T)$TSS, digits=3) 
        assign("Evaluation.results.TSS", Evaluation.results.TSS, pos=1)
    }
    	
  	if(exists("DataEvalBIOMOD") && KeepPredIndependent) assign("predind", predind, pos=1)
    assign("VarImportance", VarImportance, pos=1)
                  
    return(g.list)
}

