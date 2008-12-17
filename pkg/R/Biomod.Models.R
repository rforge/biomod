`Biomod.Models` <-
function(Model, Ids, TypeGLM, Test, No.trees, CV.tree, CV.ann, Perc025, Perc05, NbRunEval, Spline, 
                            DataSplit, Yweights, Roc, Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, VarImport)
{

    if(Model == 'GLM') {
        if(TypeGLM == "poly"){Type <- "GLMp"; cat("Model=GLM polynomial + quadratic")}
        if(TypeGLM == "quad"){Type <- "GLMq"; cat("Model=GLM quadratic \n")}
        if(TypeGLM == "simple"){Type <- "GLMs"; cat("Model=GLM linear \n")}
        if(Test == "AIC"){criteria <- 2; cat("\t Stepwise procedure using AIC criteria \n")}
        if(Test == "BIC"){criteria <- log(nrow(Ids)); cat("\t Stepwise procedure using BIC criteria \n")}
    }
    if(Model == 'GBM') {
        cat("Model=Generalised Boosting Regression \n")
        cat("\t", No.trees, "maximum different trees and lambda Fold Cross-Validation \n")
    }    
    if(Model == 'GAM') {
        cat("Model=GAM spline \n")
        cat("\t lambda Degrees of smoothing \n")
    } 
    if(Model == 'CTA') { 
        cat("Model=Classification tree \n")
        cat("\t", CV.tree, "Fold Cross-Validation \n")
    }    
    if(Model == 'ANN') {
        cat("Model=Artificial Neural Network \n")
        cat("\t", CV.ann, "Fold Cross Validation + 3 Repetitions \n")
        cat("Calibration and evaluation phase: Nb of cross-validations: ", ncol(Ids), "\n")
    }
    if(Model == 'SRE') cat("Model=Surface Range Envelop \n")  
    if(Model == 'MDA') {
        cat("Model=Mixture Discriminant Analysis \n")
        library(reshape)
    }
    if(Model == 'MARS') cat("Model=Multiple Adaptive Regression Splines \n")
    if(Model == 'RF') cat("Model=Breiman and Cutler's random forests for classification and regression \n")

    
    if(Model == 'GLM' | Model == 'GAM') Prev <- sum(DataBIOMOD[,i])/nrow(DataBIOMOD)
    if(Model == 'CTA') set.seed(123)
    if(Model == 'RF') set.seed(71)
    if(Model == 'ANN'){ 
        decay.final <- 0
        size.final <- 0
    }                     
    
    Run <- ncol(Ids)
    if(Model != 'SRE') AUC.train <- 0
    Kappa.train <- 0
    TSS.train <- 0
    g.list <- list()
        
	
  	#Cross validation
  	for(k in 1:Run){
  	
  	    assign("sp", Ids[,k], pos=1)
  	
        if(Model == 'ANN'){
            truth <- DataBIOMOD[sp,Biomod.material[["NbVar"]]+i]
            set.seed(200)
            tr <- CVnnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"), collapse="")))), data=DataBIOMOD[sp,], truth=truth, linout= F, entropy=T, skip=F, maxit=100, nifold=CV.ann)
            decay <- tr[tr[,3] == max(tr[,3]),2]
            size <- tr[tr[,3] == max(tr[,3]),1]
            model.sp <- nnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"), collapse="")))), data=DataBIOMOD[sp,], weights=Yweights[sp,i], size=size,linout=F, entropy=T, skip=F, decay=decay, maxit=100,trace=F)       
        }	
        if(Model == 'CTA'){
            temp <- rpart.control(xval=CV.tree, minbucket=5, minsplit=5,cp=0.001, maxdepth=25)
            model.sp <- rpart(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeExpSyst(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]], "CTA"),collapse="")))),DataBIOMOD[sp,], weights=Yweights[sp,i], control=temp)
            tr <- as.data.frame(model.sp$cptable)
            tr[,6] <- tr[,4] + tr[,5]
            tr <- tr[tr[,2]!=0,]
            Cp <- tr[tr[,6] == min(tr[,6]), 1]
            if(length(Cp) ==1) model.sp <- prune(model.sp, cp=Cp)
            else model.sp <- prune(model.sp, cp=Cp[2])
        }	
        if(Model == 'GAM'){ 
        		gamStart <- gam(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]),"~1",collapse=""))), data=DataBIOMOD[sp,], family=binomial, weights=Yweights[sp,i])
        		model.sp <- step.gam(gamStart, scope(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]],"s", Spline), keep=functionkeep, direction="both",trace=F, control=gam.control(maxit=50, bf.maxit=50))	
     #New   model.sp <- gam(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeGAM(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]]),collapse="")))), data=DataBIOMOD[sp,], family=binomial, weights=Yweights[sp,i] )
        }                  
        if(Model == 'GBM'){
            model.sp <- gbm.fit(x=DataBIOMOD[sp,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[sp,Biomod.material[["NbVar"]]+i], distribution="bernoulli", w=Yweights[sp,i], var.monotone=rep(0, length=Biomod.material[["NbVar"]]), n.trees=No.trees, 
                         interaction.depth=3, shrinkage=0.01, bag.fraction=0.5, train.fraction=1, verbose=F,  var.names=colnames(DataBIOMOD[1:Biomod.material[["NbVar"]]]),
                         response.name = Biomod.material[["species.names"]][i])

            best.iter <- gbm.perf(model.sp,method="OOB",plot.it=F)
        }	
        if(Model == 'GLM'){
            glmStart <- glm(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]), "~1", collapse=""))), data=DataBIOMOD[sp,], family=binomial, weights=Yweights[sp,i])
            model.sp <- stepAIC(glmStart, scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], Type), direction="both", trace=F, control=glm.control(maxit=100), k=criteria)
        }	
        if(Model == 'MARS'){
            if(is.null(Yweights)) model.sp <- mars(x=DataBIOMOD[sp,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[sp,Biomod.material[["NbVar"]]+i], degree=2)
            else model.sp <- mars(x=DataBIOMOD[sp,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[sp,Biomod.material[["NbVar"]]+i], w=Yweights[sp,i], degree=2)
        }
        if(Model == 'MDA') model.sp <- mda(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "MDA"), collapse="")))), data=DataBIOMOD[sp,], method=mars)
        if(Model == 'RF') model.sp <- randomForest(x=DataBIOMOD[sp,1:Biomod.material[["NbVar"]]], y=as.factor(DataBIOMOD[sp,Biomod.material[["NbVar"]]+i]), ntree=750, mtry=Biomod.material[["NbVar"]]/2, importance=TRUE)

          	
    		if(nrow(Ids)<nrow(DataBIOMOD)){

            if(Model == 'GLM') pred.test <- as.integer(testnull(model.sp, Prev, DataBIOMOD[-sp,]) *1000)
      			if(Model == 'GAM') pred.test <- as.integer(testnull(model.sp, Prev, DataBIOMOD[-sp,]) *1000)
      			if(Model == 'GBM') pred.test <- as.integer(predict.gbm(model.sp, DataBIOMOD[-sp,], best.iter, type='response') *1000)
            if(Model == 'CTA') pred.test <- as.integer(predict(model.sp,DataBIOMOD[-sp,], type="vector") *1000)  
            if(Model == 'SRE') pred.test <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD[sp,]$",paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[sp,1:Biomod.material[["NbVar"]]],DataBIOMOD[-sp,], Perc025, Perc05)) *1000)
            if(Model == 'ANN') pred.test <- as.integer(Rescaler2(predict(model.sp,DataBIOMOD[-sp,], type="raw"), type="range") *1000) 
      			if(Model == 'MDA') pred.test <- as.integer(Rescaler2(predict(model.sp, DataBIOMOD[-sp,1:Biomod.material[["NbVar"]]], type="post")[,2], type="range") *1000)
      			if(Model == 'MARS') pred.test <- as.integer(Rescaler2(predict(model.sp, DataBIOMOD[-sp,1:Biomod.material[["NbVar"]]]), type="range") *1000)
            if(Model == 'RF') pred.test <- as.integer(Rescaler2(predict(model.sp, DataBIOMOD[-sp,1:Biomod.material[["NbVar"]]], type="prob")[,2], type="range") *1000)
            
            if(Model == 'SRE'){ 
                Kappa.train <- Kappa.train + KappaSRE(DataBIOMOD[-sp,Biomod.material[["NbVar"]]+i], pred.test)$Kappa 
                TSS.train <- TSS.train + KappaSRE(DataBIOMOD[-sp,Biomod.material[["NbVar"]]+i], pred.test, TSS=T)$TSS 
            } 
            else{
                AUC.train <- AUC.train + somers2(pred.test, DataBIOMOD[-sp,Biomod.material[["NbVar"]]+i])["C"]
          			Kappa.train <- Kappa.train + KappaRepet(DataBIOMOD[-sp,Biomod.material[["NbVar"]]+i], pred.test)$Kappa
          			TSS.train <- TSS.train + KappaRepet(DataBIOMOD[-sp,Biomod.material[["NbVar"]]+i], pred.test, TSS=T)$TSS
      		  }
      			
        } else{
      			
            if(Model == 'GLM') pred.test <- as.integer(as.numeric(fitted(model.sp)) *1000)          # pb, ici on peut pas faire le testnull()
            if(Model == 'GBM') pred.test <- as.integer(predict.gbm(model.sp, DataBIOMOD, best.iter, type='response') *1000)
            if(Model == 'GAM') pred.test <- as.integer(testnull(model.sp, Prev, DataBIOMOD) *1000)
            if(Model == 'CTA') pred.test <- as.integer(predict(model.sp,DataBIOMOD, type="vector") *1000)
            if(Model == 'SRE') pred.test <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$",paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[,1:Biomod.material[["NbVar"]]],DataBIOMOD, Perc025, Perc05)) *1000)
            if(Model == 'ANN') TempArray <- predict(model.sp, DataBIOMOD, type="raw")
            if(Model == 'MDA') TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]], type="post")[,2]
            if(Model == 'MARS') TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]])
            if(Model == 'RF') TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]], type="prob")[,2]
            
            if(Model == 'ANN' | Model == 'MDA' | Model == 'MARS' | Model == 'RF') pred.test <- as.integer(Rescaler2(TempArray, type="range") *1000)
            
            if(Model == 'SRE'){ 
                Kappa.train <- Kappa.train + KappaSRE(DataBIOMOD[,Biomod.material[["NbVar"]]+i], pred.test)$Kappa 
                TSS.train <- TSS.train + KappaSRE(DataBIOMOD[,Biomod.material[["NbVar"]]+i], pred.test, TSS=T)$TSS 
            } 
            else{
                AUC.train <- AUC.train + somers2(pred.test, DataBIOMOD[,Biomod.material[["NbVar"]]+i])["C"]
          			Kappa.train <- Kappa.train + KappaRepet(DataBIOMOD[,Biomod.material[["NbVar"]]+i], pred.test)$Kappa
          			TSS.train <- TSS.train + KappaRepet(DataBIOMOD[,Biomod.material[["NbVar"]]+i], pred.test, TSS=T)$TSS
      		  }
        }
        if(Model == 'ANN'){
            decay.final <- decay.final + decay
		        size.final <- size.final + size
        }        
    }
	           
	  
    #mean AUC, decay and size from the calibration
    if(Model != 'SRE') AUC.train <- AUC.train/Run
    Kappa.train <- Kappa.train/Run
    TSS.train <- TSS.train/Run
    
    if(Model == 'ANN'){
        decay.final <- decay.final/Run
        size.final <- size.final/Run          
    }    

    #Final model calibrated on all the data
    
    if(nrow(Ids)<nrow(DataBIOMOD)){
        if(Model == 'ANN'){
        		model.sp <- nnet(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "NNET"),collapse="")))), 
        		               data=DataBIOMOD, weights=Yweights[,i], size=size.final,linout=F, entropy=T, skip=F, decay=decay.final, maxit=100, trace=F)
        		TempArray <- predict(model.sp, DataBIOMOD, type="raw")
        		g.pred <- as.integer(Rescaler2(as.numeric(TempArray), type="range") *1000)
  	    }
        if(Model == 'CTA'){
        		temp <- rpart.control(xval=CV.tree, minbucket=5, minsplit=5,cp=0.001, maxdepth=25)
        		model.sp <- rpart(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeExpSyst(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]],"CTA"),collapse ="")))),DataBIOMOD, weights=Yweights[,i], control=temp)
        		tr <- as.data.frame(model.sp$cptable)
        		tr[,6] <- tr[,4] + tr[,5]
        		tr <- tr[tr[,2]!=0,]
        		Cp <- tr[tr[,6] == min(tr[,6]), 1]
        		if(length(Cp)==1) model.sp <- prune(model.sp, cp=Cp)
        		else model.sp <- prune(model.sp, cp=Cp[2])
        		g.pred <- as.integer(as.numeric(predict(model.sp, DataBIOMOD, type="vector")) *1000)
        }
          if(Model == 'GAM'){
          
           	gamStart <- gam(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]), "~1", collapse=""))), data=DataBIOMOD, family=binomial, weights=Yweights[,i])
        	 	model.sp <- step.gam(gamStart, scope(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "s", Spline), keep=functionkeep, direction="both", trace=F, control=gam.control(maxit=50, bf.maxit=50))
      #New  model.sp <- gam(eval(parse(text=paste(Biomod.material[["species.names"]][i],paste(scopeGAM(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]]),collapse="")))), data=DataBIOMOD, family=binomial, weights=Yweights[,i])        	
          	g.pred <- as.integer(as.numeric(testnull(model.sp, Prev, DataBIOMOD)) *1000)
        }
        if(Model == 'GBM'){
        		model.sp <- gbm.fit(x=DataBIOMOD[,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[,Biomod.material[["NbVar"]]+i], distribution="bernoulli", w=Yweights[,i], var.monotone=rep(0, length=Biomod.material[["NbVar"]]), n.trees=No.trees, 
                         interaction.depth=3, shrinkage=0.01, bag.fraction=0.5, train.fraction=1, verbose=F,  var.names=colnames(DataBIOMOD[1:Biomod.material[["NbVar"]]]),
                         response.name = Biomod.material[["species.names"]][i])
            best.iter <- gbm.perf(model.sp,method="OOB",plot.it=F)		
        		g.pred <- as.integer(predict.gbm(model.sp, DataBIOMOD, best.iter, type="response")*1000)
        }
        if(Model == 'GLM'){
        		if(Test == "BIC"){
            		criteria <- log(nrow(DataBIOMOD))
            		cat("\t Stepwise procedure using BIC criteria \n")
        		}
        		else criteria=2
        		glmStart <- glm(eval(parse(text=paste(paste(Biomod.material[["species.names"]][i]), "~1", collapse=""))), data=DataBIOMOD, family=binomial, weights=Yweights[,i])
        		model.sp <- stepAIC(glmStart, scopeExpSyst(DataBIOMOD[1:10,1:Biomod.material[["NbVar"]]], Type), direction="both", trace=F, control=glm.control(maxit=100), k=criteria)
        		g.pred <- as.integer(as.numeric(fitted(model.sp)) *1000)   # pb, ici pas de testnull() possible
        }
        if(Model == 'MARS'){	
        	  if(is.null(Yweights)) model.sp <- mars(x=DataBIOMOD[,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[,Biomod.material[["NbVar"]]+i], degree=2)
            else model.sp <- mars(x=DataBIOMOD[,1:Biomod.material[["NbVar"]]], y=DataBIOMOD[,Biomod.material[["NbVar"]]+i], w=Yweights[,i], degree=2) 
        		TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]])
        		g.pred <- as.integer(Rescaler2(TempArray, type="range") *1000)
        }
        if(Model == 'MDA') {
        		model.sp <- mda(eval(parse(text=paste(Biomod.material[["species.names"]][i], paste(scopeExpSyst(DataBIOMOD[1:10, 1:Biomod.material[["NbVar"]]], "MDA"), collapse="")))), data=DataBIOMOD, method=mars)
        		TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]], type="post")[,2]
        		g.pred <- as.integer(Rescaler2(TempArray, type="range") *1000)
        }	
        if(Model == 'RF'){	
        		model.sp <- randomForest(x=DataBIOMOD[,1:Biomod.material[["NbVar"]]], y=as.factor(DataBIOMOD[,Biomod.material[["NbVar"]]+i]), ntree=750, importance=TRUE, mtry=Biomod.material[["NbVar"]]/2)
        		TempArray <- predict(model.sp, DataBIOMOD[,1:Biomod.material[["NbVar"]]], type="prob")[,2]
        		g.pred <- as.integer(Rescaler2(TempArray, type="range") *1000)
        }
        if(Model == 'SRE') g.pred <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$",paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[,1:Biomod.material[["NbVar"]]],DataBIOMOD, Perc025, Perc05)) *1000)
  	
  	}
    else g.pred <- pred.test #if the data have not been splited.
  	
  	#saving the models on the hard disk
  	
  	if(Model != 'SRE') eval(parse(text=paste("assign('",Biomod.material[["species.names"]][i],"_", Model,"', model.sp)", sep="")))
  	if(Model != 'SRE') eval(parse(text=paste("save(",Biomod.material[["species.names"]][i],"_", Model, ",file='", getwd(), "/models/", Biomod.material[["species.names"]][i],"_", Model,"')", sep="")))
    if(Model == 'GBM') g.list$best.iter <- best.iter
    
    
    # Save Minumum and Maxium of the calibration prediction range:
    if(Model == 'ANN' | Model == 'MDA' | Model == 'MARS' | Model == 'RF') { 
          g.list$RawPred <- c(min(TempArray), max(TempArray))
        	rm(TempArray)
    }
    
    #Evaluation of Predictor Importance in the model:
    if(VarImport > 0){
        cat("Evaluating Predictor Contributions in " , Model, "...", "\n")
        TempVarImp <- as.data.frame(matrix(data=0, nrow=1, ncol=Biomod.material[["NbVar"]]))
        names(TempVarImp) <- names(DataBIOMOD)[1:Biomod.material[["NbVar"]]]
        
        for(J in 1:Biomod.material[["NbVar"]]){
            for(K in 1:VarImport){
                TempDS <- DataBIOMOD[,1:Biomod.material[["NbVar"]]]
                TempDS[,J] <- sample(TempDS[,J])
                
                if(Model == 'ANN') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(Rescaler2(as.numeric(predict(model.sp, TempDS, type="raw")), type="range", OriMinMax=g.list$RawPred) *1000))
                if(Model == 'CTA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(as.numeric(predict(model.sp, TempDS, type="vector")) *1000))
                if(Model == 'GLM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GAM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(as.numeric(testnull(model.sp, Prev, TempDS)) *1000))
                if(Model == 'GBM') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(as.numeric(predict.gbm(model.sp, TempDS, best.iter, type='response')) *1000))
                if(Model == 'MARS') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(Rescaler2(predict(model.sp, TempDS[,1:Biomod.material[["NbVar"]]]), type="range", OriMinMax=g.list$RawPred) *1000))
                if(Model == 'MDA') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(Rescaler2(predict(model.sp, TempDS, type="post")[,2], type="range", OriMinMax=g.list$RawPred) *1000))
                if(Model == 'RF') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(Rescaler2(predict(model.sp, TempDS, type="prob")[,2], type="range", OriMinMax=g.list$RawPred) *1000))
                if(Model == 'SRE') TempVarImp[1,J] <- TempVarImp[1,J] + cor(g.pred, as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$",paste(Biomod.material[["species.names"]][i]), collapse=""))), TempDS,DataBIOMOD, Perc025, Perc05)) *1000))    
           }
        }
        VarImportance[[i]][Model,] <-  round(1 - (TempVarImp/VarImport), digits=3)
    }
	 
	 
    #Model on the independent data if any.
    if(exists("DataEvalBIOMOD")){
        if(Model == 'GLM') predtest <- as.integer(as.numeric(testnull(model.sp, Prev, DataEvalBIOMOD)) *1000)
        if(Model == 'GAM') predtest <- as.integer(as.numeric(testnull(model.sp, Prev, DataEvalBIOMOD)) *1000)
        if(Model == 'ANN') predtest <- as.integer(Rescaler2(as.numeric(predict(model.sp, DataEvalBIOMOD, type="raw")), type="range", OriMinMax=g.list$RawPred) *1000)
        if(Model == 'CTA') predtest <- as.integer(as.numeric(predict(model.sp, DataEvalBIOMOD, type="vector")) *1000)
        if(Model == 'GBM') predtest <- as.integer(predict.gbm(model.sp, DataEvalBIOMOD, best.iter, type='response') *1000)        
        if(Model == 'MARS') predtest <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]]), type="range", OriMinMax=g.list$RawPred) *1000)
        if(Model == 'MDA') predtest <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]], type="post")[,2], type="range", OriMinMax=g.list$RawPred) *1000)
        if(Model == 'RF') predtest <- as.integer(Rescaler2(predict(model.sp, DataEvalBIOMOD[,1:Biomod.material[["NbVar"]]], type="prob")[,2], type="range", OriMinMax=g.list$RawPred) *1000)
        if(Model == 'SRE') predtest <- as.integer(as.numeric(sre(eval(parse(text=paste("DataBIOMOD$", paste(Biomod.material[["species.names"]][i]), collapse=""))), DataBIOMOD[1:Biomod.material[["NbVar"]]], DataEvalBIOMOD, Perc025, Perc05)) *1000)
    }


    if(Roc && Model != 'SRE'){
        Evaluation.results.Roc[[i]][Model,] <- c(round(AUC.train, digits=3), 'none', round(somers2(g.pred, DataBIOMOD[,Biomod.material[["NbVar"]]+i])["C"], digits=3), round(CutOff.Optimised(DataBIOMOD[,Biomod.material[["NbVar"]]+i], g.pred), digits=3))  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.Roc[[i]][Model,2] <- round(somers2(predtest,DataEvalBIOMOD[,Biomod.material[["NbVar"]] +i])["C"], digits=3) 
        assign("Evaluation.results.Roc", Evaluation.results.Roc, pos=1)
	  }    
    if(Kappa){
	      Evaluation.results.Kappa[[i]][Model,] <- c(round(Kappa.train, digits=3), 'none', KappaRepet(DataBIOMOD[,Biomod.material[["NbVar"]]+i], g.pred, TSS=T)[c(1,2,4,6)]) 
        if(exists("DataEvalBIOMOD")) Evaluation.results.Kappa[[i]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,Biomod.material[["NbVar"]]+i], predtest)$Kappa, digits=3) 
        assign("Evaluation.results.Kappa", Evaluation.results.Kappa, pos=1)
    }	
    if(TSS){
	      Evaluation.results.TSS[[i]][Model,] <- c(round(TSS.train, digits=3), 'none', KappaRepet(DataBIOMOD[,Biomod.material[["NbVar"]]+i], g.pred, TSS=T)[c(1,2,4,6)])  
	      if(exists("DataEvalBIOMOD")) Evaluation.results.TSS[[i]][Model,2] <- round(KappaRepet(DataEvalBIOMOD[,Biomod.material[["NbVar"]]+i], predtest, TSS=T)$TSS, digits=3) 
        assign("Evaluation.results.TSS", Evaluation.results.TSS, pos=1)
    }
    	
  	if(exists("DataEvalBIOMOD") && KeepPredIndependent){
    		assign("g.pred", predtest, pos=1)
  	}
    
    assign("g.pred", g.pred, pos=1)
    assign("VarImportance", VarImportance, pos=1)
    
    if(Model != 'SRE') rm(model.sp)
    
    return(g.list)
}

