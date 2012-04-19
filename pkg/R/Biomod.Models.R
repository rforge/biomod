.Biomod.Models <- function (Model, Ids, PA.samp, TypeGLM, Test, No.trees, CV.tree, 
    CV.ann, quant, NbRunEval, Spline, DataSplit, Yweights, Roc, 
    Optimized.Threshold.Roc, Kappa, TSS, KeepPredIndependent, 
    VarImport) 
{
  
	
	
    if (Model == "GLM") {
        if (TypeGLM == "poly") {
            Type <- "GLMp"
            cat("Model=GLM polynomial + quadratic")
        }
        if (TypeGLM == "quad") {
            Type <- "GLMq"
            cat("Model=GLM quadratic \n")
        }
        if (TypeGLM == "simple") {
            Type <- "GLMs"
            cat("Model=GLM linear \n")
        }
        if (Test == "AIC") {
            criteria <- 2
            cat("\t Stepwise procedure using AIC criteria \n")
        }
        if (Test == "BIC") {
            criteria <- log(length(PA.samp))
            cat("\t Stepwise procedure using BIC criteria \n")
        }
    }
    if (Model == "GBM") {
        cat("Model=Generalised Boosting Regression \n")
        cat("\t", No.trees, "maximum different trees and lambda Fold Cross-Validation \n")
    }
    if (Model == "GAM") {
        cat("Model=GAM spline \n")
        cat("\t", Spline, " Degrees of smoothing \n")
    }
    if (Model == "CTA") {
        cat("Model=Classification tree \n")
        cat("\t", CV.tree, "Fold Cross-Validation \n")
    }
    if (Model == "ANN") {
        cat("Model=Artificial Neural Network \n")
        cat("\t", CV.ann, "Fold Cross Validation + 3 Repetitions \n")
        cat("Calibration and evaluation phase: Nb of cross-validations: ", 
            ncol(Ids), "\n")
    }
    if (Model == "SRE") 
        cat("Model=Surface Range Envelop \n")
    if (Model == "FDA"){
      cat("Model=Flexible Discriminant Analysis \n")
    } 

    if (Model == "MARS") 
        cat("Model=Multiple Adaptive Regression Splines \n")
    if (Model == "RF") 
        cat("Model=Breiman and Cutler's random forests for classification and regression \n")
    if (Model == "GLM" | Model == "GAM") 
        Prev <- sum(DataBIOMOD[, i + Biomod.material$NbVar])/nrow(DataBIOMOD)
    if (Model == "CTA") 
        set.seed(123)
    if (Model == "RF") 
        set.seed(71)
    if (Roc & Model != "SRE") 
        AUC.train <- 0
    if (Kappa) 
        Kappa.train <- 0
    if (TSS) 
        TSS.train <- 0
    
    SpNames <- Biomod.material$species.names
    NbVar <- Biomod.material$NbVar
    ErrorCounter <- 0
    if (Model == "GBM") 
        GBM.list <- list()
    RunWeights <- Yweights[, i]
    for (k in 1:(ncol(Ids) + 1)) {
        if (exists("model.sp")) 
            rm(model.sp)
        if (exists("g.pred")) 
            rm(g.pred)
        if (k == (ncol(Ids) + 1)) {
            if (Biomod.material$NbRepPA == 0) 
                nam <- "full"
            else nam <- paste("PA", pa, sep = "")
            calib.lines <- pred.lines <- PA.samp
        }
        else {
            if (Biomod.material$NbRepPA == 0) 
                nam <- paste("full_rep", k, sep = "")
            else nam <- paste("PA", pa, "_rep", k, sep = "")
            calib.lines <- PA.samp[Ids[, k]]
            pred.lines <- PA.samp[-Ids[, k]]
        }
      
        
        assign("calib.lines", calib.lines, pos = 1)
        if (Model == "ANN") {
            set.seed(555)
            CV_nnet = .CV.nnet(Input = data.frame(matrix(DataBIOMOD[calib.lines,
                1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), Target = DataBIOMOD[calib.lines, NbVar +
                i], nbCV = CV.ann, W = RunWeights[calib.lines])
            if (k == (ncol(Ids) + 1)) {
                if (is.null(Yweights)) 
                	model.sp <- nnet(eval(parse(text=paste(paste("DataBIOMOD[calib.lines, NbVar +",i,"]"),paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10,1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))),"GBM"),collapse="")))),
                data = DataBIOMOD[calib.lines,], size = CV_nnet[1,
                    1], rang = 0.1, decay = CV_nnet[1, 2], maxit = 200, 
                  trace = FALSE)
                else model.sp <- nnet(eval(parse(text=paste(paste("DataBIOMOD[calib.lines, NbVar +",i,"]"),paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10,1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))),"GBM"),collapse="")))),
                data = DataBIOMOD[calib.lines,], size = CV_nnet[1,
                    1], rang = 0.1, decay = CV_nnet[1, 2], weights=RunWeights[calib.lines], maxit = 200, 
                  trace = FALSE)
            }        
            else {
                if (is.null(Yweights)) 
            		try(model.sp <- nnet(eval(parse(text=paste("DataBIOMOD[calib.lines, NbVar +",i,"]",paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10,1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))),"GBM"),collapse="")))),
                data = DataBIOMOD[calib.lines,],
                size = CV_nnet[1, 1], rang = 0.1, decay = CV_nnet[1, 
                  2], maxit = 200, trace = FALSE, weights=RunWeights[calib.lines]), silent = TRUE)
                  
                else  try(model.sp <- nnet(eval(parse(text=paste(paste("DataBIOMOD[calib.lines, NbVar +",i,"]"),paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10,1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))),"GBM"),collapse="")))),
                data = DataBIOMOD[calib.lines,],
                size = CV_nnet[1, 1], rang = 0.1, decay = CV_nnet[1, 
                  2], maxit = 200, trace = FALSE, weights=RunWeights[calib.lines]), silent = TRUE)

			}  

            if (exists("model.sp")) 
                TempArray <- predict(model.sp, DataBIOMOD[PA.samp, 
                  ], type = "raw")
        }
        if (Model == "CTA") {
            temp <- rpart.control(xval = CV.tree, minbucket = 5, 
                minsplit = 5, cp = 0.001, maxdepth = 25)
            if (k == (ncol(Ids) + 1)) 
                try(model.sp <- rpart(eval(parse(text = paste("as.factor(", 
                  SpNames[i], ")", paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 
                    1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), "CTA"), collapse = "")))), DataBIOMOD[calib.lines, 
                  ], weights = RunWeights[calib.lines], control = temp))
            else try(model.sp <- rpart(eval(parse(text = paste("as.factor(", 
                SpNames[i], ")", paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), "CTA"), collapse = "")))), DataBIOMOD[calib.lines, 
                ], weights = RunWeights[calib.lines], control = temp), 
                silent = TRUE)
            if (exists("model.sp")) {
                tr <- as.data.frame(model.sp$cptable)
                tr[, 6] <- tr[, 4] + tr[, 5]
                tr <- tr[tr[, 2] != 0, ]
                Cp <- tr[tr[, 6] == min(tr[, 6]), 1]
                if (length(Cp) == 1) 
                  model.sp <- prune(model.sp, cp = Cp)
                else model.sp <- prune(model.sp, cp = Cp[2])
                g.pred <- data.frame(as.integer(as.numeric(predict(model.sp, 
                  data.frame(matrix(DataBIOMOD[PA.samp, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), type = "prob")[, 
                  2]) * 1000))
              
            }
        }
        if (Model == "GAM") {
            gamStart <- gam(eval(parse(text = paste(paste(SpNames[i]), 
                "~1", collapse = ""))), data = DataBIOMOD[calib.lines, 
                ], family = binomial, weights = RunWeights[calib.lines])
            if (k == (ncol(Ids) + 1)) 
                model.sp <- step.gam(gamStart, .scope(data.frame(matrix(DataBIOMOD[1:10, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), "s", Spline), keep = .functionkeep, 
                  direction = "both", trace = FALSE, control = gam.control(maxit = 50, 
                    bf.maxit = 50))
            else try(model.sp <- step.gam(gamStart, .scope(data.frame(matrix(DataBIOMOD[1:10, 
                1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), "s", Spline), keep = .functionkeep, 
                direction = "both", trace = FALSE, control = gam.control(maxit = 50, 
                  bf.maxit = 50)), silent = TRUE)
            if (exists("model.sp")) 
                g.pred <- data.frame(as.integer(as.numeric(.testnull(model.sp, 
                  Prev, DataBIOMOD[PA.samp, ])) * 1000))
        }
        if (Model == "GBM") {
            if (k == (ncol(Ids) + 1)) 
                model.sp <- gbm(eval(parse(text = paste(SpNames[i], 
                  paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                    "GBM"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                  ], distribution = "bernoulli", var.monotone = rep(0, 
                  length = NbVar), w = RunWeights[calib.lines], interaction.depth = 7, shrinkage = 0.001, 
                  bag.fraction = 0.5, train.fraction = 1, n.trees = No.trees, 
                  verbose = FALSE, cv.folds = 5)
            else try(model.sp <- gbm(eval(parse(text = paste(SpNames[i], 
                paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                  "GBM"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                ], distribution = "bernoulli", var.monotone = rep(0, 
                length = NbVar), w = RunWeights[calib.lines], 
                interaction.depth = 7, shrinkage = 0.001, bag.fraction = 0.5, 
                train.fraction = 1, n.trees = No.trees, verbose = FALSE, 
                cv.folds = 5), silent = TRUE)
            if (exists("model.sp")) {
                best.iter <- gbm.perf(model.sp, method = "cv", 
                  plot.it = FALSE)
                g.pred <- data.frame(as.integer(predict.gbm(model.sp, 
                  DataBIOMOD[PA.samp, ], best.iter, type = "response") * 
                  1000))
            }
        }
        if (Model == "GLM") {
            glmStart <- glm(eval(parse(text = paste(paste(SpNames[i]), 
                "~1", collapse = ""))), data = DataBIOMOD[calib.lines, 
                ], family = binomial, weights = RunWeights[calib.lines])
            if (k == (ncol(Ids) + 1)) 
                model.sp <- stepAIC(glmStart, .scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), Type), direction = "both", trace = FALSE, 
                  control = glm.control(maxit = 100), k = criteria)
            else try(model.sp <- stepAIC(glmStart, .scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 
                1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), Type), direction = "both", trace = FALSE, 
                control = glm.control(maxit = 100), k = criteria), 
                silent = TRUE)
            if (exists("model.sp")) 
                g.pred <- data.frame(as.integer(.testnull(model.sp, 
                  Prev, DataBIOMOD[PA.samp, ]) * 1000))
        }
        if (Model == "MARS") {
            if (k == (ncol(Ids) + 1)) {
                if (is.null(Yweights)) 
                  model.sp <- mars(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                    1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = DataBIOMOD[calib.lines, NbVar + 
                    i], degree = 2)
                else model.sp <- mars(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = DataBIOMOD[calib.lines, NbVar + 
                  i], degree = 2, w = RunWeights[calib.lines])
            }
            else {
                if (is.null(Yweights)) 
                  try(model.sp <- mars(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                    1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = DataBIOMOD[calib.lines, NbVar + 
                    i], degree = 2), silent = TRUE)
                else try(model.sp <- mars(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = DataBIOMOD[calib.lines, NbVar + 
                  i], degree = 2, w = RunWeights[calib.lines]), 
                  silent = TRUE)
            }
            if (exists("model.sp")) 
                TempArray <- predict(model.sp, data.frame(matrix(DataBIOMOD[PA.samp, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))))
        }
        if (Model == "FDA") {
#           counter <- 1
#           g.pred <- NA
#           while(counter < 5 | sum(is.na(g.pred)) > 0 ){
#             if(counter > 1) {cat('\n*** last modeling failed, re-runing the model')}
            if (k == (ncol(Ids) + 1)) {
                if (is.null(Yweights)) 
                  model.sp <- fda(eval(parse(text = paste(SpNames[i], 
                    paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                      "FDA"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                    ], method = mars)
                else model.sp <- fda(eval(parse(text = paste(SpNames[i], 
                  paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                    "FDA"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                  ], method = mars, weights = RunWeights[calib.lines])
            }
            else {
                if (is.null(Yweights)) 
                  try(model.sp <- fda(eval(parse(text = paste(SpNames[i], 
                    paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                      "FDA"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                    ], method = mars), silent = TRUE)
                else try(model.sp <- fda(eval(parse(text = paste(SpNames[i], 
                  paste(.scopeExpSyst(data.frame(matrix(DataBIOMOD[1:10, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                    "FDA"), collapse = "")))), data = DataBIOMOD[calib.lines, 
                  ], method = mars, weights = RunWeights[calib.lines]), silent = TRUE)
            }
            if (exists("model.sp")) 
                TempArray <- predict(model.sp, data.frame(matrix(DataBIOMOD[PA.samp, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), type = "post")[, 2]
#             counter <- counter + 1
#           }
        }
        if (Model == "RF") {
            if (k == (ncol(Ids) + 1)) 
                model.sp <- randomForest(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                  1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = as.factor(DataBIOMOD[calib.lines, 
                  NbVar + i]), ntree = 750, mtry = NbVar/2, importance = TRUE)
            else try(model.sp <- randomForest(x = data.frame(matrix(DataBIOMOD[calib.lines, 
                1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), y = as.factor(DataBIOMOD[calib.lines, 
                NbVar + i]), ntree = 750, mtry = NbVar/2, importance = TRUE), 
                silent = TRUE)
            if (exists("model.sp")) 
                g.pred <- data.frame(as.integer(predict(model.sp, 
                  data.frame(matrix(DataBIOMOD[PA.samp, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), type = "prob")[, 
                  2] * 1000))
        }
        if (Model == "SRE") 
            g.pred <- data.frame(as.integer(as.numeric(sre(eval(parse(text = paste("DataBIOMOD[calib.lines,]$", 
                paste(SpNames[i]), collapse = ""))), data.frame(matrix(DataBIOMOD[calib.lines, 
                1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), DataBIOMOD[PA.samp, ], quant)) * 1000))
        if (any(c("ANN", "FDA", "MARS") == Model)){
            if (exists("TempArray")){
              g.pred <- data.frame(as.integer(.Rescaler4(as.numeric(TempArray),
                  ref = DataBIOMOD[PA.samp, NbVar + i], run = paste(SpNames[i],
                  "_", Model, "_", nam, sep = ""), original = TRUE) *
                  1000))
            }
        }
            
        ### checking for NA in g.pred #########################
        if(exists("g.pred")){
          if(sum(is.na(g.pred[,])) > 0){ # some NA detected
            if(sum(is.na(g.pred[,])) == ncol(g.pred) * nrow(g.pred)){ # only NA 
              rm('g.pred') # model failed
            } else { # only few NAs
              na.detected <- which(is.na(g.pred))
              # remove Na's pred lines of all data
              if (k == (ncol(Ids) + 1)) {
                calib.lines <- pred.lines <- PA.samp[-na.detected]
              } else {
                if(sum( as.vector(Ids[, k]) %in% na.detected ) > 0){
                  # remove from Ids
                  Ids <- data.frame(Ids[Ids[,k] %in% na.detected,])
                }
                PA.samp <- PA.samp[-na.detected]
                calib.lines <- PA.samp[Ids[, k]]
                pred.lines <- PA.samp[-Ids[, k]]
                g.pred.with.na <- g.pred
                g.pred <- na.omit(g.pred)
              }
            }
          }
        }

            
        if (!exists("g.pred")) {
            ErrorCounter <- ErrorCounter + 1
            BM[["calibration.failures"]] <- c(BM[["calibration.failures"]], 
                paste(SpNames[i], "_", Model, "_", nam, sep = ""))
        }
        if (exists("g.pred")) {
            if (k != (ncol(Ids) + 1)) {
                if (Roc) 
                  auc.stat <- .somers2(g.pred[-Ids[, k], ], DataBIOMOD[pred.lines, 
                    NbVar + i])["C"]
                if (Kappa) 
                  kappa.stat <- KappaRepet(DataBIOMOD[pred.lines, 
                    NbVar + i], g.pred[-Ids[, k], ])$Kappa
                if (TSS) 
                  tss.stat <- KappaRepet(DataBIOMOD[pred.lines, 
                    NbVar + i], g.pred[-Ids[, k], ], TSS = TRUE)$TSS
            }
            if (k != (ncol(Ids) + 1)) {
                if (Model == "SRE") {
                  if (Kappa) 
                    Kappa.train <- Kappa.train + .KappaSRE(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[-Ids[, k], ])$Kappa
                  if (TSS) 
                    TSS.train <- TSS.train + .KappaSRE(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[-Ids[, k], ], TSS = TRUE)$TSS
                }
                else {
                  if (Roc) 
                    AUC.train <- AUC.train + auc.stat
                  if (Kappa) 
                    Kappa.train <- Kappa.train + kappa.stat
                  if (TSS) 
                    TSS.train <- TSS.train + tss.stat
                }
            }
            else {
                if (Model == "SRE") {
                  if (Kappa) 
                    Kappa.final <- .KappaSRE(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[, ])$Kappa
                  if (TSS) 
                    TSS.final <- .KappaSRE(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[, ], TSS = TRUE)$TSS
                }
                else {
                  if (Roc) 
                    AUC.final <- .somers2(g.pred[, ], DataBIOMOD[pred.lines, 
                      NbVar + i])["C"]
                  if (Kappa) 
                    Kappa.final <- KappaRepet(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[, ])$Kappa
                  if (TSS) 
                    TSS.final <- KappaRepet(DataBIOMOD[pred.lines, 
                      NbVar + i], g.pred[, ], TSS = TRUE)$TSS
                }
            }
            if (k != (ncol(Ids) + 1)) {
                if (Roc && Model != "SRE") 
                  Evaluation.results.Roc[[paste(SpNames[i], "_", 
                    nam, sep = "")]][Model, ] <- c(round(auc.stat, 
                    digits = 3), "none", round(.somers2(g.pred[, 
                    ], DataBIOMOD[PA.samp, NbVar + i])["C"], 
                    digits = 3), round(CutOff.Optimised(DataBIOMOD[PA.samp, 
                    NbVar + i], g.pred[, ]), digits = 3))
                if (Kappa) 
                  Evaluation.results.Kappa[[paste(SpNames[i], 
                    "_", nam, sep = "")]][Model, ] <- c(round(kappa.stat, 
                    digits = 3), "none", KappaRepet(DataBIOMOD[PA.samp, 
                    NbVar + i], g.pred[, ])[c(1, 2, 4, 6)])
                if (TSS) 
                  Evaluation.results.TSS[[paste(SpNames[i], "_", 
                    nam, sep = "")]][Model, ] <- c(round(tss.stat, 
                    digits = 3), "none", KappaRepet(DataBIOMOD[PA.samp, 
                    NbVar + i], g.pred[, ], TSS = TRUE)[c(1, 2, 
                    4, 6)])
            }
            if (Model == "GBM") 
                eval(parse(text = paste("GBM.list$", nam, " <- best.iter", 
                  sep = "")))
                
            # save g.pred
            if (k == (ncol(Ids) + 1)){
              if(exists('g.pred.with.na')){ # na in g.pred
                Array[, Model, 1, pa] <- g.pred.with.na[,]
              } else{
                Array[, Model, 1, pa] <- g.pred[, ]
              }
            } else {
              if(exists('g.pred.with.na')){ # na in g.pred
                Array[, Model, (k + 1), pa] <- g.pred.with.na[,]
              } else {
                Array[, Model, (k + 1), pa] <- g.pred[, ]
              }
            }

           # if (Model != "SRE") 
           #    eval(parse(text = paste("assign('", SpNames[i], 
           #      "_", Model, "_", nam, "', model.sp, pos=1)", sep = "")))
          #  if (Model != "SRE") 
           #     eval(parse(text = paste("save(", SpNames[i], 
            #      "_", Model, "_", nam, ",file='", getwd(), "/models/", 
             #     SpNames[i], "_", Model, "_", nam, "', compress='xz')", sep = "")))
           if (Model != "SRE"){ 
            	eval(parse(text = paste(SpNames[i], "_", Model, "_", nam, "= model.sp", sep = "")))       
                eval(parse(text = paste("save(", SpNames[i], 
                  "_", Model, "_", nam, ",file='", getwd(), "/models/", 
                  SpNames[i], "_", Model, "_", nam, "', compress='xz')", sep = "")))
           		eval(parse(text = paste("rm(", SpNames[i], "_", Model, "_", nam, ")", sep = "")))
           }       
                 
        }
    }
    assign("Array", Array, pos = 1)
    assign("BM", BM, pos = 1)
    if (Model == "GBM") 
        assign("GBM.list", GBM.list, pos = 1)
    if (Roc & Model != "SRE") 
        AUC.train <- AUC.train/(ncol(Ids) - ErrorCounter)
    if (Kappa) 
        Kappa.train <- Kappa.train/(ncol(Ids) - ErrorCounter)
    if (TSS) 
        TSS.train <- TSS.train/(ncol(Ids) - ErrorCounter)
    if (VarImport > 0) {
        cat("Evaluating Predictor Contributions in ", Model, 
            "...", "\n")
        TempVarImp <- as.data.frame(matrix(data = 0, nrow = 1, 
            ncol = Biomod.material$NbVar))
        names(TempVarImp) <- names(DataBIOMOD)[1:Biomod.material$NbVar]
        for (J in 1:Biomod.material$NbVar) {
            for (K in 1:VarImport) {
                TempDS <- DataBIOMOD[PA.samp, 1:Biomod.material$NbVar]
                TempDS[, J] <- sample(TempDS[, J])
                if (Model == "ANN") {
                  set.seed(555)
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(.Rescaler4(as.numeric(predict(model.sp, 
                    TempDS, type = "raw")), ref = DataBIOMOD[PA.samp, 
                    NbVar + i], run = paste(SpNames[i], "_ANN_", 
                    nam, sep = "")) * 1000))
                }    
                if (Model == "CTA") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(predict(model.sp, 
                    TempDS, type = "prob")[, 2]) * 1000))
                if (Model == "GLM") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(.testnull(model.sp, 
                    Prev, TempDS)) * 1000))
                if (Model == "GAM") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(.testnull(model.sp, 
                    Prev, TempDS)) * 1000))
                if (Model == "GBM") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(predict.gbm(model.sp, 
                    TempDS, best.iter, type = "response")) * 
                    1000))
                if (Model == "MARS") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(.Rescaler4(predict(model.sp, 
                    data.frame(matrix(TempDS[, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames)))), ref = DataBIOMOD[PA.samp, 
                    NbVar + i], run = paste(SpNames[i], "_MARS_", 
                    nam, sep = "")) * 1000))
                if (Model == "FDA") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(.Rescaler4(predict(model.sp, 
                    TempDS, type = "post")[, 2], ref = DataBIOMOD[PA.samp, 
                    NbVar + i], run = paste(SpNames[i], "_FDA_", 
                    nam, sep = "")) * 1000))
                if (Model == "RF") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(predict(model.sp, 
                    TempDS, type = "prob")[, 2]) * 1000))
                if (Model == "SRE") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(sre(eval(parse(text = paste("DataBIOMOD[PA.samp,]$", 
                    paste(SpNames[i]), collapse = ""))), TempDS, 
                    DataBIOMOD[PA.samp, ], quant)) * 1000))
            }
        }
        VarImportance[[i]][Model, ] <- round(as.numeric(1 - (TempVarImp/VarImport)), 
            digits = 3)
    }
    assign("VarImportance", VarImportance, pos = 1)
    if (exists("DataEvalBIOMOD")) {
        if (Model == "GLM") 
            predind <- .testnull(model.sp, Prev, DataEvalBIOMOD)
        if (Model == "GAM") 
            predind <- .testnull(model.sp, Prev, DataEvalBIOMOD)
        if (Model == "ANN") {
            set.seed(555)
            predind <- predict(model.sp, DataEvalBIOMOD, type = "raw")
        }
        if (Model == "CTA") 
            predind <- predict(model.sp, data.frame(matrix(DataEvalBIOMOD[, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                type = "prob")[, 2]
        if (Model == "GBM") 
            predind <- predict.gbm(model.sp, DataEvalBIOMOD, 
                best.iter, type = "response")
        if (Model == "MARS") 
            predind <- predict(model.sp, data.frame(matrix(DataEvalBIOMOD[, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))))
        if (Model == "FDA") 
            predind <- predict(model.sp, data.frame(matrix(DataEvalBIOMOD[, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                type = "post")[, 2]
        if (Model == "RF") 
            predind <- predict(model.sp, data.frame(matrix(DataEvalBIOMOD[, 1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                type = "prob")[, 2]
        if (Model == "SRE") 
            predind <- as.integer(as.numeric(sre(eval(parse(text = paste("DataBIOMOD$", 
                paste(SpNames[i]), collapse = ""))), data.frame(matrix(DataBIOMOD[1:NbVar],ncol=Biomod.material$NbVar, dimnames=list(NULL,Biomod.material$VarNames))), 
                DataEvalBIOMOD, quant)) * 1000)
        if (any(c("ANN", "FDA", "MARS") == Model)) 
            predind <- .Rescaler4(predind, ref = DataBIOMOD[PA.samp, 
                NbVar + i], run = paste(SpNames[i], "_", Model, 
                "_", nam, sep = ""))
        predind <- as.integer(as.numeric(predind) * 1000)
    }
    if (Roc && Model != "SRE") {
        Evaluation.results.Roc[[paste(Biomod.material$species.names[i], 
            "_", nam, sep = "")]][Model, ] <- c(round(AUC.train, 
            digits = 3), "none", round(.somers2(g.pred[, ], DataBIOMOD[PA.samp, 
            NbVar + i])["C"], digits = 3), round(CutOff.Optimised(DataBIOMOD[PA.samp, 
            NbVar + i], g.pred[, ]), digits = 3))
        if (exists("DataEvalBIOMOD")) 
            Evaluation.results.Roc[[paste(Biomod.material$species.names[i], 
                "_", nam, sep = "")]][Model, 2] <- round(.somers2(predind, 
                DataEvalBIOMOD[, NbVar + i])["C"], digits = 3)
        assign("Evaluation.results.Roc", Evaluation.results.Roc, 
            pos = 1)
    }
    if (Kappa) {
        Evaluation.results.Kappa[[paste(Biomod.material$species.names[i], 
            "_", nam, sep = "")]][Model, ] <- c(round(Kappa.train, 
            digits = 3), "none", KappaRepet(DataBIOMOD[PA.samp, 
            NbVar + i], g.pred[, ])[c(1, 2, 4, 6)])
        if (exists("DataEvalBIOMOD")) 
            Evaluation.results.Kappa[[paste(Biomod.material$species.names[i], 
                "_", nam, sep = "")]][Model, 2] <- round(KappaRepet(DataEvalBIOMOD[, 
                NbVar + i], predind)$Kappa, digits = 3)
        assign("Evaluation.results.Kappa", Evaluation.results.Kappa, 
            pos = 1)
    }
    if (TSS) {
        Evaluation.results.TSS[[paste(Biomod.material$species.names[i], 
            "_", nam, sep = "")]][Model, ] <- c(round(TSS.train, 
            digits = 3), "none", KappaRepet(DataBIOMOD[PA.samp, 
            NbVar + i], g.pred[, ], TSS = TRUE)[c(1, 2, 4, 6)])
        if (exists("DataEvalBIOMOD")) 
            Evaluation.results.TSS[[paste(Biomod.material$species.names[i], 
                "_", nam, sep = "")]][Model, 2] <- round(KappaRepet(DataEvalBIOMOD[, 
                NbVar + i], predind, TSS = TRUE)$TSS, digits = 3)
        assign("Evaluation.results.TSS", Evaluation.results.TSS, 
            pos = 1)
    }
    if (exists("DataEvalBIOMOD") && KeepPredIndependent) 
        assign("predind", predind, pos = 1)
    assign("Biomod.material", Biomod.material, pos = 1)
 
}