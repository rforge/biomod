`Migration` <-
function(CurrentPred=NULL, FutureProj=NULL, X=NULL, Y=NULL, MaxMigr=NULL, Pred.Save=NULL)
{
    #if(dim(CurrentPred)[1] != dim(FutureProj)[1])
    #stop("\n\n Number of pixel should be the same in CurrentPred, FutureProj, X and Y")
    #if(dim(CurrentPred)[2] != dim(FutureProj)[2])
    #stop("\n\n Species number should be the same in CurrentPred and FutureProj")
    
    if(sum(CurrentPred[,1] >1) >=1 | sum(FutureProj[,1] >1) >=1) stop("\n Data inputs should be binary")
    
    PredMigr <- as.data.frame(matrix(0, nrow=nrow(CurrentPred)))    
    
    if(length(MaxMigr) == 1 | length(MaxMigr) == ncol(CurrentPred)){
        if(length(MaxMigr) == 1) MaxMigr <- rep(MaxMigr, ncol(CurrentPred))  #duplicate to nb of columns extent
    
        for(j in 1:ncol(CurrentPred)){
            if(sum(FutureProj[CurrentPred[, j] == 0, j]) != 0) 
            PredMigr[, j] <- Dispersal.Limit(CurrentPred[, j], FutureProj[, j], X, Y, MaxMigr[j])
            else PredMigr[, j] <- FutureProj[, j]
        }    
    } else stop("\n Incorrect values in MaxMigr, does not match the number of columns of data input")
       
    #Saving the results; checking if object already created   
    if(!is.null(Pred.Save)){
        if(exists(Pred.Save, where=1)) {
            New.Pred.Save <- unique(Pred.Save, pos=1)
            assign(New.Pred.Save, PredMigr, pos=1)
            warning(paste("\n Prediction with dispersal saved in", New.Pred.Save, "\n -> object overwritten \n"))
        }
        else assign(Pred.Save, PredMigr, pos=1)
    }
}

