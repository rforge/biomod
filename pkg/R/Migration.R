`Migration` <-
function(CurrentPred=NULL, FutureProj=NULL, X=NULL, Y=
NULL, MaxMigr=NULL, Pred.Save=NULL)
{
    #if(dim(CurrentPred)[1] != dim(FutureProj)[1])
    #stop("\n\n Number of pixel should be the same in CurrentPred, FutureProj, X and Y"
    # )
    #if(dim(CurrentPred)[2] != dim(FutureProj)[2])
    # stop("\n\n Species number should be the same in CurrentPred and FutureProj"
    # )
    NrowCur <- nrow(CurrentPred)
    NcolCur <- ncol(CurrentPred)
    NmaxMigr <- length(MaxMigr)
    PredMigr <- as.data.frame(matrix(0, nrow=NrowCur))
    if(NmaxMigr == 1) {
        for(j in 1:NcolCur)
            if(sum(FutureProj[CurrentPred[, j] == 0, j]) != 0) {
                PredMigr[, j] <- Dispersal.Limit(CurrentPred[
                    , j], FutureProj[, j], X, Y, MaxMigr)
            }
            else PredMigr[, j] <- FutureProj[, j]
    }
    else {
        if(NmaxMigr == NcolCur) {
            for(j in 1:NcolCur) {
                if(sum(FutureProj[CurrentPred[, j] == 0, j]) !=
                    0) {
                    PredMigr[, j] <- Dispersal.Limit(
                        CurrentPred[, j], FutureProj[
                        , j], X, Y, MaxMigr[j])
                }
                else PredMigr[, j] <- FutureProj[, j]
            }
        }
        else stop("\n\nDispersal limit not executed. NbPixel should be either a number or a vector in the same order than species column of CurrentPred and FuturPred")
    }
    if(!is.null(Pred.Save)) {
        if(exists(Pred.Save, where=1)) {
            New.Pred.Save <- unique(Pred.Save, pos=1)
            assign(New.Pred.Save, PredMigr, pos=1)
            warning(paste("\n\nPrediction with dispersal saved in",
                New.Pred.Save))
        }
        else assign(Pred.Save, PredMigr, pos=1)
    }
}

