`Dispersal.Limit` <-
function(CurrentPred=NULL, FutureProj=NULL, X=NULL, Y=NULL, MaxMigr=NULL)
{
    Data <- FutureProj - 2 * CurrentPred
    subX <- X[CurrentPred == 1]
    subY <- Y[CurrentPred == 1]
    subX1 <- X[Data == 1]
    subY1 <- Y[Data == 1]
    Cs <- sqrt(apply(array(subX1), 1, function(x, subX)
    {
        (subX - x)^2
    }
    , subX=subX) + apply(array(subY1), 1, function(x, subY)
    {
        (subY - x)^2
    }
    , subY=subY))
    Cs1 <- apply(Cs, 2, min) <= MaxMigr
    Data[Data == 1][Cs1 != T] <- 0
    Data[Data == -1] <- 1
    Data[Data == -2] <- 0
    cat(".")
    return(Data)
}

