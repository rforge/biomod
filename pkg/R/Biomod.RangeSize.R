`Biomod.RangeSize` <-
function(CurrentPred=NULL, FutureProj=NULL, 
SpChange.Save=NULL)
{
    #Function that estimate the number of pixel gain, loss, stable (present and absent) by species
    Value <- c(-2, 0, -1, 1)
    CompteurSp <- function(Data, Value)
    {
        if(is.data.frame(Data)) {
            N <- dim(Data)[2]
            Compt <- as.data.frame(matrix(0, ncol=4, nrow=dim(Data)[2]))
            i <- 1
            while(i <= N) {
                Compt[i, 1] <- length(Data[Data[,i] == Value[1], i])
                Compt[i, 2] <- length(Data[Data[,i] == Value[2], i])
                Compt[i, 3] <- length(Data[Data[,i] == Value[3], i])
                Compt[i, 4] <- length(Data[Data[,i] == Value[4], i])
                i <- i + 1
            }
        }
        return(Compt)
    }
    Diff.By.Pixel <- as.data.frame(FutureProj - 2 * CurrentPred)
    Compt.By.Species <- as.data.frame(CompteurSp(Diff.By.Pixel, Value))
    Compt.By.Species[, 5] <- (100 * Compt.By.Species[, 1])/(Compt.By.Species[, 1] + Compt.By.Species[,3])
    Compt.By.Species[, 6] <- (100 * Compt.By.Species[, 4])/(Compt.By.Species[, 1] + Compt.By.Species[,3])
    Compt.By.Species[, 7] <- Compt.By.Species[, 6] - Compt.By.Species[, 5]
    Compt.By.Species[, 8] <- Compt.By.Species[, 1] + Compt.By.Species[,3]
    Compt.By.Species[, 9] <- Compt.By.Species[,3]
    Compt.By.Species[, 10] <- Compt.By.Species[, 4] + Compt.By.Species[,3]
    dimnames(Compt.By.Species) <- list(colnames(CurrentPred), c("Disa","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", 
        "SpeciesRangeChange", "CurrentRangeSize", 
        "FutureRangeSize.0Disp", "FutureRangeSize.1Disp"))
    Output <- c()
    Output <- list(Compt.By.Species=Compt.By.Species, Diff.By.Pixel=Diff.By.Pixel)
    assign(SpChange.Save, Output, pos=1)
    invisible(Output)
}

