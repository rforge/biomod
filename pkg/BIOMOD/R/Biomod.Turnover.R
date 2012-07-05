`Biomod.Turnover` <-
function(CurrentPred=NULL, FutureProj=NULL, Turnover.Save=NULL)
{
    #handle exceptions here
    #------------------------------------------------------------------------------------
    #Function to estimate number of species loss, gain and stable per pixel
    Value <- c(-2, 0, -1, 1)
    Compteur <- function(Data, Value)
    {
        N <- dim(Data)[2]
        Compteur <- as.data.frame(matrix(0, ncol=4, nrow=dim(Data)[1], dimnames=list(seq(dim(Data)[1]), c("disapearing", "stable 0", "stable 1", "gaining"))))
        i <- 1
        while(i <= N) {
            #cat(i,"\t")
            Compteur[Data[,i] == Value[1], 1] <- Compteur[Data[,i] == Value[1], 1] + 1
            Compteur[Data[,i] == Value[2], 2] <- Compteur[Data[,i] == Value[2], 2] + 1
            Compteur[Data[,i] == Value[3], 3] <- Compteur[Data[,i] == Value[3], 3] + 1
            Compteur[Data[,i] == Value[4], 4] <- Compteur[Data[,i] == Value[4], 4] + 1
            i <- i + 1
        }
        return(Compteur)
    }
    Diff.By.Pixel <- as.data.frame(FutureProj - 2 * CurrentPred)
    Compt.By.Pixel <- as.data.frame(Compteur(Diff.By.Pixel, Value))
    SR.Current <- (Compt.By.Pixel[, 1] + Compt.By.Pixel[,3])
    Compt.By.Pixel[,5] <- (100 * Compt.By.Pixel[, 1])/SR.Current
    Compt.By.Pixel[,6] <- (100 * Compt.By.Pixel[, 4])/SR.Current
    Compt.By.Pixel[,7] <- (100 * (Compt.By.Pixel[, 1] + Compt.By.Pixel[, 4]))/(SR.Current + Compt.By.Pixel[, 4])
    Compt.By.Pixel[,8] <- SR.Current
    Compt.By.Pixel[,9] <- Compt.By.Pixel[,3]
    Compt.By.Pixel[,10] <- (Compt.By.Pixel[,3] + Compt.By.Pixel[, 4])
    dimnames(Compt.By.Pixel) <- list(seq(1:dim(Compt.By.Pixel)[1]), c(
        "Loss", "Stable0", "Stable1", "Gain", "PercLoss", "PercGain",
        "Turnover", "CurrentSR", "FutureSR.NoDisp", "FutureSR.FullDisp"))
  assign(Turnover.Save, Compt.By.Pixel, pos=1)
    invisible(Compt.By.Pixel)
}

