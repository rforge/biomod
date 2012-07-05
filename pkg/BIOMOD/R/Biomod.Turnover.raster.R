`Biomod.Turnover.raster` <-
function(CurrentPred=NULL, FutureProj=NULL, Turnover.Save=NULL)
{
    #Function to estimate number of species loss, gain and stable per pixel
    
    if(class(CurrentPred)=="RasterLayer") CurrentPred <- stack(CurrentPred) 
    if(class(FutureProj)=="RasterLayer") FutureProj <- stack(FutureProj)
    if(length(FutureProj@layers) != length(CurrentPred@layers)) stop("\n The two inputs (CurrentPred and FutureProj) are of different length")

    
    #Create blank raster
    BLANK.ras <- CurrentPred@layers[[1]]
    BLANK.ras[!is.na(BLANK.ras)] <- 0



    CurrentSR <- sum(CurrentPred)
    FutureSR.FullDisp <- sum(FutureProj)

    Gain <- Loss <- CurrentSR  -  FutureSR.FullDisp
    Loss[Loss<0] <- 0
    Gain[Gain>0] <- 0
    Gain <- abs(Gain)
    FutureSR.NoDisp <- CurrentSR - Loss
        
    #Stable 0s and 1s    
    Stable1 <- Stable0 <- BLANK.ras
    for(i in 1:length(CurrentPred@layers)){
        Cur <- CurrentPred@layers[[i]]
        Fut <- FutureProj@layers[[i]]
        Stable0 <- Stable0 + (Cur==0 & Fut==0) 
        Stable1 <- Stable1 + (Cur==1 & Fut==1)
    }
    
    #percentages
    PercLoss <- Loss / CurrentSR *100
    PercGain <- Gain / CurrentSR *100
    Turnover <- (Loss + Gain) / (CurrentSR + Gain) *100    
    
    
       
    StacK <- stack(Loss, Stable0, Stable1, Gain, PercLoss, PercGain, Turnover, CurrentSR, FutureSR.NoDisp, FutureSR.FullDisp)
    layerNames(StacK) <- c("Loss", "Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "Turnover", "CurrentSR", "FutureSR.NoDisp", "FutureSR.FullDisp")
    
    if(is.null(Turnover.Save)) Turnover.Save <- "NoName_Turnover"
    assign(Turnover.Save, StacK, pos=1)
}

