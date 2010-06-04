`Biomod.RangeSize.raster` <-
function(CurrentPred=NULL, FutureProj=NULL, SpChange.Save=NULL)
{
    #Function to estimate number of species loss, gain and stable per pixel

    if(class(CurrentPred)=="RasterLayer") CurrentPred <- stack(CurrentPred) 
    if(class(FutureProj)=="RasterLayer") FutureProj <- stack(FutureProj)
    if(length(FutureProj@layers) != length(CurrentPred@layers)) stop("\n The two inputs (CurrentPred and FutureProj) are of different length")

    CBS <- matrix(ncol=10, nrow=length(CurrentPred@layers), dimnames=list(layerNames(CurrentPred), 
      c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")))

    
    sp.stack <- stack()  
    for(i in 1:length(CurrentPred@layers)){
        #DiffByPixel
        Cur <- CurrentPred@layers[[i]]
        Fut <- FutureProj@layers[[i]]
        Ras <- Fut - 2 * Cur
        sp.stack <- stack(sp.stack, Ras)
        
        #ComptBySpecies
        CBS[i, 1] <- length(which(Ras[]==-2))
        CBS[i, 2] <- length(which(Ras[]== 0))
        CBS[i, 3] <- length(which(Ras[]==-1))
        CBS[i, 4] <- length(which(Ras[]== 1))
        
        CBS[i, 5] <- round(CBS[i,1] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
        CBS[i, 6] <- round(CBS[i,4] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
        CBS[i, 7] <- round((CBS[i,3]+CBS[i,4]) / (CBS[i,1]+CBS[i,3]) *100 -100, digits=3)       
        
        CBS[i, 8] <- CBS[i,1]+CBS[i,3]                                        
        CBS[i, 9] <- CBS[i,3]                                                 
        CBS[i, 10] <- CBS[i,3]+CBS[i,4]                                       
    }

    if(is.null(SpChange.Save)) SpChange.Save <- "NoName"
    assign(paste(SpChange.Save, "_Compt.By.Species", sep=""), CBS, pos=1)
    #layerNames(sp.stack) <- Biomod.material$species.names
    assign(paste(SpChange.Save, "_Diff.By.Pixel", sep=""), sp.stack, pos=1)
    return(CBS)
}

