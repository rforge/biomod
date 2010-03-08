`Rescaler3` <-
function(dataToRescale, OriMinMax=NULL)
{
	  if(is.null(OriMinMax)) OriMinMax <- range(dataToRescale)	
    RescaledData <- (dataToRescale - min(OriMinMax)) / (max(OriMinMax)-min(OriMinMax))
    
    RescaledData[which(RescaledData>1)] <- 1
	  RescaledData[which(RescaledData<0)] <- 0
	 
    return(RescaledData)
}
