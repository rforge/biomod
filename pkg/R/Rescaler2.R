`Rescaler2` <-
function(DataToRescale, type="range", OriMinMax=NULL)
{
	if(!is.null(OriMinMax)){
		#cat("Safe MODE was used: Original MinMax: ", OriMinMax, "  To Rescale MinMax: ", min(DataToRescale), max(DataToRescale), "\n")
		# Warn the user if he gets outside of the calibration range.
		MinMaxRange <- abs(OriMinMax[2] - OriMinMax[1])
		if(min(DataToRescale) < min(DataToRescale) - MinMaxRange*0.1) warning("Prediction Range is more than 10% smaller than calibration range: Prediction quality might be affected! \n")
		if(max(DataToRescale) > max(DataToRescale) + MinMaxRange*0.1) warning("Prediction Range is more than 10% larger than calibration range: Prediction quality might be affected! \n")
		DataToRescale <- c(DataToRescale, OriMinMax)
		rm(MinMaxRange)
	}
	
	RescaledData <- DataToRescale
	if(length(unique(DataToRescale))==1){
		if(DataToRescale[1]>1) RescaledData[1:length(DataToRescale)] <- 1
		if(DataToRescale[1]<0) RescaledData[1:length(DataToRescale)] <- 0
	}
	else RescaledData <- rescaler(DataToRescale, type=type)
	
	# remove the two values we added at the end of the vector:
	if(!is.null(OriMinMax)) RescaledData <- RescaledData[-c(length(DataToRescale)-1, length(DataToRescale))]
	
	return(RescaledData)
}

