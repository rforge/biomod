`Rescaler4` <-
function(dataToRescale, ref=NULL, run, original=FALSE)
{
    #preparing data
    #homogenize the format accross original predictions and new projections 
    if(!class(dataToRescale)[1]=='RasterLayer'){
        DataF <- as.data.frame(dataToRescale)         
        colnames(DataF) <- "DataF"                  
    } else{
        layerNames(dataToRescale) <-"DataF"
	      DataF <- stack(dataToRescale) 
    }
    
    #Creating or loading the rescaling model
    if(original){ 
        Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial")
        eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(), "/models/rescaling_models/Rmod_", run, "')", sep=""))) 
    } else
        eval(parse(text=paste("load('", getwd(), "/models/rescaling_models/Rmod_", run, "')", sep="")))
	 	
    #make the rescaling prediction
    if(!class(dataToRescale)[1]=='RasterLayer') RescaledData <- predict(Rescaling_GLM, DataF, type="response") 
	  if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
	   
    
    #Just in case
    RescaledData[which(RescaledData>1)] <- 1
	  RescaledData[which(RescaledData<0)] <- 0
       
    return(RescaledData)
}
