# `.Rescaler4` <-
# function(dataToRescale, ref=NULL, run, original=FALSE)
# {
#     #preparing data
#     #homogenize the format accross original predictions and new projections 
#     if(!class(dataToRescale)[1]=='RasterLayer'){
#         DataF <- as.data.frame(dataToRescale)         
#         colnames(DataF) <- "DataF"                  
#     } else{
#         layerNames(dataToRescale) <-"DataF"
# 	      DataF <- stack(dataToRescale) 
#     }
#     
#     #Creating or loading the rescaling model
#     if(original){ 
#         Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial")
#         eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(), "/models/rescaling_models/Rmod_", run, "', compress='xz')", sep=""))) 
#     } else
#         eval(parse(text=paste("load('", getwd(), "/models/rescaling_models/Rmod_", run, "')", sep="")))
# 	 	
#     #make the rescaling prediction
#     if(!class(dataToRescale)[1]=='RasterLayer') RescaledData <- predict(Rescaling_GLM, DataF, type="response") 
# 	if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
# 	   
#    
#     return(RescaledData)
# }

.Rescaler5 <-
function(dataToRescale, ref=NULL, name, original=FALSE)
{
#     #preparing data
#     #homogenize the format accross original predictions and new projections 
#     if(!class(dataToRescale)[1]=='RasterLayer'){
        DataF <- as.data.frame(dataToRescale)         
        colnames(DataF) <- "DataF"                  
#     } else{
#         layerNames(dataToRescale) <-"DataF"
#         DataF <- stack(dataToRescale) 
#     }
#     
    #Creating or loading the rescaling model
  if(original){
      if(! file.exists(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/rescaling_models/", sep=""))){
        dir.create(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/rescaling_models/", sep=""), showWarnings=F)
      }
      Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial", mustart = rep(0.5,length(ref)))
      eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(),"/",
                            unlist(strsplit(name,'_'))[1], "/models/rescaling_models/",
                            name, "_rescaled' , compress='",ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
                            ,"')", sep=""))) 
    } else{
      eval(parse(text=paste("load('", getwd(),"/",unlist(strsplit(name,'_'))[1],
                            "/models/rescaling_models/",name,"_rescaled')", sep="")))
    }
    #make the rescaling prediction
    RescaledData <- predict(Rescaling_GLM, DataF, type="response") 
# 	  if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
	   
    return(RescaledData)
}