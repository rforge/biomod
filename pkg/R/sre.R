`sre` <- function(Response=NULL, Explanatory=NULL, NewData=NULL, Quant=0.025, na.rm=TRUE)
{
    if(Quant>=0.5 | Quant<0) stop("\n settings in Quant should be a value between 0 and 0.5 ")
    quants <- c(0+Quant, 1-Quant)
    
    
    #Check the variables from calib to predict
    #%%%%%%%%%%%
    
    
    #Calibration on point data, projection on points or rasters
    if(class(Explanatory)[1]!='RasterStack'){
        Response <- as.data.frame(Response)                                                        #to get the colnames if only a vector is given
        NbVar <- dim(Explanatory)[2]
        #storing output
        if(class(NewData)[1]!='RasterStack') Pred <- as.data.frame(matrix(0, nr=dim(NewData)[1], nc=dim(Response)[2], dimnames=list(seq(dim(NewData)[1]), colnames(Response))))
        
        #for multiple species(not active for rasters for now -> overwritting TF)
        for(i in 1:dim(Response)[2]){
            ref <- Explanatory[Response[,i]==1,]
            
            #object for storing the True/False values (now +1s)
            if(class(NewData)[1]=='RasterStack') { TF <- NewData@layers[[1]] ; TF <- TF<TF@data@min }                 #set raster to FALSE for all values (exept NAs)
            else TF <- rep(0, dim(NewData)[1])
            
            
            for(j in 1:NbVar){
                Q <- quantile(ref[,j], probs=quants, na.rm=na.rm)
                
                if(class(NewData)[1]!='RasterStack'){
                    TF <- TF + (NewData[,names(ref)[j]]>=Q[1])                                          #add the T/F values to vector (+T = +1)
                    TF <- TF + (NewData[,names(ref)[j]]<=Q[2])
                } else{                                                                             #NewData = raster
                    TFmin <- NewData@layers[[which(NewData@layernames==names(ref)[j])]]>=Q[1]       #produce T/F rasterlayers
                    TFmax <- NewData@layers[[which(NewData@layernames==names(ref)[j])]]<=Q[2]       
                    TFmin["TRUE"] <- 1                                                              #convert them to 0s and 1s
                    TFmax["TRUE"] <- 1
                    TF <- TF + TFmin + TFmax                                                        #add to TF for storage accross variables
                }      
            }
            TF[TF!=(NbVar*2)] <- 0                                                                  #convert to binary
            TF[TF==(NbVar*2)] <- 1                                                                  #important to set to 0 first, then the 1s
            if(class(TF)[1]!='RasterLayer') Pred[,i] <- TF else Pred <- TF                                          #store in matrix if point data
        }
    }
    
    #Calibration on RasterStack, projection on rasters  
    #else{ 
        
    #    Q <- quantile(belalp.stk@layers[[1]], na.rm=T, probs=quants)
        
    #}
            
            

    
    if(class(NewData)[1]!='RasterStack' & dim(Response)[2]==1) Pred <- Pred[[1]]   #return a vector, not a data frame
    return(Pred)
}
