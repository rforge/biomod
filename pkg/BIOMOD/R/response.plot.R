`response.plot` <-
function(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE){
  
  if(inherits(Data,"Raster")){
    cat("\n   > Extracting raster infos..")
    DataTmp <- matrix(0,ncol=nlayers(Data), nrow=100)
    colnames(DataTmp) <- layerNames(Data)
    maxVal <- maxValue(Data)
    minVal <- minValue(Data)
    for(i in 1:ncol(DataTmp)){
      DataTmp[,i] <- seq(minVal[i],maxVal[i],length.out=100)
    }
    Data <- DataTmp
    rm(list=c('maxVal','minVal','DataTmp'))
    
  }

    if(sum(show.variables > ncol(Data)) > 0) stop("columns wanted in show.variables do not match the data \n")

    NbVar <- ncol(Data)
    NbVarShow <- length(show.variables)
    
    if(plot==F) temp <- array(0, dim=c(nrow(Data), 2, NbVarShow), dimnames=list(NULL, c("Var", "Pred"), colnames(Data)[show.variables]))
    
    
    #consider if factorial variables :     
    Xp  <- as.data.frame(matrix(NA, ncol=NbVar, nrow=nrow(Data), dimnames=list(NULL, colnames(Data))))
    for(i in 1:NbVar){
        if(is.numeric(Data[,i])) { Xp[,i] <- mean(Data[,i])
        } 
        else { 
          Xp[, i] <- as.factor(rep(names(which.max(summary(Data[, i]))), nrow(Data)))
        	levels(Xp[,i]) <- levels(Data[, i])	 	
        }
    }   
      
    if(substr(class(model)[1],1,4)=="nnet" ) if(sum(search()=="package:nnet")==0) library(nnet)
    if(class(model)[1]=="rpart") if(sum(search()=="package:rpart")==0) library(rpart)
    if(class(model)[1]=="mars" | class(model)[1]=="fda") if(sum(search()=="package:mda")==0) library(mda)
    if("randomForest" %in% class(model)) if(sum(search()=="package:randomForest")==0) library(randomForest,  verbose=FALSE)
    
    if(plot){
    if(save.file=="pdf") pdf(paste(name, "pdf", sep="."))
    if(save.file=="jpeg") jpeg(paste(name, "jpeg", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="tiff") tiff(paste(name, "tiff", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="postscript") postscript(paste(name, "eps", sep="."))
    
    #plotting window
    W.width <- ceiling(sqrt(NbVarShow))
    W.height <- ceiling(NbVarShow/W.width)
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
    
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves ", class(model)[1], sep=""),col="#4c57eb")
    par(mar = c(2,2,3.5,1))

	}
    for(i in 1:NbVar){ if(sum(i==show.variables) > 0){
    
            #consider if factorial variables :
            if(!is.factor(Data[,i])){  
                xr <- range(Data[,i])
                Xp1 <- Xp
                Xp1[,i] <- seq(xr[1], xr[2],  len=nrow(Data))
            } else {
                Xp1 <- Xp
                Nrepcat <- floor(nrow(Data)/length(levels(Data[,i])))
                Xp1[,i] <- as.factor(c(rep(levels(Data[,i])[1], nrow(Data)-(Nrepcat*length(levels(Data[,i])))), rep(levels(Data[,i]), each=Nrepcat)))
        
            }
        
            if(class(model)[1]=="glm" | class(model)[1]=="gam") Xf <- predict(model, as.data.frame(Xp1), type="response")
            if(class(model)[1]=="gbm") Xf <-  predict.gbm(model, as.data.frame(Xp1), model$n.trees, type="response") 
            if(class(model)[1]=="rpart") Xf <- as.numeric(predict(model, Xp1, type="vector"))
            if(substr(class(model)[1],1,4)=="nnet" ) Xf <- as.numeric(predict(model, as.data.frame(Xp1), type="raw"))
            if(class(model)[1]=="mars") Xf <- as.numeric(predict(model, as.data.frame(Xp1)))
            if(class(model)[1]=="fda") Xf <- predict(model, as.data.frame(Xp1), type="post")[,2]
            if("randomForest" %in% class(model)) Xf <- predict(model, as.data.frame(Xp1), type="prob")[,2]
      
      
            #rescaling preds (not possible to use rescaling_GLM -> no info on calib data)
            if(class(model)[1]=="mars" | substr(class(model)[1],1,4)=="nnet"  | class(model)[1]=="fda" ){ 
                OriMinMax <- range(Xf)	
                Xf <- (Xf - min(OriMinMax)) / (max(OriMinMax)-min(OriMinMax))
                Xf[Xf<0]<-0
                Xf[Xf>1]<-1            
            }
#             cat("no")
			if(plot) {
				plot(Xp1[ ,i], Xf, ylim=c(0,1), xlab="", ylab="", type="l", main=names(Data)[i])
				rug(Data[ ,i])
			}	
			else{ 
				temp[,1,i] <-Xp1[ ,i]; temp[,2,i] <- Xf
			}     
    }}# i loop for variables
   
   
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
    if(plot==F) return(temp)
   
  #  if(substr(class(model)[1],1,4)=="nnet" )  detach(package:nnet)
   # if(class(model)[1]=="rpart") detach(package:rpart)
   # if(class(model)[1]=="mars" | class(model)[1]=="fda") detach(package:mda)
   # if(class(model)[1]=="randomForest") detach(package:randomForest)            
} 

 