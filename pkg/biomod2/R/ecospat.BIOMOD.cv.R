##' @name ecospat.BIOMOD.cv
##' @aliases ecospat.BIOMOD.cv
##' 
##' @title Custom models cross-validation procedure
##' 
##' @description This function creates a DataSplitTable which could be used to evaluate models in Biomod with repeated
##'   k-fold cross-validation (cv) or stratified cv instead of repeated split sample runs
##' 
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param k               number of bins/partitions for k-fold cv
##' @param stratified.cv   logical. run a stratified cv 
##' @param stratify        stratification method of the cv. Could be "x", "y", "both" (default), "block" or the name of a predictor for environmental stratified cv.
##' @param balance         make balanced particions for "presences" (default) or "absences" (resp. pseudo-absences or background).
##' @param repetition      number of repetitions of k-fold cv (1 if stratified.cv=TRUE)
##' @param do.full.models  if true, models calibrated and evaluated with the whole dataset are done
##' 
##' @details
##'   Stratified cv could be used to test for model overfitting and for assessing transferability in geographic and environmental space. 
##'   If balance = "presences" presences are divided (balanced) equally over the particions (e.g. Fig. 1b in Muscarelly et al. 2014).
##'   Pseudo-Absences will however be unbalanced over the particions especially if the presences are clumped on an edge of the study area.
##'   If balance = "absences" absences (resp. Pseudo-Absences or background) are divided (balanced) as equally as possible for the particions
##'   (geographical balanced bins given that absences are spread over the study area equally, approach similar to Fig. 1 in Wenger et Olden 2012).
##'   Presences will however be unbalanced over the particians. Be careful: If the presences are clumped on an edge of the study area it is possible that all presences are in one bin.
##' 
##' @return
##' DataSplitTable matrix with k*repetition (+ 1 for Full models if  do.full.models = TRUE) columns for BIOMOD_Modeling function.
##' Stratification "x" and "y" was described in Wenger and Olden 2012. While Stratification "y" uses k partitions along the y-gradient, "x" does the same for the x-gradient and "both" combines them.
##' Stratification "block" was described in Muscarella et al. 2014. For bins of equal number are partitioned (bottom-left, bottom-right, top-left and top-right).
##' 
##' @author Frank Breiner
##' 

ecospat.BIOMOD.cv <- function(data, k=5,repetition=5, do.full.models = TRUE, stratified.cv=FALSE, stratify="both", balance="pres"){
 
  DataSplitTable.y <-  DataSplitTable.x <-   DataSplitTable <- NULL
  if(stratified.cv){
    repetition <- 1
    if(balance == "absences"){balance <- data@data.species==1| data@data.species==0
    }else{balance <- data@data.species==1}
  if(stratify == "x" | stratify == "both"){
    DataSplitTable.x <- matrix(NA,nrow(data@coord),k)
    bands <- quantile(data@coord[balance,1],  probs = seq(0,100,100/k)/100)
    bands[1] <- bands[1]-1
    bands[k+1] <- bands[k+1]+1
    for(i in 1:k){
    DataSplitTable.x[,i] <- data@coord[,1] >= bands[i] & data@coord[,1] < bands[i+1]
  }
    if(stratify == "x"){DataSplitTable <- DataSplitTable.x}
  }
  if(stratify == "y" | stratify == "both"){
    DataSplitTable.y <- matrix(NA,nrow(data@coord),k)
    bands <- quantile(data@coord[balance,2],  probs = seq(0,100,100/k)/100)
    bands[1] <- bands[1]-1
    bands[k+1] <- bands[k+1]+1
    for(i in 1:k){
      DataSplitTable.y[,i] <- data@coord[,2] >= bands[i] & data@coord[,2] < bands[i+1]
    }
      if(stratify == "y"){DataSplitTable <- DataSplitTable.y}
  }
  if(stratify == "both"){
    DataSplitTable  <- cbind(DataSplitTable.x, DataSplitTable.y)
  }
  if(stratify == "block"){
    DataSplitTable <- matrix(NA,nrow(data@coord),4)
    DataSplitTable[,1] <- data@coord[,1] < median(data@coord[balance,1]) & data@coord[,2] < median(data@coord[balance,2]) # bottom-left
    DataSplitTable[,2] <- data@coord[,1] < median(data@coord[balance,1]) & data@coord[,2] > median(data@coord[balance,2]) # top-left
    DataSplitTable[,3] <- data@coord[,1] > median(data@coord[balance,1]) & data@coord[,2] < median(data@coord[balance,2]) # bottom-right
    DataSplitTable[,4] <- data@coord[,1] > median(data@coord[balance,1]) & data@coord[,2] > median(data@coord[balance,2]) # top-rigth
  }
  if(stratify != "block" & stratify != "x" & stratify != "y" & stratify != "both"){
    DataSplitTable2 <- matrix(NA,nrow(data@coord),k)
    bands <- quantile(data@data.env.var[balance,stratify],  probs = seq(0,100,100/k)/100)
    bands[1] <- bands[1]-1
    bands[k+1] <- bands[k+1]+1
    for(i in 1:k){
      DataSplitTable2[,i] <- data@data.env.var[balance,stratify] <= bands[i] | data@data.env.var[balance,stratify] > bands[i+1]
    }
  }
  }else{
  for(rep in 1:repetition){ 
   fold <- dismo::kfold(data@data.species,by=data@data.species,k=k) 
  for(i in 1:k){
    DataSplitTable <- cbind(DataSplitTable,fold!=i)
  }
  }
  }
  colnames(DataSplitTable) <- paste("RUN",1:(k*repetition),sep="")
  if(do.full.models == TRUE){
    DataSplitTable <- cbind(DataSplitTable,T)
    colnames(DataSplitTable)[k*repetition+1] <- "Full"
  }
  return(DataSplitTable)  
}
  