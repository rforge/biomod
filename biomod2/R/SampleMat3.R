SampleMat3 <- function(pres, abs, dataSplit, nbRun = 1){
  # return a matrix with the rowId in a first col and '1' (calib) or '2' (eval) for nbRun
  pres <- as.vector(pres)
  abs <- as.vector(abs)
  
  nbPresEval <- round(length(pres) * dataSplit/100)
  nbAbsEval <- round(length(abs) * dataSplit/100)
  
#   cat ("nbPres = ", nbPresEval,"\tnbAbs = ", nbAbsEval)
  
  mat.out <- matrix(2, nrow=length(pres)+length(abs), ncol=nbRun+1)
  colnames(mat.out) <- c('Ids',paste('Run',1:nbRun,sep=''))
  
  mat.out[,1] = c(pres,abs)
  
  for (r in 2:ncol(mat.out)){
    mat.out[sample(1:length(pres),nbPresEval),r] <- 1
    mat.out[sample((length(pres)+1) :nrow(mat.out) ,nbAbsEval),r] <- 1
  }
  
  return(mat.out[order(mat.out[,1]),])

}