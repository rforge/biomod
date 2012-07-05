.Prepare.Maxent.WorkDir <- function(Data, xy, calibLines, RunName=NULL, VarImport=0, evalData=NULL, evalxy=NULL ){
  cat('\n\tCreating Maxent Temp Proj Data..')
  dir.create(paste(getwd(),'/MaxentTmpData', sep=""), showWarnings=FALSE)
  dir.create(paste(getwd(),'/',colnames(Data)[1],'/models/',RunName,'_MAXENT',sep=''), showWarnings=FALSE)
  
  if(is.null(RunName)) RunName <- colnames(Data)[1]
  # Presences Data
  presLines <- which((Data[,1]==1) & calibLines)
  absLines <- which((Data[,1]==0) & calibLines)
  Sp_swd <- cbind(rep(RunName,length(presLines)),
                      xy[presLines,],
                      Data[presLines,2:ncol(Data)])
  colnames(Sp_swd) <- c('specie','X','Y',colnames(Data)[2:ncol(Data)])
  write.table(Sp_swd, file=paste(getwd(),"/MaxentTmpData/Sp_swd.csv",sep=""), quote=FALSE, row.names=FALSE, sep=",")
  
  # Background Data
  # keep only 0 of calib lines
  Back_swd <- cbind(rep("background",length(absLines)),xy[absLines,],Data[absLines,2:ncol(Data)])
  colnames(Back_swd)  <- c("background",colnames(Back_swd)[-1])
  write.table(Back_swd, file=paste(getwd(),"/MaxentTmpData/Back_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  # Prediction Data
  dir.create(paste(getwd(),'/MaxentTmpData/Pred', sep=""), showWarnings=FALSE)
  
  Pred_swd <- cbind(rep("predict",nrow(xy)),xy,Data[,2:ncol(Data)])
  colnames(Pred_swd)  <- c("predict",colnames(Back_swd)[-1])
  write.table(Pred_swd, file=paste(getwd(),"/MaxentTmpData/Pred/Pred_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  # dealing with variable importances stuff
  if( VarImport > 0){
    for( vari in colnames(Data)[-1] )
      for(vi in 1:VarImport){
        proj_tmp <- Pred_swd
        proj_tmp[,1] <- rep(paste(vari,'_',vi,sep=""),nrow(proj_tmp))
        proj_tmp[,vari] <- sample(proj_tmp[,vari])
        write.table(proj_tmp, file=paste(getwd(),"/MaxentTmpData/Pred/",vari,'_',vi,"_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
      }
  }
  
  # dealing with independent evaluation data
  if(!is.null(evalData)){
    Pred_eval_swd <- cbind(rep("predictEval",nrow(evalxy)),evalxy,evalData[,2:ncol(evalData)])
    colnames(Pred_eval_swd)  <- c("predict",colnames(Back_swd)[-1])
    write.table(Pred_eval_swd, file=paste(getwd(),"/MaxentTmpData/Pred/Pred_eval_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  }
}

.Delete.Maxent.WorkDir <- function(){
  cat('\n\tRemoving Maxent Temp Data..')
#   system('rm -rf MaxentTmpData')
  unlink('MaxentTmpData', recursive = TRUE, force = TRUE)
}

.Prepare.Maxent.Proj.WorkDir <- function(Data, xy){
  cat('\n\tCreating Maxent Temp Proj Data..')
  dir.create(paste(getwd(),'/MaxentTmpData', sep=""), showWarnings=FALSE)
  
  # Proj Data
  Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
  colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
  write.table(Proj_swd, file=paste(getwd(),"/MaxentTmpData/Proj_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
}

.Prepare.Maxent.Proj.Raster.WorkDir <- function(Data){
  cat('\n\tCreating Maxent Temp Proj Data..')
  dir.create(paste(getwd(),'/MaxentTmpData/Proj', sep=""), showWarnings=FALSE, recursive=TRUE)
  
  # Proj Data
  for(i in 1:nlayers(Data)){
    writeRaster(Data[[i]], filename=paste(getwd(),'/MaxentTmpData/Proj/',layerNames(Data)[i],'.asc',sep=''),
                format='ascii', overwrite=TRUE)
  }
}
