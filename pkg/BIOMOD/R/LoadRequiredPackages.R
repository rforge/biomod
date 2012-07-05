.LoadRequiredPackages <- function(biomod.material, raster.req = FALSE){
  algo <- biomod.material$algo[biomod.material$algo.choice]
  
 
  
  if("ANN" %in% algo){
    require(nnet, quietly=TRUE)
  }

  if("CTA" %in% algo){
    require(rpart, quietly=TRUE)
  }

  if("GAM" %in% algo){
    require(gam, quietly=TRUE)
    require(MASS, quietly=TRUE)
  }
  
  if("GLM" %in% algo){
    require(MASS, quietly=TRUE)
  }

  if("GBM" %in% algo){
    require(gbm, quietly=TRUE)
  }

  if("MARS" %in% algo | "FDA" %in% algo){
    require(mda, quietly=TRUE)
  }
  
  if("RF" %in% algo){
    require(randomForest, quietly=TRUE)
  }
  
  if(raster.req){
    require(foreign, quietly=TRUE)
    require(sp, quietly=TRUE)
    require(rgdal, quietly=TRUE)
    require(raster, quietly=TRUE)
    require(maptools, quietly=TRUE)
  }

}