# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD objects definition
# Damien Georges
# 09/02/2012
# v2.0
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This file defines the BIOMOD objects and all their methods 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# We choose here to create monospecific objects to make all ppROC:::rocedures and parallelising easier
require(sp, quietly=TRUE)
require(raster, quietly=TRUE)
require(rasterVis, quietly=TRUE)

# if('mgcv' %in% rownames(installed.packages())){
#   require(mgcv, quietly=TRUE)
# } else{
#   require(gam, quietly=TRUE)
# }
# 
# require(rpart, quietly=TRUE)
# require(mda, quietly=TRUE)

# 1. The BIOMOD.formated.data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this object is the basic one

# 1.1 Class Definition
setClass("BIOMOD.formated.data",
         representation(sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
#                         data.counting = "matrix",
                        data.env.var = "data.frame",
                        data.mask = "RasterStack",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ return(TRUE) } )

# 1.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data" ) ) {
  setGeneric( "BIOMOD.formated.data", 
              def = function(sp, env, ...){
  	                  standardGeneric( "BIOMOD.formated.data" )
                      } )
# }

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='data.frame' ), 
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE, data.mask=NULL ){
    if(is.null(data.mask)) data.mask <- stack()
    
    
    if(is.null(eval.sp)){
      BFD <- new('BIOMOD.formated.data', 
                 coord=xy, 
                 data.species=sp, 
                 data.env.var=env, 
                 sp.name=sp.name,
                 data.mask=data.mask,
                 has.data.eval=FALSE)
    } else{
      BFDeval <- BIOMOD.formated.data(sp=eval.sp,
                                      env=eval.env,
                                      xy=eval.xy,
                                      sp.name=sp.name)
      
      if(nlayers(BFDeval@data.mask)>0){
        data.mask.tmp <- try(addLayer(data.mask,BFDeval@data.mask))
        if( !inherits(data.mask.tmp,"try-error")){
          data.mask <- data.mask.tmp
          names(data.mask) <- c("calibration","validation")
        }
      }
      
      BFD <- new('BIOMOD.formated.data', 
                 coord=xy, 
                 data.species=sp, 
                 data.env.var=env, 
                 sp.name=sp.name,
                 data.mask=data.mask,
                 has.data.eval=TRUE,
                 eval.coord = BFDeval@coord,
                 eval.data.species = BFDeval@data.species,
                 eval.data.env.var = BFDeval@data.env.var )
      
                 
      rm('BFDeval')
    }
    if(na.rm){
      rowToRm <- unique(unlist(lapply(BFD@data.env.var,function(x){return(which(is.na(x)))})))
      if(length(rowToRm)){
        cat("\n\t\t\t! Some NAs have been automaticly removed from your data")
        BFD@coord <- BFD@coord[-rowToRm,]
        BFD@data.species <- BFD@data.species[-rowToRm]
        BFD@data.env.var <- BFD@data.env.var[-rowToRm,,drop=FALSE]
      }
      if(BFD@has.data.eval){
        rowToRm <- unique(unlist(lapply(BFD@eval.data.env.var,function(x){return(which(is.na(x)))})))
        if(length(rowToRm)){
          cat("\n\t\t\t! Some NAs have been automaticly removed from your evaluation data")
          BFD@eval.coord <- BFD@eval.coord[-rowToRm,]
          BFD@eval.data.species <- BFD@eval.data.species[-rowToRm]
          BFD@eval.data.env.var <- BFD@eval.data.env.var[-rowToRm,,drop=FALSE]
        }      
      }
    }
    
    
    
    # count data occutances
#     BFD@data.counting <- matrix(c(sum(BFD@data.species, na.rm=TRUE),sum(BFD@data.species==0, na.rm=TRUE)),
#                             ncol=1,nrow=2, dimnames=list(c("nb.pres","nb.abs"),c("data.species") ) )
#     
#     if(BFD@has.data.eval){
#       BFD@data.counting <- cbind(BFD@data.counting,c(sum(BFD@eval.data.species, na.rm=TRUE),sum(BFD@eval.data.species==0, na.rm=TRUE)))
#       colnames(BFD@data.counting)[ncol(BFD@data.counting)] <- "eval.data.species"
#     }
    
    return(BFD)
	}
)

setMethod('BIOMOD.formated.data', signature(sp='data.frame'), 
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
    if(ncol(sp) > 1 ){
      stop("Invalid response variable")
    }
    sp <- as.numeric(unlist(sp))
    BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
    return(BFD)
  }
)

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='matrix' ), 
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
    env <- as.data.frame(env)
    BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
    return(BFD)
  }
)
          
setMethod('BIOMOD.formated.data', signature(sp='numeric', env='RasterStack' ), 
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
    categorial_var <- names(env)[raster:::is.factor(env)]
    
    # take the same eval environemental variables than calibrating ones 
    if(!is.null(eval.sp)){
      if(is.null(eval.env)){
        eval.env <- as.data.frame(extract(env,eval.xy))
        if(length(categorial_var)){
          for(cat_var in categorial_var){
            eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
          }
        }
      }
    }
    
    if(is.null(xy)) xy <- as.data.frame(coordinates(env))
      
    data.mask = reclassify(raster:::subset(env,1,drop=T), c(-Inf,Inf,-1))
    data.mask[cellFromXY(data.mask,xy[which(sp==1),])] <- 1
    data.mask[cellFromXY(data.mask,xy[which(sp==0),])] <- 0
    data.mask <- stack(data.mask)
    names(data.mask) <- sp.name
    
    env <- as.data.frame(extract(env,xy))
    
    if(length(categorial_var)){
      for(cat_var in categorial_var){
        env[,cat_var] <- as.factor(env[,cat_var])
      }
    }
    
    BFD <- BIOMOD.formated.data(sp,env,xy,sp.name,eval.sp, eval.env, eval.xy, na.rm=na.rm, data.mask=data.mask)
    .bmCat('Done')
    return(BFD)
  }
)

# 1.3 Other Functions
if( !isGeneric( "plot" ) ) {
  setGeneric( "plot", 
              def = function(x, ...){
  	                  standardGeneric( "plot" )
                      } )
}

setMethod('plot', signature(x='BIOMOD.formated.data'),
          function(x,coord=NULL,col=NULL){
            if(nlayers(x@data.mask)>0){
              
              ## define the breaks of the color key
              my.at <- seq(-1.5,1.5,by=1)
              ## the labels will be placed vertically centered
              my.labs.at <- seq(-1,1,by=1)
              ## define the labels
              my.lab = c("undifined","absences","presences")
              
              levelplot(x@data.mask, at=my.at, margin=T, col.regions=c("lightgrey","red4","green4"),
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))
              
            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'grey')
              }
              
              # plot data
              # all points (~mask) 
              
              plot(x=x@coord[,1], y=x@coord[,2], col=col[3], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name, sep=""), pch=20 )
              # presences 
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)
              
            }
            

        
          })

setMethod('show', signature('BIOMOD.formated.data'),
          function(object){
            .bmCat("'BIOMOD.formated.data'")
            cat("\nsp.name = ", object@sp.name, fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            .bmCat()
          })

# 2. The BIOMOD.formated.data.PA =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this class inherits from BIOMOD.formated.data and have one more slot 'PA' giving PA selected

# 2.1 Class Definition
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy='character', PA = 'data.frame'),
         validity = function(object){
           return(TRUE)
           })

# 2.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data.PA" ) ){
#   setGeneric( "BIOMOD.formated.data.PA", 
#               def = function(sp, env, PA.NbRep, ...){
#                 standardGeneric( "BIOMOD.formated.data.PA" )
#               })
# }

# setMethod('BIOMOD.formated.data.PA',signature(sp='ANY',
#                                               env='ANY',
#                                               PA.NbRep='integer'),

BIOMOD.formated.data.PA <-  function(sp, env, xy, sp.name,
                                     eval.sp=NULL, eval.env=NULL, eval.xy=NULL,
                                     PA.NbRep=1,
                                     PA.strategy='random',
                                     PA.nb.absences = NULL,
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     PA.table = NULL,
                                     na.rm=TRUE){
  
  if(inherits(env,'Raster')){
    categorial_var <- names(env)[raster:::is.factor(env)]
  }  else categorial_var <- NULL
  
  
  # take the same eval environemental variables than calibrating ones 
  if(!is.null(eval.sp)){
    if(is.null(eval.env)){
      if(inherits(env,'Raster')){
        eval.env <- as.data.frame(extract(env,eval.xy))
        if(length(categorial_var)){
          for(cat_var in categorial_var){
            eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
          }
        }
      } else{
        stop("No evaluation explanatory variable given")
      }
    }
  }
                                     
  # convert sp in spatial obj
  if(is.numeric(sp)){
    if(is.null(xy)){
      sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp),match.ID=FALSE)
    } else{
      sp <- SpatialPointsDataFrame(data.matrix(xy), data.frame(sp),match.ID=FALSE)
    }
    
  }

  pa.data.tmp <- .pseudo.absences.sampling(sp = sp,
                                          env = env,
                                          nb.repet = PA.NbRep,
                                          strategy = PA.strategy,
                                          nb.points = PA.nb.absences,
                                          distMin = PA.dist.min, 
                                          distMax = PA.dist.max,
                                          quant.SRE = PA.sre.quant,
                                          PA.table = PA.table)
    
  if(!is.null(pa.data.tmp)){
    
    if(length(categorial_var)){
      for(cat_var in categorial_var){
        pa.data.tmp$env[,cat_var] <- as.factor(pa.data.tmp$env[,cat_var])
      }
    }
    
    
    if(na.rm){
      rowToRm <- unique(unlist(lapply(pa.data.tmp$env,function(x){return(which(is.na(x)))})))
      if(length(rowToRm)){
        cat("\n\t\t\t! Some NAs have been automaticly removed from your data")
        pa.data.tmp$xy <- pa.data.tmp$xy[-rowToRm,]
        pa.data.tmp$sp <- pa.data.tmp$sp[-rowToRm]
        pa.data.tmp$env <- pa.data.tmp$env[-rowToRm,,drop=FALSE]
        pa.data.tmp$pa.tab <- pa.data.tmp$pa.tab[-rowToRm,,drop=FALSE]
      }      
    }
    
    
    # data counting
#     pa.data.tmp$data.counting <- apply(pa.data.tmp$pa.tab,2,function(x){nbPres <- sum(pa.data.tmp$sp[x],na.rm=T) ; return(c(nbPres,sum(x)-nbPres))})
#     colnames(pa.data.tmp$data.counting) <- colnames(pa.data.tmp$pa.tab)
    
    
      
    BFD <- BIOMOD.formated.data(sp=pa.data.tmp$sp,
                                env=pa.data.tmp$env,
                                xy=as.data.frame(pa.data.tmp$xy),
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy,
                                na.rm=na.rm) # because check is already done
    
    if(inherits(env,'Raster')){
      
      ## create data.mask for ploting
      data.mask.tmp <- reclassify(raster:::subset(env,1), c(-Inf,Inf,-1))
      data.mask <- stack(data.mask.tmp)
      xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1),]
      xy_abs <- pa.data.tmp$xy[which(pa.data.tmp$sp==0),]
      if(nrow(xy_pres)){
        data.mask[cellFromXY(data.mask.tmp, xy_pres)] <- 1
      }
      if(nrow(xy_abs)){
        data.mask[cellFromXY(data.mask.tmp, xy_abs)] <- 0
      }
      names(data.mask) <- "input_data"
      
      ## add eval data
      if(BFD@has.data.eval){
        ### TO DO
        
      }
      
      for(pa in 1:ncol(as.data.frame(pa.data.tmp$pa.tab))){
        data.mask.tmp2 <- data.mask.tmp
        
        xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1 & as.data.frame(pa.data.tmp$pa.tab)[,pa]) ,]
        xy_abs <- pa.data.tmp$xy[which( (pa.data.tmp$sp!=1 | is.na(pa.data.tmp$sp)) & as.data.frame(pa.data.tmp$pa.tab)[,pa]) ,]
        
        if(nrow(xy_pres)){
          id_pres <- cellFromXY(data.mask.tmp, xy_pres)
          data.mask.tmp2[id_pres] <- 1
        }
        
        if(nrow(xy_abs)){
          id_abs <- cellFromXY(data.mask.tmp, xy_abs)
          data.mask.tmp2[id_abs] <- 0
        }
        
        data.mask <- addLayer(data.mask, data.mask.tmp2)
      }
      
      names(data.mask) <- c("input_data", colnames(as.data.frame(pa.data.tmp$pa.tab)))
      
    } else{
      data.mask <- stack()
    }
    

      
    BFDP <- new('BIOMOD.formated.data.PA',
                sp.name = BFD@sp.name,
                coord = BFD@coord,
#                 data.counting = cbind(BFD@data.counting,pa.data.tmp$data.counting) ,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                data.mask = data.mask,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA = as.data.frame(pa.data.tmp$pa.tab),
                PA.strategy = PA.strategy)
    
    rm(list='BFD')
  } else {
    cat("\n   ! PA selection not done", fill=.Options$width)
      
    BFDP <- BIOMOD.formated.data(sp=sp,
                                env=env,
                                xy=xy,
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy)
  
  }
  
  rm(list = "pa.data.tmp" )
  
  .bmCat('Done')
  return(BFDP)

}


# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.formated.data.PA'),
          function(x,coord=NULL,col=NULL){
            if(nlayers(x@data.mask)>0){
              
              ## define the breaks of the color key
              my.at <- seq(-1.5,1.5,by=1)
              ## the labels will be placed vertically centered
              my.labs.at <- seq(-1,1,by=1)
              ## define the labels
              my.lab = c("undifined","absences","presences")
              
              levelplot(x@data.mask, at=my.at, margin=T, col.regions=c("lightgrey","red4","green4"),
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))
              
            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'orange', 'grey')
              }
              
              # plot data
              par(mfrow=c(.CleverCut(ncol(x@PA)+1)))
              # all points (~mask)
              plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name," original data", sep=""), pch=20 )
              # presences 
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)
              
              # PA data
              for(i in 1:ncol(x@PA)){
                # all points (~mask)
                plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                     main = paste(x@sp.name," Pseudo Absences ", i, sep=""), pch=20 )
                # presences 
                points(x=x@coord[(x@data.species == 1) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 1) & x@PA[,i],2],
                       col=col[1],pch=18)
                # true absences
                points(x=x@coord[(x@data.species == 0) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 0) & x@PA[,i],2],
                       col=col[2],pch=18)
                # PA
                points(x=x@coord[is.na(x@data.species) & x@PA[,i],1],
                       y=x@coord[is.na(x@data.species) & x@PA[,i],2],
                       col=col[3],pch=18)
              }
            } })

setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object){
            .bmCat("'BIOMOD.formated.data.PA'")
            cat("\nsp.name = ", object@sp.name,fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            cat("\n\n", ncol(object@PA), 'Pseudo Absences dataset available (', colnames(object@PA),") with ",
                sum(object@PA[,1], na.rm=T) - sum(object@data.species, na.rm=TRUE), 'absences in each (true abs + pseudo abs)', fill=.Options$width)
            .bmCat()
          })


setClass("BIOMOD.Model.Options", 
         representation(GLM = "list", 
                        GBM = "list",
                        GAM = "list",
                        CTA = "list",
                        ANN = "list",
                        SRE = "list",
                        FDA = "list",
                        MARS = "list",
                        RF = "list",
                        MAXENT = "list"),
         
         prototype(GLM = list( type = 'quadratic',
                               interaction.level = 0,
                               myFormula = NULL,
                               test = 'AIC',
                               family = binomial(link='logit'),
                               mustart = 0.5,
                               control = glm.control(maxit = 50)),
                   
                   GBM = list(  distribution = 'bernoulli',
                                n.trees = 1000,
                                interaction.depth = 7,
                                n.minobsinnode = 5,
                                shrinkage = 0.001,
                                bag.fraction = 0.5,
                                train.fraction = 1,
                                cv.folds = 3,
                                keep.data = FALSE,
                                verbose = FALSE,
#                                 class.stratify.cv = 'bernoulli',
                                perf.method = 'cv'),
                   
                   GAM = list( algo = "GAM_mgcv",
                               type = "s_smoother",
                               k = NULL,
                               interaction.level = 0,
                               myFormula = NULL,
                               family = binomial(link='logit'),
                               control = list(epsilon = 1e-06, trace = FALSE ,maxit = 100),
                               method = "GCV.Cp",
                               optimizer=c("outer","newton"),
                               select=FALSE,
                               knots=NULL,
                               paraPen=NULL), 
                   
                   CTA = list(method = 'class',
                              parms = 'default',
#                               control = rpart.control(xval = 5, minbucket = 5, minsplit = 5,
#                                                       cp = 0.001, maxdepth = 25),
                              control = list(xval = 5, minbucket = 5, minsplit = 5,
                                                      cp = 0.001, maxdepth = 25),                              
                              cost = NULL ),
                                                      
                   ANN = list(NbCV = 5,
                              rang = 0.1,
                              maxit = 200),
                   
                   SRE = list(quant = 0.025),
                   
                   FDA = list(method = 'mars'),
                   
                   MARS = list(degree = 2,
                               penalty = 2,
                               thresh = 0.001,
                               prune = TRUE),
                   
                   RF = list(do.classif = TRUE,
                             ntree = 500,
                             mtry = 'default',
                             nodesize = 5,
                             maxnodes= NULL),
                   
                   MAXENT = list(path_to_maxent.jar = getwd(),
                                 maximumiterations = 200,
                                 visible = FALSE,
                                 linear = TRUE,
                                 quadratic = TRUE,
                                 product = TRUE,
                                 threshold = TRUE,
                                 hinge = TRUE,
                                 lq2lqptthreshold = 80,
                                 l2lqthreshold = 10,
                                 hingethreshold = 15,
                                 beta_threshold = -1.0,
                                 beta_categorical = -1.0,
                                 beta_lqp = -1.0,
                                 beta_hinge = -1.0,
                                 defaultprevalence = 0.5)
                   
                   ),
         validity = function(object){
           test <- TRUE
           ## GLM ##
           if(!(object@GLM$type %in% c('simple','quadratic','polynomial','user.defined'))){ cat("\nGLM$type must be 'simple',  'quadratic', 'polynomial' or 'user.defined'"); test <- FALSE}
           if(!is.numeric(object@GLM$interaction.level)){ cat("\nGLM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GLM$interaction.level < 0 | object@GLM$interaction.level%%1!=0){ cat("\nGLM$interaction.level must be a positive integer"); test <- FALSE }
           }
           if(!is.null(object@GLM$myFormula)) if(class(object@GLM$myFormula) != "formula"){ cat("\nGLM$myFormula must be NULL or a formula object"); test <- FALSE }
           if(!(object@GLM$test %in% c('AIC','BIC','none'))){ cat("\nGLM$test must be 'AIC','BIC' or 'none'"); test <- FALSE}
           fam <- 'none'
           if(class(object@GLM$family) != "family"){ cat("\nGLM$family must be a valid family object"); test <- FALSE }
           if(!is.list(object@GLM$control)){cat("\nGLM$control must be a list like that returned by glm.control"); test <- FALSE}
           
           
           ## GBM ##
           if(! object@GBM$distribution %in% c("bernoulli","huberized", "multinomial", "adaboost")){cat("\nGBM$distribution must be 'bernoulli', 'huberized', 'multinomial'or  'adaboost'") }

           if(!is.numeric(object@GBM$n.trees)){ cat("\nGBM$n.trees must be a integer"); test <- FALSE } else{
             if(object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees){ cat("\nGBM$n.trees must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$interaction.depth)){ cat("\nGBM$interaction.depth must be a integer"); test <- FALSE } else{
             if(object@GBM$interaction.depth < 0 | object@GBM$interaction.depth%%1!=0){ cat("\nGBM$interaction.depth must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$n.minobsinnode)){ cat("\nGBM$n.minobsinnode must be a integer"); test <- FALSE } else{
             if(object@GBM$n.minobsinnode < 0 | object@GBM$n.minobsinnode%%1!=0){ cat("\nGBM$n.minobsinnode must positive "); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$shrinkage)){ cat("\nGBM$shrinkage must be a numeric"); test <- FALSE } else{
             if(object@GBM$shrinkage < 0 ){ cat("\nGBM$shrinkage must positive "); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$bag.fraction)){ cat("\nGBM$bag.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$bag.fraction < 0 | object@GBM$bag.fraction > 1){ cat("\nGBM$bag.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$train.fraction)){ cat("\nGBM$train.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$train.fraction < 0 | object@GBM$train.fraction > 1){ cat("\nGBM$train.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$cv.folds)){ cat("\nGBM$cv.folds must be a integer"); test <- FALSE } else{
             if(object@GBM$cv.folds < 0 | object@GBM$cv.folds%%1!=0){ cat("\nGBM$cv.folds must be a positive integer"); test <- FALSE }
           }
            
           if(!is.logical(object@GBM$keep.data)){ cat("\nGBM$keep.data must be a logical"); test <- FALSE }
           
           if(!is.logical(object@GBM$verbose)){ cat("\nGBM$verbose must be a logical"); test <- FALSE }

#            if(! object@GBM$class.stratify.cv %in% c("bernoulli", "multinomial")){cat("\nGBM$class.stratify.cv must be 'bernoulli' or 'multinomial'") }
           
           if(! object@GBM$perf.method %in% c('OOB', 'test', 'cv')){cat("\nGBM$perf.method must be 'OOB', 'test' or 'cv'"); test <- FALSE }
           
           
           ## GAM ##
           if(! object@GAM$algo %in% c('GAM_mgcv','GAM_gam', 'BAM_mgcv')){cat("\nGAM$algo must be 'GAM_mgcv','GAM_gam' or  'BAM_mgcv'"); test <- FALSE }
           
           if(! object@GAM$type %in% c('s_smoother','s', 'lo', 'te')){cat("\nGAM$type must be c('s_smoother','s', 'lo' or 'te'"); test <- FALSE }

           if(! is.null(object@GAM$k)){
             if(! is.numeric(object@GAM$k)  ){ cat("\nGAM$k must be a integer"); test <- FALSE } else{
               if(object@GAM$k < -1 | object@GAM$k%%1!=0){ cat("\nGAM$k must be > -1"); test <- FALSE }
             }
           }

           
           if(!is.numeric(object@GAM$interaction.level)){ cat("\nGAM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GAM$interaction.level < 0 | object@GAM$interaction.level%%1!=0){ cat("\nGAM$interaction.level must be a positive integer"); test <- FALSE }
           }
           
           if(!is.null(object@GAM$myFormula)) if(class(object@GAM$myFormula) != "formula"){ cat("\nGAM$myFormula must be NULL or a formula object"); test <- FALSE }
           
           if(class(object@GAM$family) != "family"){ cat("\nGAM$family must be a valid family object"); test <- FALSE }
           
           if(!is.list(object@GAM$control)){cat("\nGAM$control must be a list like that returned by gam.control"); test <- FALSE}
           if(! object@GAM$method %in% c('GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML', 'P-ML')){cat("\nGAM$method must be 'GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML'or 'P-ML'"); test <- FALSE}
           
           if(sum(! object@GAM$optimizer %in% c('perf','outer', 'newton', 'bfgs', 'optim', 'nlm', 'nlm.fd')) > 0 ){cat("\nGAM$optimizer bad definition (see ?mgcv::gam)") ; test <- FALSE}
           
           if(!is.logical(object@GAM$select)){ cat("\nGAM$select must be a logical"); test <- FALSE }

#            knots=NULL,
#            paraPen=NULL
           
           
           ## CTA ##
           if(! object@CTA$method %in% c( 'anova', 'poisson', 'class', 'exp')){cat("\nCTA$method must be 'anova', 'poisson', 'class' or 'exp'"); test <- FALSE }
           
           #parms = 'default',

           if(!is.list(object@CTA$control)){cat("\nCTA$control must be a list like that returned by rpart.control"); test <- FALSE}                           
           if(length(object@CTA$cost)){
             if(!is.numeric(object@CTA$cost)){cat("\nCTA$cost must be a non negative cost vector"); test <- FALSE}
           }
           
           
           
           ## ANN ##
           if(!is.numeric(object@ANN$NbCV)){ cat("\nANN$NbCV must be a integer"); test <- FALSE } else{
             if(object@ANN$NbCV < 0 | object@ANN$NbCV%%1!=0){ cat("\nANN$NbCV must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@ANN$rang)){ cat("\nANN$rang must be a numeric"); test <- FALSE } else{
             if(object@ANN$rang < 0 ){ cat("\nANN$rang must be positive"); test <- FALSE }
           }
           
           if(!is.numeric(object@ANN$maxit)){ cat("\nANN$maxit must be a integer"); test <- FALSE } else{
             if(object@ANN$maxit < 0 | object@ANN$maxit%%1!=0){ cat("\nANN$maxit must be a positive integer"); test <- FALSE }
           }
           
           
           
           ## FDA ##
           if(! object@FDA$method %in% c( 'polyreg', 'mars', 'bruto')){cat("\nFDA$method must be 'polyreg', 'mars' or 'bruto'"); test <- FALSE }
           
           
           
           ## SRE ##
           if(!is.numeric(object@SRE$quant)){ cat("\nSRE$quant must be a numeric"); test <- FALSE } else{
             if(object@SRE$quant >= 0.5 | object@SRE$quant < 0){ cat("\nSRE$quant must between 0 and 0.5"); test <- FALSE }
           }
           
           
           
           ## MARS ##
           if(!is.numeric(object@MARS$degree)){ cat("\nMARS$degree must be a integer"); test <- FALSE } else{
             if(object@MARS$degree < 0 | object@MARS$degree%%1!=0){ cat("\nMARS$degree must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@MARS$penalty)){ cat("\nMARS$penalty must be a integer"); test <- FALSE } else{
             if(object@MARS$penalty < 0 | object@MARS$penalty%%1!=0){ cat("\nMARS$penalty must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@MARS$thresh)){ cat("\nMARS$thresh must be a numeric"); test <- FALSE } else{
             if(object@MARS$thresh < 0 ){ cat("\nMARS$thresh must be positive"); test <- FALSE }
           }
           
           if(!is.logical(object@MARS$prune)){ cat("\nMARS$prune must be a logical"); test <- FALSE }
           
           
           
           ## RF ##
           if(!is.logical(object@RF$do.classif)){ cat("\nRF$do.classif must be a logical"); test <- FALSE }
           
           if(!is.numeric(object@RF$ntree)){ cat("\nRF$ntree must be a integer"); test <- FALSE } else{
             if(object@RF$ntree < 0 | object@RF$ntree%%1!=0){ cat("\nRF$ntree must be a positive integer"); test <- FALSE }
           }
           
           if( object@RF$mtry != 'default'){
             if(!is.numeric(object@RF$mtry)){ cat("\nRF$mtry must be a integer"); test <- FALSE } else{
               if(object@RF$mtry < 0 | object@RF$mtry%%1!=0){ cat("\nRF$mtry must be a positive integer"); test <- FALSE }
             }
           }
           
           if(!is.numeric(object@RF$nodesize)){ cat("\nRF$nodesize must be a integer"); test <- FALSE } else{
             if(object@RF$nodesize < 0 | object@RF$nodesize%%1!=0){ cat("\nRF$nodesize must be a positive integer"); test <- FALSE }
           }
           
           if(length(object@RF$maxnodes)){
             if(!is.numeric(object@RF$maxnodes)){ cat("\nRF$maxnodes must be a integer"); test <- FALSE } else{
               if(object@RF$maxnodes < 0 | object@RF$maxnodes%%1!=0){ cat("\nRF$maxnodes must be a positive integer"); test <- FALSE }
             }             
           }

           
           
           ## MAXENT ##
           if(!is.character(object@MAXENT$path_to_maxent.jar)){ cat("\nMAXENT$path_to_maxent.jar must be a character"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$maximumiterations)){ cat("\nMAXENT$maximumiterations must be a integer"); test <- FALSE } else{
             if(object@MAXENT$maximumiterations < 0 | object@MAXENT$maximumiterations%%1!=0){ cat("\nMAXENT$maximumiterations must be a positive integer"); test <- FALSE }
           }
           
           if(!is.logical(object@MAXENT$visible)){ cat("\nMAXENT$visible must be a logical"); test <- FALSE }
           
           if(!is.logical(object@MAXENT$linear)){ cat("\nMAXENT$linear must be a logical"); test <- FALSE }
           
           if(!is.logical(object@MAXENT$quadratic)){ cat("\nMAXENT$quadratic must be a logical"); test <- FALSE }
           
           if(!is.logical(object@MAXENT$product)){ cat("\nMAXENT$product must be a logical"); test <- FALSE }
           
           if(!is.logical(object@MAXENT$threshold)){ cat("\nMAXENT$threshold must be a logical"); test <- FALSE }
           
           if(!is.logical(object@MAXENT$hinge)){ cat("\nMAXENT$hinge must be a logical"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$lq2lqptthreshold)){ cat("\nMAXENT$lq2lqptthreshold must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$l2lqthreshold)){ cat("\nMAXENT$l2lqthreshold must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$lq2lqptthreshold)){ cat("\nMAXENT$lq2lqptthreshold must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$hingethreshold)){ cat("\nMAXENT$hingethreshold must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$beta_threshold)){ cat("\nMAXENT$beta_threshold must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$beta_categorical)){ cat("\nMAXENT$beta_categorical must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$beta_lqp)){ cat("\nMAXENT$beta_lqp must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$beta_hinge)){ cat("\nMAXENT$beta_hinge must be a numeric"); test <- FALSE }
           
           if(!is.numeric(object@MAXENT$defaultprevalence)){ cat("\nMAXENT$defaultprevalence must be a numeric"); test <- FALSE }
                    
           return(test)
         })

setMethod('show', signature('BIOMOD.Model.Options'),
          function(object){
            .bmCat(" 'BIOMOD.Model.Options' ")
            cat("\n")

            ## GLM options
            cat("\nGLM = list( type = '", object@GLM$type, "',", sep="")
            cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep="")
            cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="")
            cat("\n            test = '", object@GLM$test, "',", sep="")
            cat("\n            family = ", object@GLM$family$family,"(link = '",object@GLM$family$link,"'),", sep="")
            cat("\n            mustart = ", object@GLM$mustart, ",", sep="")
            cat("\n            control = glm.control(", .print.control(object@GLM$control), ") ),", sep="", fill=.Options$width)
            
            ## GBM options
            cat("\n")
            cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep="")
            cat("\n            n.trees = ", object@GBM$n.trees, ",", sep="")
            cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep="")
            cat("\n            n.minobsinnode = ", object@GBM$n.minobsinnode, ",", sep="")
            cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep="")
            cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep="")
            cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep="")
            cat("\n            cv.folds = ", object@GBM$cv.folds, ",", sep="")
            cat("\n            keep.data = ", object@GBM$keep.data, ",", sep="")
            cat("\n            verbose = ", object@GBM$verbose, ",", sep="")
#             cat("\n            class.stratify.cv = '", object@GBM$class.stratify.cv, "',", sep="")
            cat("\n            perf.method = '", object@GBM$perf.method, "'),", sep="")
            
            ## GAM options
            cat("\n")
            cat("\nGAM = list( algo = '", object@GAM$algo, "',", sep="")
            cat("\n            type = '", object@GAM$type, "',", sep="")
            cat("\n            k = ", ifelse(length(object@GAM$k) < 1,'NULL',object@GAM$k), ",", sep="")
            cat("\n            interaction.level = ", object@GAM$interaction.level, ",", sep="")
            cat("\n            myFormula = ", ifelse(length(object@GAM$myFormula) < 1,'NULL',paste(object@GAM$myFormula[2],object@GAM$myFormula[1],object@GAM$myFormula[3])), ",", sep="")
            cat("\n            family = ", object@GAM$family$family,"(link = '",object@GAM$family$link,"'),", sep="")
            
            if(object@GAM$algo=='GAM_mgcv'){
              cat("\n            method = '", object@GAM$method, "',", sep="")
              cat("\n            optimizer = c('", paste(object@GAM$optimizer,collapse="','"), "'),", sep="")
              cat("\n            select = ", object@GAM$select, ",", sep="")
              cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1,'NULL',"'user.defined'"), ",", sep="")
              cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1,'NULL',"'user.defined'"), ",", sep="")
            }
            
            cat("\n            control = list(", .print.control(object@GAM$control), ") ),", sep="", fill=.Options$width)
            
            

            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep="")
            cat("\n            parms = '", object@CTA$parms, "',", sep="")
            cat("\n            cost = ", ifelse(length(object@CTA$cost)<1,'NULL',object@CTA$cost), ",", sep="")
            cat("\n            control = list(", .print.control(object@CTA$control), ") ),", sep="", fill=.Options$width)
            
            ## ANN options
            cat("\n")
            cat("\nANN = list( NbCV = ", object@ANN$NbCV, ",", sep="")
            cat("\n            rang = ", object@ANN$rang, ",", sep="")
            cat("\n            maxit = ", object@ANN$maxit, "),", sep="")
            
            ## SRE options
            cat("\n")
            cat("\nSRE = list( quant = ", object@SRE$quant, "),", sep="")
                                                      
            ## FDA options
            cat("\n")
            cat("\nFDA = list( method = '", object@FDA$method, "'),", sep="")
                   
            ## MARS options
            cat("\n")
            cat("\nMARS = list( degree = ", object@MARS$degree, ",", sep="")
            cat("\n             penalty = ", object@MARS$penalty, ",", sep="")
            cat("\n             thresh = ", object@MARS$thresh, ",", sep="")
            cat("\n             prune = ", object@MARS$prune, "),", sep="")
            
            ## RF options
            cat("\n")
            cat("\nRF = list( do.classif = ", object@RF$do.classif, ",", sep="")
            cat("\n           ntree = ", object@RF$ntree, ",", sep="")
            cat("\n           mtry = '", object@RF$mtry, "',", sep="")
            cat("\n           nodesize = ", object@RF$nodesize, ",", sep="")
            cat("\n           maxnodes = ", ifelse(length(object@RF$maxnodes) < 1,'NULL',object@RF$maxnodes), "),", sep="")
                   
            ## MAXENT options
            cat("\n")
            cat("\nMAXENT = list( path_to_maxent.jar = '", object@MAXENT$path_to_maxent.jar, "',", sep="")
            cat("\n               maximumiterations = ", object@MAXENT$maximumiterations, ",", sep="")
            cat("\n               visible = ", object@MAXENT$visible, ",", sep="")
            cat("\n               linear = ", object@MAXENT$linear, ",", sep="")
            cat("\n               quadratic = ", object@MAXENT$quadratic, ",", sep="")
            cat("\n               product = ", object@MAXENT$product, ",", sep="")
            cat("\n               threshold = ", object@MAXENT$threshold, ",", sep="")
            cat("\n               hinge = ", object@MAXENT$hinge, ",", sep="")
            cat("\n               lq2lqptthreshold = ", object@MAXENT$lq2lqptthreshold, ",", sep="")
            cat("\n               l2lqthreshold = ", object@MAXENT$l2lqthreshold, ",", sep="")
            cat("\n               hingethreshold = ", object@MAXENT$hingethreshold, ",", sep="")
            cat("\n               beta_threshold = ", object@MAXENT$beta_threshold, ",", sep="")
            cat("\n               beta_categorical = ", object@MAXENT$beta_categorical, ",", sep="")
            cat("\n               beta_lqp = ", object@MAXENT$beta_lqp, ",", sep="")
            cat("\n               beta_hinge = ", object@MAXENT$beta_hinge, ",", sep="")
            cat("\n               defaultprevalence = ", object@MAXENT$defaultprevalence, ")", sep="")

            .bmCat()
          })

.print.control <- function(ctrl){
  out <-  paste(names(ctrl)[1], " = ", ctrl[[1]], sep="")
  
  if(length(ctrl) > 1){
    i=2
    while(i <= length(ctrl)){
      if(is.list(ctrl[[i]])){
        out <- c(out, paste(", ", names(ctrl)[i], " = list(",
                            paste(names(ctrl[[i]]), "=",unlist(ctrl[[i]]), sep="", collapse=", "),")", sep=""))
#         i <- i+length(ctrl[[i]])
        i <- i+1
      } else {
        out <- c(out, paste(", ", names(ctrl)[i], " = ", ctrl[[i]], sep=""))
        i <- i+1
      }
    }    
  }
#   return(toString(out))
  return(out)
  
}
####################################################################################################
### BIOMOD Storing Results Objects #################################################################
####################################################################################################
setClass("BIOMOD.stored.data",
         representation(inMemory = 'logical',
                        link = 'character'),
         prototype(inMemory=FALSE,
                   link = ''),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.array",
         contains = "BIOMOD.stored.data",
         representation(val = 'array'),
         prototype(val = array()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.raster.stack",
         contains = "BIOMOD.stored.data",
         representation(val = 'RasterStack'),
         prototype(val = stack()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.files",
         contains = "BIOMOD.stored.data",
         representation(val = 'character'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.formated.data",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.formated.data'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.models.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.Model.Options'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.models.out",
         representation(modeling.id = 'character', 
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.computed = 'character',
                        models.failed = 'character',
                        has.evaluation.data = 'logical',
                        rescal.all.models = 'logical',
                        models.evaluation = 'BIOMOD.stored.array',
                        variables.importances = 'BIOMOD.stored.array',
                        models.prediction = 'BIOMOD.stored.array',
                        models.prediction.eval = 'BIOMOD.stored.array',
                        formated.input.data = 'BIOMOD.stored.formated.data',
                        calib.lines = 'BIOMOD.stored.array',
                        models.options = 'BIOMOD.stored.models.options',
                        link = 'character'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   sp.name='',
                   expl.var.names = '',
                   models.computed='',
                   models.failed='',
                   has.evaluation.data=FALSE,
                   rescal.all.models=TRUE, 
                   models.evaluation = new('BIOMOD.stored.array'),
                   variables.importances = new('BIOMOD.stored.array'),
                   models.prediction = new('BIOMOD.stored.array'),
                   models.prediction.eval = new('BIOMOD.stored.array'),
                   formated.input.data = new('BIOMOD.stored.formated.data'),
                   calib.lines = new('BIOMOD.stored.array'),
                   models.options = new('BIOMOD.stored.models.options'),
                   link=''),
         validity = function(object){
           return(TRUE)
           })

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setMethod('show', signature('BIOMOD.models.out'),
          function(object){
            .bmCat("BIOMOD.models.out")
            cat("\nModeling id :", object@modeling.id, fill=.Options$width)
            cat("\nSpecies modeled :", object@sp.name, fill=.Options$width)
            cat("\nConsidered variables :", object@expl.var.names, fill=.Options$width)
            
            cat("\n\nComputed Models : ", object@models.computed, fill=.Options$width)
            cat("\n\nFailed Models : ", object@models.failed, fill=.Options$width)
            .bmCat()
          })


setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

### GETTEURS ###
setGeneric("getModelsPrediction",
           function(obj,...){
             standardGeneric("getModelsPrediction")
           })

setMethod("getModelsPrediction", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            if(!as.data.frame){
              if(obj@models.prediction@inMemory ){
                return(obj@models.prediction@val)
              } else{
                if(obj@models.prediction@link != ''){
#                   load(obj@models.prediction@link)
#                   return(models.prediction)
                  
                  return(get(load(obj@models.prediction@link)))
                } else{ return(NULL) }
              }              
            } else {
              if(obj@models.prediction@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                return(mod.pred)
              } else{
                if(obj@models.prediction@link != ''){
#                   load(obj@models.prediction@link)
#                   mod.pred <- as.data.frame(models.prediction)
                  mod.pred <- as.data.frame(get(load(obj@models.prediction@link)))                  
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                   function(x){
                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                     }))
                return(mod.pred)
                } else{ return(NULL) }
              }
              
            }
          }
          )




setGeneric("getModelsPredictionEval",
           function(obj,...){
             standardGeneric("getModelsPredictionEval")
           })

setMethod("getModelsPredictionEval", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            if(!as.data.frame){
              if(obj@models.prediction.eval@inMemory ){
                return(obj@models.prediction.eval@val)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  load(obj@models.prediction.eval@link)
                  return(models.prediction.eval)
                } else{ return(NULL) }
              }              
            } else {
              if(obj@models.prediction.eval@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction.eval@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                return(mod.pred)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  load(obj@models.prediction.eval@link)
                  mod.pred <- as.data.frame(models.prediction.eval)
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                   function(x){
                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                     }))
                return(mod.pred)
                } else{ return(NULL) }
              }
              
            }
          }
          )


setGeneric("getModelsEvaluations",
           function(obj){
             standardGeneric("getModelsEvaluations")
           })

setMethod("getModelsEvaluations", "BIOMOD.models.out",
          function(obj){
            if(obj@models.evaluation@inMemory ){
              return(obj@models.evaluation@val)
            } else{
              if(obj@models.evaluation@link != ''){
#                 load(obj@models.evaluation@link)
#                 return(models.evaluation)
                return(get(load(obj@models.evaluation@link)))                
              } else{ return(NA) }
            }
          }
          )


setGeneric("getModelsVarImport",
           function(obj){
             standardGeneric("getModelsVarImport")
           })

setMethod("getModelsVarImport", "BIOMOD.models.out",
          function(obj){
            if(obj@variables.importances@inMemory ){
              return(obj@variables.importances@val)
            } else{
              if(obj@variables.importances@link != ''){
#                 load(obj@variables.importances@link)
#                 return(variables.importances)
                return(get(load(obj@variables.importances@link)))
              } else{ return(NA) }
            }
          }
          )


setGeneric("getModelsOptions",
           function(obj){
             standardGeneric("getModelsOptions")
           })

setMethod("getModelsOptions", "BIOMOD.models.out",
          function(obj){
            if(obj@models.options@inMemory ){
              return(obj@models.options@val)
            } else{
              if(obj@models.options@link != ''){
#                 load(obj@models.options@link)
#                 return(models.options)
                return(get(load(obj@models.options@link)))                
              } else{ return(NA) }
            }
          }
          )

setGeneric("getModelsInputData",
           function(obj, ...){
             standardGeneric("getModelsInputData")
           })

setMethod("getModelsInputData", "BIOMOD.models.out",
          function(obj, subinfo = NULL){
            if(is.null(subinfo)){
              if(obj@formated.input.data@inMemory ){
                return(obj@formated.input.data@val)
              } else{
                if(obj@formated.input.data@link != ''){
                  load(obj@formated.input.data@link)
                  return(data)
                } else{ return(NA) }
              }              
            } else if(subinfo == 'MinMax'){
              return(apply(getModelsInputData(obj)@data.env.var,2, function(x){
                  if(is.numeric(x)){
                    return( list(min = min(x,na.rm=T), max = max(x, na.rm=T) ) )
                  } else if(is.factor(x)){
                    return(list(levels = levels(x)))
                  }
                }) )
            } else if(subinfo == 'expl.var'){
              return(getModelsInputData(obj)@data.env.var)
            } else if(subinfo == 'expl.var.names'){
              return(obj@expl.var.names)
            } else if(subinfo == 'resp.var'){
              return(getModelsInputData(obj)@data.species)
            } else{
              stop("Unknow subinfo tag")
            }

          }
          )

setGeneric("getModelsBuiltModels",
           function(obj){
             standardGeneric("getModelsBuiltModels")
           })

setMethod("getModelsBuiltModels", "BIOMOD.models.out",
          function(obj){
            return(obj@models.computed)
          }
          )


setGeneric("RemoveProperly",
           function(obj,obj.name=deparse(substitute(obj))){
             standardGeneric("RemoveProperly")
           })

setMethod("RemoveProperly", "BIOMOD.models.out",
          function(obj, obj.name=deparse(substitute(obj))){
            cat("\n\t> Removing .BIOMOD_DATA files...")
            unlink(file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id), recursive=T, force=TRUE)
            cat("\n\t> Removing models...")
            unlink(file.path(obj@sp.name, "models", obj@modeling.id), recursive=T, force=TRUE)
            cat("\n\t> Removing object hard drive copy...")
            unlink(obj@link, recursive=T, force=TRUE)
            cat("\n\t> Removing object from memory...")
            rm(list=obj.name,envir=sys.frame(-2))
            cat("\nCompleted!")
          })


####################################################################################################
### BIOMOD Storing Projection Objects ##############################################################
####################################################################################################
# setClass("BIOMOD.projection",
#          representation(proj.names = 'character',
#                         sp.name = 'character',
#                         expl.var.names = 'character',
#                         models.computed = 'character',
#                         models.failed = 'character',
#                         models.thresholds = 'BIOMOD.stored.array',
#                         models.prediction = 'BIOMOD.stored.array',
#                         formated.input.data = 'BIOMOD.stored.formated.data',
#                         calib.lines = 'BIOMOD.stored.array',
#                         models.options = 'BIOMOD.stored.models.options'),
#          prototype(sp.name='',
#                    expl.var.names = '',
#                    models.computed='',
#                    models.failed='',
#                    models.evaluation = new('BIOMOD.stored.array'),
#                    variables.importances = new('BIOMOD.stored.array'),
#                    models.prediction = new('BIOMOD.stored.array'),
#                    formated.input.data = new('BIOMOD.stored.formated.data'),
#                    calib.lines = new('BIOMOD.stored.array'),
#                    models.options = new('BIOMOD.stored.models.options')),
#          validity = function(object){
#            return(TRUE)
#            })

setClass("BIOMOD.projection.out",
         representation(proj.names = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.projected = 'character',
                        scaled.models = 'logical',
                        modeling.object = 'BIOMOD.stored.data',
                        modeling.object.id = 'character',
                        type = 'character',
                        proj = 'BIOMOD.stored.data',
                        xy.coord = 'matrix'),
         prototype(proj.names = '',
                   sp.name='',
                   expl.var.names='',
                   models.projected='',
                   scaled.models=TRUE,
                   modeling.object.id='',
                   type='',
                   xy.coord = matrix()),
         validity = function(object){
           return(TRUE)
           })

setGeneric("getProjection",
           function(obj, ...){
             standardGeneric("getProjection")
           })

setMethod("getProjection", "BIOMOD.projection.out",
          function(obj, model = NULL, as.data.frame = FALSE){
            if(!as.data.frame & is.null(model)){
              if(obj@proj@inMemory ){
                return(obj@proj@val)
              } else {
                if( grepl(".RData", obj@proj@link) ){
                  return(get(load(obj@proj@link)))
                } else if(grepl(".grd", obj@proj@link) | grepl(".img", obj@proj@link)){
                  return(stack(obj@proj@link))
                } else {
                  cat("\n!! You have to load your projections by yourself !!")
                  return(NA) 
                }
              } 
            } else if(as.data.frame){
              if(obj@proj@inMemory ){
                proj <- as.data.frame(obj@proj@val)
                names(proj) <- unlist(lapply(strsplit(names(proj),".", fixed=TRUE), 
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                return(proj)
              } else{
                if(obj@proj@link != ''){
                  load(obj@proj@link)
                  project <- as.data.frame(proj)
                  names(project) <- unlist(lapply(strsplit(names(project),".", fixed=TRUE), 
                                   function(x){
                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                     }))
                return(project)
                } else{ return(NA) }
              }
            }
              
          }
          )

# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.projection.out'),
          function(x,col=NULL, str.grep=NULL){
            if(class(x@proj) == "BIOMOD.stored.raster.stack"){
              ## define the breaks of the color key
              my.at <- seq(0,1000,by=100)
              ## the labels will be placed vertically centered
              my.labs.at <- seq(0,1000,by=250)
              ## define the labels
              my.lab <- seq(0,1000,by=250)
              ## define colors
              my.col <- colorRampPalette(c("red4","orange4","yellow4","green4"))(100)
              
              
              if(is.null(str.grep)){
                levelplot(getProjection(x),
                          at=my.at, margin=T, col.regions=my.col,
                          main=paste(x@sp.name,x@proj.names,"projections"),
                          colorkey=list(labels=list(
                            labels=my.lab,
                            at=my.labs.at)))
                
              } else if(length(grep(str.grep, x@models.projected,value=T))>0){
                levelplot(raster:::subset(getProjection(x), grep(str.grep, x@models.projected,value=T)),
                          at=my.at, margin=T, col.regions=my.col,
                          main=paste(x@sp.name,x@proj.names,"projections"),
                          colorkey=list(labels=list(
                            labels=my.lab,
                            at=my.labs.at)))
                
              } else{ stop("invalid str.grep arg")}
              
            } else if(class(x@proj) == "BIOMOD.stored.array"){
              if(ncol(x@xy.coord) != 2){
                cat("\n ! Impossible to plot projections because xy coordinates are not available !")
              } else {
                if(is.null(str.grep)){
                  multiple.plot(Data = getProjection(x,as.data.frame=T), coor = x@xy.coord)
                } else if(length(grep(str.grep, x@models.projected,value=T))>0){
                  multiple.plot(Data = getProjection(x, model=grep(str.grep, x@models.projected,value=T) ,as.data.frame=T), coor = x@xy.coord)
                } else{ stop("invalid str.grep arg")}               
              }

              
            } else {cat("\n !  Biomod Projection plotting issue !", fill=.Options$width)}

          })

setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            .bmCat("'BIOMOD.projection.out'")
            cat("\nProjection directory :", paste(object@sp.name,"/",object@proj.names, sep=""), fill=.Options$width)
            cat("\n")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodeling id :", object@modeling.object.id ,"(",object@modeling.object@link ,")", fill=.Options$width)
            cat("\nmodels projected :", toString(object@models.projected), fill=.Options$width)

            .bmCat()
          })

setGeneric("free",
           function(obj){
             standardGeneric("free")
           })

setMethod(f='free', 
                 signature='BIOMOD.projection.out',
                 definition = function(obj){
                   if(inherits(obj@proj,"BIOMOD.stored.array")){
                     obj@proj@val = array()
                   } else if(inherits(obj@proj,"BIOMOD.stored.raster.stack")){
                     obj@proj@val = stack()
                   } else{
                     obj@proj@val = NULL
                   }
                   obj@proj@inMemory = FALSE
                   
                   return(obj)
                 })

####################################################################################################
### BIOMOD Storing Ensemble Modeling Objects #######################################################
####################################################################################################           
setClass("BIOMOD.EnsembleModeling.out",
         representation(sp.name = 'character',
                        expl.var.names = 'character',
                        models.out.obj = 'BIOMOD.stored.models.out',
                        eval.metric = 'character',
                        eval.metric.quality.threshold = 'numeric',
                        em.computed = 'character',
                        em.by = 'character',
#                         em.models.kept = 'list',
#                         em.prediction = 'BIOMOD.stored.array',
#                         em.evaluation = 'BIOMOD.stored.array',
                        em.res = 'list',
                        em.ci.alpha = 'numeric',
                        em.weight = 'list',
                        em.bin.tresh = 'list'),
         prototype( sp.name = '',
                    expl.var.names = '',
                    models.out.obj = new('BIOMOD.stored.models.out'),
                    eval.metric = '',
                    eval.metric.quality.threshold = NULL,
#                     em.computed = '',
#                     em.models.kept = NULL,
#                     em.prediction = NULL,
#                     em.evaluation = NULL,
                    em.res = list(),
                    em.ci.alpha = 0.05,
                    em.weight = list(),
                    em.bin.tresh = list()),
         validity = function(object){
           return(TRUE)
           })


setMethod('show', signature('BIOMOD.EnsembleModeling.out'),
          function(object){
            .bmCat("'BIOMOD.EnsembleModeling.out'")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill=.Options$width)

            .bmCat()
          })


setGeneric("getEMalgos",
           function(obj, model){
             standardGeneric("getEMalgos")
           })

setMethod("getEMalgos", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            if(is.character(model) | is.numeric(model)){
              return(obj@em.res[[model]]$em.algo)
            } else{
              return(NULL)
            }

          }
          )

setGeneric("getEMkeptModels",
           function(obj, model){
             standardGeneric("getEMkeptModels")
           })

setMethod("getEMkeptModels", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            if(is.character(model) | is.numeric(model)){
              return(obj@em.res[[model]]$em.models.kept)
            } else{
              return(NULL)
            }

          }
          )

setGeneric("getEMeval",
           function(obj, ...){
             standardGeneric("getEMeval")
           })

setMethod("getEMeval", "BIOMOD.EnsembleModeling.out",
          function(obj, model=NULL, met=NULL){
            if(is.null(model)){
              model <- obj@em.computed
            }
            if(is.character(model) | is.numeric(model)){
              lout <- list()
              for(mod in model){
                if(is.null(met)){
                  lout[[mod]] <- obj@em.res[[mod]]$em.cross.validation[,,,drop=F]
                } else if(!is.null(meth)){
                  lout[[mod]] <- (obj@em.res[[mod]]$em.cross.validation[met,,,drop=F])
                } 
              }
              return(lout)
            } else{
              return(NULL)
            }

          }
          )

setGeneric("getEMbuiltModels",
           function(obj){
             standardGeneric("getEMbuiltModels")
           })

setMethod("getEMbuiltModels", "BIOMOD.EnsembleModeling.out",
          function(obj){
            return(obj@em.computed)
          })

####################################################################################################
### BIOMOD Storing Ensemble Forecasting Objects ####################################################
####################################################################################################

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# if( !isGeneric( ".Models.prepare.data" ) ) {
  setGeneric( ".Models.prepare.data", 
              def = function(data, ...){
                      standardGeneric( ".Models.prepare.data" )
                      } )
# }

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data'),
          function(data, NbRunEval, DataSplit, Yweights=NULL, Prevalence=NULL, do.full.models=TRUE, DataSplitTable=NULL){
            list.out <- list()
            name <- paste(data@sp.name,'_AllData',sep="")
            xy <- data@coord
            dataBM <- data.frame(cbind(data@data.species,data@data.env.var))
            colnames(dataBM)[1] <- data@sp.name
            
            # dealing with evaluation data
            if(data@has.data.eval){
              evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var[,,drop=FALSE]))
              colnames(evalDataBM)[1] <- data@sp.name
              eval.xy <- data@eval.coord
            } else{ evalDataBM <- eval.xy <- NULL }
            
            ### Calib/Valid lines
            if(!is.null(DataSplitTable)){
              cat("\n*** DataSplitTable is not NULL")
              calibLines <- DataSplitTable
              colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
            } else {
              if(NbRunEval == 0){ # take all available data
                calibLines <- matrix(rep(TRUE,length(data@data.species)),ncol=1)
                colnames(calibLines) <- '_Full'
              } else {
                calibLines <- .SampleMat(data@data.species, DataSplit, NbRunEval)                    
                if(do.full.models){
                  calibLines <- cbind(calibLines, rep(TRUE,length(data@data.species)))
                  colnames(calibLines)[NbRunEval+1] <- '_Full'
                }
              }
            }
              
            if(is.null(Yweights)){ # 1 for all points
              if(!is.null(Prevalence)){
                cat("\n\t> Automatic weights creation to rise a", Prevalence,"prevalence")
                Yweights <- .automatic_weights_creation(data@data.species ,prev=Prevalence)
              } else{
                cat("\n\t> No weights : all observations will have the same weight")
                Yweights <- rep(1,length(data@data.species))
              }
              
            }
            list.out[[name]] <- list(name=name,
                                   xy=xy,
                                   dataBM=dataBM,
                                   calibLines=calibLines,
                                   Yweights = Yweights,
                                   evalDataBM = evalDataBM,
                                   eval.xy = eval.xy)
            return(list.out)
          })

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data.PA'),
          function(data, NbRunEval, DataSplit, Yweights=NULL, Prevalence=NULL, do.full.models=TRUE, DataSplitTable=NULL){
            list.out <- list()
            formal_weights <- Yweights
            for(pa in 1:ncol(data@PA)){
              Yweights <- formal_weights
              name <- paste(data@sp.name,"_",colnames(data@PA)[pa],sep="")
              xy <- data@coord[data@PA[,pa],]
              resp <- data@data.species[data@PA[,pa]] # response variable (with pseudo absences selected)
              resp[is.na(resp)] <- 0
              dataBM <- data.frame(cbind(resp,
                                         data@data.env.var[data@PA[,pa],,drop=FALSE]))
              colnames(dataBM)[1] <- data@sp.name
              
              ### Calib/Valid lines
              if(!is.null(DataSplitTable)){
                cat("\n*** DataSplitTable is not NULL")
                if(length(dim(DataSplitTable))==2){
                  calibLines <- DataSplitTable
                } else {
                  calibLines <- asub(DataSplitTable,pa,3,drop=TRUE)
                }
                colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
                calibLines[which(!data@PA[,pa]),] <- NA
              } else{
                if(NbRunEval == 0){ # take all available data
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=1)
                  calibLines[data@PA[,pa],1] <- TRUE
                  colnames(calibLines) <- '_Full'
                } else {
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=NbRunEval)
                  sampled.mat <- .SampleMat(data@data.species[data@PA[,pa]], DataSplit, NbRunEval)
                  calibLines[data@PA[,pa],] <- sampled.mat
                  colnames(calibLines) <- colnames(sampled.mat)
                  if(do.full.models){
                    calibLines <- cbind(calibLines, rep(NA,length(data@data.species)))
                    calibLines[data@PA[,pa],NbRunEval+1] <- TRUE
                    colnames(calibLines)[NbRunEval+1] <- '_Full'
                  }                
                }
              }
              

              
              # dealing with evaluation data
              if(data@has.data.eval){
                evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var))
                colnames(evalDataBM)[1] <- data@sp.name
                eval.xy <- data@eval.coord
              } else{ evalDataBM <- eval.xy <- NULL }

              if(is.null(Yweights)){ # prevalence of 0.5... may be parametrize
                if(is.null(Prevalence)) Prevalence <- 0.5
                
                cat("\n\t\t\t! Weights where automaticly defined for", name, "to rise a", Prevalence, "prevalence !")
                
                
                Yweights <- rep(NA, length(data@data.species))
                Yweights[data@PA[,pa]] <- .automatic_weights_creation(as.numeric(dataBM[,1]) ,prev=Prevalence)#, subset=data@PA[,pa])
              } else{
                # remove useless weights
                Yweights[!data@PA[,pa]] <- NA
              }
              
              list.out[[name]] <- list(name=name,
                                     xy=xy,
                                     dataBM=dataBM,
                                     calibLines=calibLines,
                                     Yweights = Yweights,
                                     evalDataBM = evalDataBM,
                                     eval.xy = eval.xy)
            }
            return(list.out)
          })


.automatic_weights_creation <- function(resp,prev=0.5, subset=NULL){
  if(is.null(subset)) subset<- rep(TRUE, length(resp))
  
  nbPres <- sum(resp[subset], na.rm=TRUE)
  nbAbsKept <- sum(subset, na.rm=T) - sum(resp[subset], na.rm=TRUE) # The number of true absences + pseudo absences to maintain true value of prevalence
  Yweights <- rep(1,length(resp))
  
  if(nbAbsKept > nbPres){ # code absences as 1
    Yweights[which(resp>0)] <- (prev * nbAbsKept) / (nbPres * (1-prev))
  } else{ # code presences as 1
    Yweights[which(resp==0 | is.na(resp))] <- (nbPres * (1-prev)) / (prev * nbAbsKept)
  }
  Yweights = round(Yweights[])
  Yweights[!subset] <- 0
  
  return(Yweights)
}