# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD objects definition
# Damien Georges
# 09/02/2012
# v2.0
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This file defines the BIOMOD objects and all their methods 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# We choose here to create monospecific objects to make all procedures and parallelising easier
require(sp, quietly=TRUE)
require(raster, quietly=TRUE)
require(gam, quietly=TRUE)
require(rpart, quietly=TRUE)
require(mda, quietly=TRUE)

# 1. The BIOMOD.formated.data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this object is the basic one

# 1.1 Class Definition
setClass("BIOMOD.formated.data",
         representation(sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        data.env.var = "data.frame",
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
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL ){
    
    if(is.null(eval.sp)){
      BFD <- new('BIOMOD.formated.data', 
                 coord=xy, 
                 data.species=sp, 
                 data.env.var=env, 
                 sp.name=sp.name,
                 has.data.eval=FALSE)
    } else{
      BFDeval <- BIOMOD.formated.data(sp=eval.sp,
                                      env=eval.env,
                                      xy=eval.xy,
                                      sp.name=sp.name)
      
      BFD <- new('BIOMOD.formated.data', 
                 coord=xy, 
                 data.species=sp, 
                 data.env.var=env, 
                 sp.name=sp.name,
                 has.data.eval=TRUE,
                 eval.coord = BFDeval@coord,
                 eval.data.species = BFDeval@data.species,
                 eval.data.env.var = BFDeval@data.env.var )
                 
      rm('BFDeval')
    }
    return(BFD)
	}
)

# setMethod('BIOMOD.formated.data', signature(sp='data.frame'), 
#   function(sp,env,xy=NULL,sp.name=NULL){
#     if(ncol(sp) > 1 ){
#       stop("Invalid response variable")
#     }
#     sp <- as.numeric(unlist(sp))
#     BFD <- BIOMOD.formated.data(sp,env,xy,sp.name)
#     return(BFD)
#   }
# )
# 
# setMethod('BIOMOD.formated.data', signature(sp='numeric', env='matrix' ), 
#   function(sp,env,xy=NULL,sp.name=NULL){
#     env <- as.data.frame(env)
#     BFD <- BIOMOD.formated.data(sp,env,xy,sp.name)
#     return(BFD)
#   }
# )
          
setMethod('BIOMOD.formated.data', signature(sp='numeric', env='RasterStack' ), 
  function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL){
    # take the same eval environemental variables than calibrating ones 
    if(!is.null(eval.sp)){
      if(is.null(eval.env)){
        eval.env <- as.data.frame(extract(env,eval.xy))
      }
    }
    
    if(!is.null(xy)){
      env <- as.data.frame(extract(env,xy))
    } else{
      xy <- as.data.frame(coordinates(env))
      env <- as.data.frame(extract(env,xy))
    }
    BFD <- BIOMOD.formated.data(sp,env,xy,sp.name,eval.sp, eval.env, eval.xy)
    return(BFD)
  }
)

# setMethod('BIOMOD.formated.data', signature(sp='numeric', env='SpatialPointsDataFrame' ), 
#   function(sp,env,xy=NULL,sp.name=NULL){
#     if(!is.null(xy)){
#       env <- merge(x=data.frame(env),y=as.data.frame(xy), by.x = c(dim(env)[2]+1, dim(env)[2]+2))
#       xy <- env[,((ncol(env)-1):(ncol(env)))]
#       env <- env[,1:(ncol(env)-2)]
#     } else{
#       xy <- data.frame(env)[,((dim(env)[2]+1):(dim(env)[2]+2))]
#       env <- data.frame(env)[,1:dim(env)[2]]
#     }
#     BFD <- BIOMOD.formated.data(sp,env,xy,sp.name)
#     return(BFD)
#   }
# )

# setMethod('BIOMOD.formated.data', signature(sp='RasterLayer'), 
#   function(sp,env,xy=NULL,sp.name=NULL){
#     # NOTE : if species occurences are given into 'RasterLayer', the NA's ddont correspund to no info
#     # so we have to convert it
#     if(is.null(xy)) xy <- as.data.frame(coordinates(sp)[which(!is.na(sp[])),])
#     sp <- as.numeric(na.omit(sp[]))
#     # convert sp to have the good code for pres, abs and ndata
#     sp[which(sp==-1)] <- NA
#     BFD <- BIOMOD.formated.data(sp,env,xy,sp.name)
#     return(BFD)
#   }
# )
#     
# setMethod('BIOMOD.formated.data', signature(sp='SpatialPointsDataFrame'), 
#   function(sp,env,xy=NULL,sp.name=NULL){
#     if(dim(sp)[2]>1){
#       stop("BIOMOD.formated.data is monospecific")
#     }
# #     if(is.null(xy)) xy <- coordinates(sp) #as.data.frame(coordinates(sp))
# #     sp <- merge(x=data.frame(sp),y=as.data.frame(xy), by.x = c(dim(env)[2]+1, dim(env)[2]+2))
# #     xy <- sp[,((ncol(sp)-1):ncol(sp))]
# #     sp <- sp[,1]
#     BFD <- BIOMOD.formated.data(unlist(sp@data),env,coordinates(sp),sp.name)
#     return(BFD)
#   }
# )
          



# 1.3 Other Functions
if( !isGeneric( "plot" ) ) {
  setGeneric( "plot", 
              def = function(x, ...){
  	                  standardGeneric( "plot" )
                      } )
}

setMethod('plot', signature(x='BIOMOD.formated.data'),
          function(x,coord=NULL,col=NULL){
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
        
          })

setMethod('show', signature('BIOMOD.formated.data'),
          function(object){
            cat("\n-=-=-=- 'BIOMOD.formated.data' -=-=-=-")
            cat("\nsp.name = ", object@sp.name)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset')
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n')
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :")
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset')
              cat("\n\n")
              print(summary(object@eval.data.env.var))
            }
            
            cat("\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
          })

# 2. The BIOMOD.formated.data.PA =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this class inherits from BIOMOD.formated.data and have one more slot 'PA' giving PA selected

# 2.1 Class Definition
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA = 'data.frame'),
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
                                     PA.sre.quant = 0.025){
  
  # take the same eval environemental variables than calibrating ones 
  if(!is.null(eval.sp)){
    if(is.null(eval.env)){
      if(inherits(env,'Raster')){
        eval.env <- as.data.frame(extract(env,eval.xy))
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

  pa.data.tmp <- pseudo.absences.sampling(sp = sp,
                                          env = env,
                                          nb.repet = PA.NbRep,
                                          strategy = PA.strategy,
                                          nb.points = PA.nb.absences,
                                          distMin = PA.dist.min, 
                                          distMax = PA.dist.max,
                                          quant.SRE = PA.sre.quant )
  if(!is.null(pa.data.tmp)){
      
    BFD <- BIOMOD.formated.data(sp=pa.data.tmp$sp,
                                env=pa.data.tmp$env,
                                xy=as.data.frame(pa.data.tmp$xy),
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy)
    
      
    BFDP <- new('BIOMOD.formated.data.PA',
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA = as.data.frame(pa.data.tmp$pa.tab) )
    
    rm(list='BFD')
  } else {
    cat("\n   ! PA selection not done")
      
    BFDP <- BIOMOD.formated.data(sp=pa.data.tmp$sp,
                                env=pa.data.tmp$env,
                                xy=as.data.frame(pa.data.tmp$xy),
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy)
  
  }
  
  rm(list = "pa.data.tmp" )
  
  return(BFDP)

}


# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.formated.data.PA'),
          function(x,coord=NULL,col=NULL){
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
              # PA
              points(x=x@coord[x@PA[,i],1],
                     y=x@coord[x@PA[,i],2],
                     col=col[3],pch=18)
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

setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object){
            cat("\n-=-=-=- 'BIOMOD.formated.data.PA' -=-=-=-")
            cat("\nsp.name = ", object@sp.name)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset')
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n')
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :")
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset')
              cat("\n\n")
              print(summary(object@eval.data.env.var))
            }
            
            cat("\n\n", ncol(object@PA), 'Pseudo Absences dataset available (', colnames(object@PA),") with ",
                sum(object@PA[,1], na.rm=T) - sum(object@data.species, na.rm=TRUE), 'absences in each (true abs + pseudo abs)')
            cat("\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
          })

# # 3. The BIOMOD.formated.data.indep =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# # this class inherits from BIOMOD.formated.data and have same additional plots concerning 
# # independent data        
# 
# # 3.1 Class Definition
# setClass("BIOMOD.formated.data.indep",
#          contains = "BIOMOD.formated.data",
#          representation(indep.coord = 'data.frame',
#                         indep.data.species = 'numeric',
#                         indep.data.env.var = 'data.frame'),
#          validity = function(object){
#            if(length(colnames(object@data.env.var) %in% colnames(object@indep.data.env.var))!=length(colnames(object@data.env.var))){
#              stop("Independent environemental variables must be the same than calibration ones")
#            }
#            return(TRUE)
#            })
# 
# BIOMOD.formated.data.indep = function(sp, env, xy, sp.name,
#                                       indep.sp, indep.env, indep.xy){
#   BFD1 <- BIOMOD.formated.data(sp,env,xy, sp.name)
#   BFD2 <- BIOMOD.formated.data(indep.sp,indep.env,indep.xy, sp.name)
#   
#   BFDI <- new('BIOMOD.formated.data.indep',
#               sp.name = BFD1@sp.name,
#               coord = BFD1@coord,
#               data.env.var = BFD1@data.env.var,
#               data.species = BFD1@data.species,
#               indep.coord = BFD2@coord,
#               indep.data.env.var = BFD2@data.env.var,
#               indep.data.species = BFD2@data.species)
#   rm(BFD1,BFD2)
#   return(BFDI)
# }
# 
# # 3.3 Others functions
# setMethod('plot', signature(x='BIOMOD.formated.data.indep'),
#           function(x,coord=NULL,col=NULL){
#             # coordinates checking
#             if(is.null(coord)){
#               if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
#                 stop("coordinates are required to plot your data")
#               } else {
#                 coord <- x@coord
#               }
#             }
#             
#             # colors checking
#             if(is.null(col) | length(col) < 3){
#               col = c('green', 'red', 'grey')
#             }
#             
#             # plot data
#             par(mfrow=c(.CleverCut(2)))
#             # all points (~mask)
#             plot(x=x@coord[,1], y=x@coord[,2], col=col[3], xlab = 'X', ylab = 'Y',
#                  main = paste(x@sp.name," Calibration data", sep=""), pch=20 )
#             # presences 
#             points(x=x@coord[which(x@data.species == 1),1],
#                    y=x@coord[which(x@data.species == 1),2],
#                    col=col[1],pch=18)
#             # true absences
#             points(x=x@coord[which(x@data.species == 0),1],
#                    y=x@coord[which(x@data.species == 0),2],
#                    col=col[2],pch=18)
#             
#             # Indep data
#             # all points (~mask)
#             plot(x=x@indep.coord[,1], y=x@indep.coord[,2], col=col[3], xlab = 'X', ylab = 'Y',
#                  main = paste(x@sp.name," Validation data", sep=""), pch=20 )
#             # presences 
#             points(x=x@indep.coord[which(x@indep.data.species == 1),1],
#                    y=x@indep.coord[which(x@indep.data.species == 1),2],
#                    col=col[1],pch=18)
#             # true absences
#             points(x=x@indep.coord[which(x@indep.data.species == 0),1],
#                    y=x@indep.coord[which(x@indep.data.species == 0),2],
#                    col=col[2],pch=18)
# 
#           })
# 
# 
# 
# setClass("BIOMOD.formated.data.PA.indep",
#          contains = "BIOMOD.formated.data.PA",
#          representation(indep.coord = 'data.frame',
#                         indep.data.species = 'numeric',
#                         indep.data.env.var = 'data.frame'),
#          validity = function(object){
#            return(TRUE)
#            })
# 
# BIOMOD.formated.data.PA.indep <- function(sp, env, xy, sp.name,
#                                           PA.NbRep=1,
#                                           PA.strategy='random',
#                                           PA.distance = 0,
#                                           PA.nb.absences = NULL,
#                                           indep.sp, indep.env, indep.xy){
#   BFDP <- BIOMOD.formated.data.PA(sp, env, xy, sp.name,
#                                    PA.NbRep,
#                                    PA.strategy,
#                                    PA.distance,
#                                    PA.nb.absences)
#   BFD <- BIOMOD.formated.data(indep.sp,indep.env,indep.xy, sp.name)
#   
#   BFDPI <- new('BIOMOD.formated.data.PA.indep',
#               sp.name = BFDP@sp.name,
#               coord = BFDP@coord,
#               data.env.var = BFDP@data.env.var,
#               data.species = BFDP@data.species,
#               PA = BFDP@PA,
#               indep.coord = BFD@coord,
#               indep.data.env.var = BFD@data.env.var,
#               indep.data.species = BFD@data.species)
#   rm(BFDP,BFD)
#   return(BFDPI)
# }
# 
# setMethod('plot', signature(x='BIOMOD.formated.data.PA.indep'),
#           function(x,coord=NULL,col=NULL){
#             # coordinates checking
#             if(is.null(coord)){
#               if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
#                 stop("coordinates are required to plot your data")
#               } else {
#                 coord <- x@coord
#               }
#             }
#             
#             # colors checking
#             if(is.null(col) | length(col) < 3){
#               col = c('green', 'red', 'orange', 'grey')
#             }
#             
#             # plot data
#             par(mfrow=c(.CleverCut(ncol(x@PA)+2)))
#             # all points (~mask)
#             plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
#                  main = paste(x@sp.name," Calibration data", sep=""), pch=20 )
#             # presences 
#             points(x=x@coord[which(x@data.species == 1),1],
#                    y=x@coord[which(x@data.species == 1),2],
#                    col=col[1],pch=18)
#             # true absences
#             points(x=x@coord[which(x@data.species == 0),1],
#                    y=x@coord[which(x@data.species == 0),2],
#                    col=col[2],pch=18)
#             
#             # PA data
#             for(i in 1:ncol(x@PA)){
#               # all points (~mask)
#               plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
#                    main = paste(x@sp.name," Pseudo Absences ", i, sep=""), pch=20 )
#               # PA
#               points(x=x@coord[x@PA[,i],1],
#                      y=x@coord[x@PA[,i],2],
#                      col=col[3],pch=18)
#               # presences 
#               points(x=x@coord[which(x@data.species == 1),1],
#                      y=x@coord[which(x@data.species == 1),2],
#                      col=col[1],pch=18)
#               # true absences
#               points(x=x@coord[which(x@data.species == 0),1],
#                      y=x@coord[which(x@data.species == 0),2],
#                      col=col[2],pch=18)
#             }
#             
#             # all points (~mask)
#             plot(x=x@indep.coord[,1], y=x@indep.coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
#                  main = paste(x@sp.name," Validation data", sep=""), pch=20 )
#             # presences 
#             points(x=x@indep.coord[which(x@indep.data.species == 1),1],
#                    y=x@indep.coord[which(x@indep.data.species == 1),2],
#                    col=col[1],pch=18)
#             # true absences
#             points(x=x@indep.coord[which(x@indep.data.species == 0),1],
#                    y=x@indep.coord[which(x@indep.data.species == 0),2],
#                    col=col[2],pch=18)
#             
#           })
# 


####################################################################################################
### Tests ##########################################################################################
####################################################################################################


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
                               family = 'binomial',
                               mustart = 0.5,
                               control = glm.control(maxit = 50)),
                   
                   GBM = list(  distribution = 'bernoulli',
                                interaction.depth = 7,
                                shrinkage = 0.001,
                                bag.fraction = 0.5,
                                train.fraction = 1,
                                n.trees = 500,
                                cv.folds = 5),
                   
                   GAM = list( spline = 2,
                               test = 'AIC',
                               family = 'binomial',
                               control = gam::gam.control(maxit = 50, bf.maxit = 50)),
                   
                   CTA = list(method = 'default',
                              parms = 'default',
                              control = rpart.control(xval = 5, minbucket = 5, minsplit = 5,
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
                             ntree = 50,
                             mtry = 'default'),
                   
                   MAXENT = list(maximumiterations = 200)
                   
                   ))

setMethod('show', signature('BIOMOD.Model.Options'),
          function(object){
            cat("\n-=-=-=- 'BIOMOD.Model.Options' -=-=-=-")
            cat("\n")

            ## GLM options
            cat("\nGLM = list( type = '", object@GLM$type, "',", sep="")
            cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep="")
            cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="")
            cat("\n            test = '", object@GLM$test, "',", sep="")
            cat("\n            family = '", object@GLM$family, "',", sep="")
            cat("\n            mustart = ", object@GLM$mustart, ",", sep="")
            cat("\n            control = glm.control(", .print.control(object@GLM$control), ") ),", sep="")
            
            ## GBM options
            cat("\n")
            cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep="")
            cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep="")
            cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep="")
            cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep="")
            cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep="")
            cat("\n            n.trees = ", object@GBM$n.trees, ",", sep="")
            cat("\n            cv.folds = ", object@GBM$cv.folds, "),", sep="")
            
            ## GAM options
            cat("\n")
            cat("\nGAM = list( spline = ", object@GAM$spline, ",", sep="")
            cat("\n            test = '", object@GAM$test, "',", sep="")
            cat("\n            family = '", object@GAM$family, "',", sep="")
            cat("\n            control = gam::gam.control(", .print.control(object@GAM$control), ") ),", sep="")

            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep="")
            cat("\n            parms = '", object@CTA$parms, "',", sep="")
            cat("\n            cost = ", ifelse(length(object@CTA$cost)<1,'NULL',object@CTA$cost), ",", sep="")
            cat("\n            control = rpart.control(", .print.control(object@CTA$control), ") ),", sep="")
            
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
            cat("\n           mtry = '", object@RF$mtry, "'),", sep="")
                   
            ## MAXENT options
            cat("\n")
            cat("\nMAXENT = list( maximumiterations = ", object@MAXENT$maximumiterations, "),", sep="")

            cat("\n\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
          })

.print.control <- function(ctrl){
  out <- c()
  for (i in 1:length(ctrl)){
    out <- c(out, paste(names(ctrl)[i], " = ", ctrl[[i]], sep=""))
  }
  return(toString(out))
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
         representation(sp.name = 'character',
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
                        models.options = 'BIOMOD.stored.models.options'),
         prototype(sp.name='',
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
                   models.options = new('BIOMOD.stored.models.options')),
         validity = function(object){
           return(TRUE)
           })

setMethod('show', signature('BIOMOD.models.out'),
          function(object){
            cat("\n\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
            cat("\nBIOMOD.models.out")
            cat("\nSpecie modelised :", object@sp.name)
            cat("\nConsidered variables :", object@expl.var.names)
            
            cat("\n\nComputed Models : ", object@models.computed)
            cat("\n\nFailed Models : ", object@models.failed)
            cat("\n\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
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
                  load(obj@models.prediction@link)
                  return(models.prediction)
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
                  load(obj@models.prediction@link)
                  mod.pred <- as.data.frame(models.prediction)
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
                load(obj@models.evaluation@link)
                return(models.evaluation)
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
                load(obj@variables.importances@link)
                return(variables.importances)
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
                load(obj@models.options@link)
                return(models.options)
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
                        rescaled.models = 'logical',
                        type = 'character',
                        proj = 'BIOMOD.stored.data'),
         prototype(proj.names = '',
                   sp.name='',
                   expl.var.names='',
                   models.projected='',
                   rescaled.models=TRUE,
                   type=''),
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
              } else{
                if(obj@proj@link != ''){
                  load(obj@proj@link)
                  return(proj)
                } else{ return(NA) }
              }            
            } else if(as.data.frame){
              if(obj@models.prediction@inMemory ){
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
              if(is.null(str.grep)){
                plot(getProjection(x))
              } else if(length(grep(str.grep, x@models.projected,value=T))>0){
                plot(raster::subset(getProjection(x), grep(str.grep, x@models.projected,value=T)))
              } else{ stop("invalid str.grep arg")}
              
            } else if(class(x@proj) == "BIOMOD.stored.raster.stack"){
              cat("will be able soon!")
            } else {cat("\n !  Biomod Projection plotting issue !")}

          })

setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            cat("\n-=-=-=- 'BIOMOD.projection.out' -=-=-=-")
            cat("\nProjection directory :", paste(object@sp.name,"/",object@proj.names, sep=""))
            cat("\n")
            cat("\nsp.name :", object@sp.name)
            cat("\nexpl.var.names :", object@expl.var.names)
            cat("\n")
            cat("\nmodels projected :", toString(object@models.projected))

            cat("\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
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
            cat("\n-=-=-=- 'BIOMOD.EnsembleModeling.out' -=-=-=-")
            cat("\nsp.name :", object@sp.name)
            cat("\nexpl.var.names :", object@expl.var.names)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed))

            cat("\n-=-=-=--=-=-=--=-=-=--=-=-=--=-=-=--=-\n")
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
          function(obj, model=NULL, meth=NULL){
            if(is.null(model)){
              model <- obj@em.computed
            }
            if(is.character(model) | is.numeric(model)){
              lout <- list()
              for(mod in model){
                if(is.null(meth)){
                  lout[[mod]] <- obj@em.res[[mod]]$em.cross.validation
                } else if(!is.null(meth)){
                  lout[[mod]] <- (obj@em.res[[mod]]$em.cross.validation[meth,,])
                } 
              }
              return(lout)
            } else{
              return(NULL)
            }

          }
          )

setGeneric("getEMbuiltModels",
           function(obj, ...){
             standardGeneric("getEMbuiltModels")
           })

setMethod("getEMbuiltModels", "BIOMOD.EnsembleModeling.out",
          function(obj){
            return(obj@em.computed)
          })

####################################################################################################
### BIOMOD Storing Ensemble Forecasting Objects ####################################################
####################################################################################################
