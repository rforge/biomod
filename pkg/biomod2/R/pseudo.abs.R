`pseudo.absences.sampling` <-
function(sp, env, nb.repet=1, strategy='random', distMin=0, distMax=NULL, nb.points=NULL, quant.SRE = 0){
  
  # 1. Parameters checking
  args <- .check.params.pseudo.absences.sampling(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE)
  
  sp <- args$sp
  env <- args$env
  nb.repet <- args$nb.repet
  strategy <- args$strategy
  distMin <- args$distMin
  distMax <- args$distMax
  nb.points <- args$nb.points
  quant.SRE <- args$quant.SRE
  
  rm("args")
  
  if( nb.repet == 0 | nb.points <= 0){
    out <- NULL
  } else {
    out <- switch(strategy,
                   random = random.pseudo.abs.selection( sp, env, nb.points, nb.repet ),
                   sre = sre.pseudo.abs.selection(sp, env, quant.SRE, nb.points, nb.repet),
                   disk = disk.pseudo.abs.selection(sp, env, distMin, distMax, nb.points, nb.repet))
  }

  return(out)
  
#   # 2. Check if NA are present in sp or not to determine which dataset to use 
#   if(sum(is.na(sp@data)) > 0 ){ # PA will be taken into response variable
#     cat("\n*** PA selection")
#     pa.tab <- switch(strategy,
#                      random = random.pseudo.abs.selection(data=sp, nb.points=nb.points, nb.repet=nb.repet),
#                      sre = sre.pseudo.abs.selection(sp),
#                      disk = disk.pseudo.abs.selection(sp))
#     .arranging.pa.table()
#   } else{ # PA will be taken into explanatory variables
#     if(inherits(env, 'Raster')){ # Raster env var case
#       
#     } else if(inherits(env, 'SpatialPoints')){ # spatial data.frame case
#       
#     }
# 
#   }
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.check.params.pseudo.absences.sampling <- function(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE){
  cat("\n   > Pseudo Absences Selection checkings...")
  
  # define here the implemented strategies
  availableStrategies <- c("random", "sre", "disk")
  
  # 1. sp input checking
  if(is.vector(sp)){
    sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp))
  }
  
  if(!(inherits(sp, 'SpatialPoints'))){
    stop("specie input must be a SpatialPointsDataFrame object")
  }
  
  # 2. env input checking
  if(is.matrix(env) | is.data.frame(env)){
    if(nrow(env) != nrow(sp)){
      stop("Species and Explanatory must have same dimentions")
    }
    env <- SpatialPointsDataFrame(coordinates(sp), as.data.frame(env))
  }
  
  if(!inherits(env, 'SpatialPoints') & !inherits(env, 'Raster')){
    stop("Explanatory variables input must be a SpatialPointsDataFrame or a RasterStack object")
  }
  
  # 3. Strategy checking
  if( ! (strategy %in% c("random", "sre")) ){
    if( ( sum(abs(coordinates(sp))) == 0 ) | !( strategy %in% availableStrategies ) ){ # no coordinates or unknow strategy
      strategy <- "random"
      cat("\n   ! Random strategy was automaticly selected (that can be due to points coordinates lack or unavailable strategy choosen)")
    }
  }
      
  # 4. Nb points checking
  if(is.null(nb.points)){
    stop("You must give the number of pseudo absences you want")
  } else{
    nbTrueAbs <- .get.nb.true.abs(sp)
    if(nbTrueAbs >= nb.points){
      cat("\n    ! There is more 'true absences' than desired pseudo absences. No pseudo absences selection done.")
      nb.points = 0
#       #### Return a flag that tell to function that no PA selected
#       return(NULL)
    } else { 
      nb.points = nb.points - nbTrueAbs
      }
  }
      
  # 4. Nb repetition checking
  
  # 5. Distances checking
  if(!is.null(distMax) & !is.null(distMin)){
    if(distMin >= distMax){
      stop("distMin >= distMax")
    }
  }
  
  # 6. SRE quantil checking
  if(strategy == 'SRE'){
    if( quant.SRE >= 0.5 | quant.SRE <0 ){
      stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
    }
  }
  
  # 7. return checked params
  return(list(sp = sp,
              env = env,
              nb.repet = nb.repet,
              strategy = strategy,
              distMin = distMin,
              distMax = distMax,
              nb.points = nb.points,
              quant.SRE = quant.SRE))
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.get.nb.true.abs <- function(sp){
  if(is.vector(sp)) return(sum(sp==0, na.rm=TRUE)) 
  
  if(inherits(sp, 'SpatialPoints')) return(sum(sp@data==0, na.rm=TRUE))
  
  if(inherits(sp, 'Raster')) return(sum(sp[]==0, na.rm=TRUE))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# .is.some.na.in.data <- function(sp){
#   if(is.vector(sp)){
#     if(sum(is.na(sp)) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }    
#   }
#   
#   if(inherits(sp, 'SpatialPoints')){
#     if(sum(is.na(sp[,1])) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
#   
#   if(inherits(sp, 'Raster')){
#     if(sp@data@min >= 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.nb.available.pa.cells <- function(data){
  if(is.vector(data)){ return(sum(is.na(data))) }
  
  if(is.data.frame(data) | is.matrix(data)){ return(sum(is.na(data))) }

  if(inherits(data, 'SpatialPoints')){ return(sum(is.na(data@data))) }
  
  if(inherits(data, 'Raster')){ return( sum(na.omit(data[]) == -1) )}
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.rand.pseudo.abs.selection <- function(data, nb.points){
  if(is.vector(data)){ return(sample(1:length(data), nb.points, replace=FALSE)) }

  if(inherits(data, 'SpatialPoints')){ return(sample(1:nrow(data@data), nb.points, replace=FALSE))}
  
  if(inherits(data, 'Raster')){ return(sort(sampleRandom(x=data, size=nb.points, cells=T)[,"cell"]))}
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "random.pseudo.abs.selection", 
              def = function(sp,env, ...){
                      standardGeneric( "random.pseudo.abs.selection" )
                      } )
# }

setMethod('random.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function( sp, env, nb.points, nb.repet ){
            cat("\n   > random pseudo absences selection")
            
            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                pa.tab[sample(x=cand.cells,size=nb.points,replace=FALSE),j] <- TRUE
              }
              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.tab))
            } else {
              cat("\nUnsuported case yet!")
              return(NULL)
            }
          })

setMethod('random.pseudo.abs.selection', signature(env="RasterStack"),
          function( sp, env, nb.points, nb.repet ){
            cat("\n   > random pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                pa.tab[sample(x=cand.cells,size=nb.points,replace=FALSE),j] <- TRUE
              }
              env <- as.data.frame(extract(env, coordinates(sp), method='bilinear'))
              
              return(list(xy = coordinates(sp),
                          sp = as.numeric(unlist(sp@data, use.names=FALSE)),
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              # create a mask
              mask <- raster::subset(env,1)
              mask <- reclass(mask, c(-Inf,Inf,-1))
              
              # remove presences and true absences from our raster
              mask[cellFromXY(mask,coordinates(sp))] <- NA
              
              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
              }
              
              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol=nb.repet, nrow=nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                pa.tab.tmp[,j] <- sampleRandom(x=mask, size=nb.points, cells=T)[,"cell"]
              }
              
              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
              }
              
              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask, selected.cells))
              sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))), use.names=FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
                             pa.tab)
              
              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))              
              
            }
          })
  
# random.pseudo.abs.selection <- function(data, nb.points=1000, nb.repet=1){
#   nb.points.max <- .nb.available.pa.cells(data)
#   if(nb.points > nb.points.max){ nb.repet <- 0 } # all available pa cell will be selected
#   
#   if(inherits(data, 'SpatialPoints')){ # convert sp.data.frame into vector
#     data <- as.vector(data@data)
#   }
#   
#   if(is.vector(data)){
#     if( nb.repet < 1 ){
#       pa.tab <- data.frame(matrix(1:length(data),ncol=1, nrow=length(data), dimnames=list(NULL,c('PA.all.abs'))))
#     } else {
#       pa.tab <- sapply(1:nb.repet, function(i){.rand.pseudo.abs.selection(data, nb.points)})
#       colnames(pa.tab) <- paste('PA',1:nb.repet,sep="")
#     }
#   }
#   
#   if(inherits(data, 'Raster')){
#     if( nb.repet < 1 ){
#       pa.tab <- data.frame(matrix(which(data[]==-1),ncol=1, nrow=nb.points.max, dimnames=list(NULL,c('PA.all.abs'))))
#     } else {
#       pa.tab <- sapply(1:nb.repet, function(i){.rand.pseudo.abs.selection(data, nb.points)})
#       colnames(pa.tab) <- paste('PA',1:nb.repet,sep="")
#     }
#   }
#     
#   return(pa.tab)
# 
#   ## return xy, sp, env, pa.tab
#   
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "sre.pseudo.abs.selection", 
              def = function(sp,env, ...){
                      standardGeneric( "sre.pseudo.abs.selection" )
                      } )
# }
setMethod('sre.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")

            # 1. calculate sre to determine availables 
            mask <- sre(Response = sp, Explanatory = env, NewData = env@data, Quant = quant.SRE)
            
            # removing cells in envelops, presences and absences
            mask[mask[] == 0] <- NA
            mask[which(as.vector(sp@data)==1),1] <- 1
            mask[which(as.vector(sp@data)==0),1] <- 0

            # 2. Check if NA are present in sp or not to determine which dataset to use
#             if(.nb.available.pa.cells(mask) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(mask)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(mask))
              for(j in 1:ncol(pa.tab)){
                pa.tab[sample(x=cand.cells,size=nb.points,replace=FALSE),j] <- TRUE
              }
              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.tab))

#             }
          })



setMethod('sre.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")
            
            # 1. calculate sre to determine availables 
            mask <- sre(Response = sp, Explanatory = env, NewData = env, Quant = quant.SRE) 
            
            # removing cells in envelops, presences and absences
            mask[mask[]==1] <- NA
            mask[cellFromXY(mask,coordinates(sp)[which(as.vector(sp@data)==1),])] <- NA
            mask[cellFromXY(mask,coordinates(sp)[which(as.vector(sp@data)==0),])] <- NA
            
            mask <- reclass(mask, c(-Inf,Inf,-1))
            
            # checking of nb candidates
            nb.cells <- .nb.available.pa.cells(mask)
            if(nb.cells <= nb.points){
              nb.repet <- 1
              nb.points <- nb.cells
              cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
            }
            
            # select cells into raster
            pa.tab.tmp <- matrix(NA, ncol=nb.repet, nrow=nb.points)
            for( j in 1:ncol(pa.tab.tmp)){
              pa.tab.tmp[,j] <- sampleRandom(x=mask, size=nb.points, cells=T)[,"cell"]
            }
            
            # puting cells in good format
            selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
            pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
            colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
            for( j in 1:ncol(pa.tab)){
              pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
            }
            
            # puting presences, true absences and pseudo absences together
            xy <- rbind(coordinates(sp)[which(!is.na(as.vector(sp@data))),], xyFromCell(mask, selected.cells))
            sp <- as.numeric(unlist(c(na.omit(as.vector(sp@data)), rep(NA,length(selected.cells))), use.names=FALSE))
            env <- extract(env, xy)
          
            pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
                           pa.tab)
            
            return(list(xy = xy,
                        sp = sp,
                        env = as.data.frame(env),
                        pa.tab = as.data.frame(pa.tab)))              
            
          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setGeneric( "disk.pseudo.abs.selection", 
              def = function(sp,env, ...){
                      standardGeneric( "disk.pseudo.abs.selection" )
                      } )

setMethod('disk.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")
            
            # 1. determining selectable area
            coor <- coordinates(sp)
            pres <- which(sp@data[,1]==1)
            true.abs <- which(sp@data[,1]==0)
            tmp.abs <- which(is.na(sp@data[,1])) #(1:ncol(sp@data))[-c(pres,true.abs)]
            outside <- rep(0, length(abs))
            inside <- rep(0, length(abs))
            

            for(i in 1:length(pres)){
              # removing points too close from presences
              inside <- inside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) > distMin )
              # keeping points not to far from presences
              outside <- outside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) < distMax )
            }
            selected.abs <- tmp.abs[ (inside == length(pres)) & (outside > 0) ]
        			         
            # 2. adding presences and true absences and selecting randomly pseudo absences
            
            return(random.pseudo.abs.selection( sp[c(pres, true.abs, selected.abs),],
                                                env[c(pres, true.abs, selected.abs),],
                                                nb.points, nb.repet ))
            

          })

setMethod('disk.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")
              
            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
                                                data = as.data.frame(extract(env,coordinates(sp),method='bilinear')))
                                                
              return(disk.pseudo.abs.selection(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask
              mask <- maskInside <- maskOutside <- reclass(raster::subset(env,1), c(-Inf,Inf,0))
              pres.xy <- coordinates(sp[which(sp@data[,1]==1),])
              
#               inside <- unique(unlist(extract(mask, pres.xy, buffer = distMin * pointDistance(c(0, 0), c(0, 1), longlat=TRUE), cellnumbers=TRUE)$cells))
#               outside <- unique(unlist(extract(mask, pres.xy, buffer = distMax * pointDistance(c(0, 0), c(0, 1), longlat=TRUE), cellnumbers=TRUE)$cells))
#               cat("\n*** length(inside) = ", length(inside))
#               cat("\n*** length(outside) = ", length(outside))
#               
             
              # to convert longitudinal degrees into metters
              coef.conversion <- ifelse(grep("longlat",env@crs@projargs), 111319.5, 1)
#               coef.conversion <- 1
              ## progress bar
              cat("\n")
              pb <- txtProgressBar(min = 0, max = nrow(pres.xy), initial = 0, char = "=-",width = 20,  style = 3, file = "")
              for(i in 1:nrow(pres.xy)){
                setTxtProgressBar(pb,i)
                maskInside <- maskInside + (distanceFromPoints(mask, pres.xy[i,]) > (distMin * coef.conversion))
                maskOutside <- maskOutside + (distanceFromPoints(mask, pres.xy[i,]) <= (distMax * coef.conversion))
              }
              
              maskInside <- maskInside == nrow(pres.xy)
              maskOutside <- maskOutside > 0
                          
              mask <- maskInside * maskOutside
              mask[mask==0] <- NA
              mask <- (-1) * mask
            
              # remove presences and true absences from our raster
              mask[cellFromXY(mask,coordinates(sp))] <- NA
              
              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All availables cells have been selected (", nb.points, "pseudo absences selected )")
              }
              
              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol=nb.repet, nrow=nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                pa.tab.tmp[,j] <- sampleRandom(x=mask, size=nb.points, cells=T)[,"cell"]
              }
              
              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
              }
              
              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask, selected.cells))
              sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))), use.names=FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
                             pa.tab)
              
              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))              
              
            } 
          })

# 
#     coor <- data.biomod@coord
#     pres <- which(data.biomod@data.species==1)
#     true.abs <- which(data.biomod@data.species==0)
#     abs <- (1:length(data.biomod@data.species))[-c(pres,true.abs)]
# # 
#     	for(i in 1:length(pres))
#   			out <- out + ( sqrt((coor[abs,1]-coor[pres[i],1])^2 + (coor[abs,2]-coor[pres[i],2])^2) > distance)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
# .arranging.pa.table(pa.data, pa.tab, sp.data=NULL, xy=NULL){
# 
#   # transforming sp.data into vector if it's not
#   if(!is.null(sp.data)){ # that means that PA were chosed into explanatories data
#     if(inherits(sp.data, 'SpatialPoints')){
#       xy <- coordinates(sp.data)
#       sp.data <- sp.data@data
#     }
#     if(inherits(sp.data, 'Raster')){
#       xy <- rbind(xyFromCell(sp.data, Which(sp.data >= 1), cells=TRUE), xyFromCell(sp.data, Which(sp.data == 0)))
#       sp.data.tmp <- rep(0,nrow(xy))
#       sp.data.tmp[1:length(Which(sp.data >= 1, cells=TRUE))] <- 1
#       sp.data <- sp.data.tmp
#       rm('sp.data.tmp')
#     }
#   }
#   
#   # getting PA selected
# 
# }

`pseudo.abs_v2` <-
function(data.biomod, nb.repet=1, strategy='random', distance=0, nb.points=NULL)
{	
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# 'pseudo.abs' function permits Pseudo Absences selection. Several strategy for selecting the 
# absences are available ('random', 'per', 'squares', 'circles' or 'sre')
#
#  data.biomod <- BIOMOD.formated.data object
#  nb_repet <- nb of different PA set you want
#  strategy <- one of random, circles, per, sre, squares strategy
#  distance <- the distance min between selectable absences and presence
# Damien Georges, nov. 2011
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  cat("\n-=-=-=- BIOMOD Pseudo-Absences Sampling -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")
  # 1. args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  #   1.1 Validity of given arguments
  
  # check length of all args are tue same
#   if (length(strategy) != length(distance) || (!is.null(nb.points) && length(strategy) != 
#     length(nb.points)) )
#     stop("Please give as many Pseudo Absences sampling strategies as number of Pseudo Absences points
#          to select and distances ")
  
  # if the number of pseudo absences to select is defined, it must be a single integer (same for all
  # species) or a vector of Nb species integer (a different nuber by specie)
#   if( !is.null(nb.points) && (sum(length(nb.points) != c(1, length(data.biomod['species'])) ) == 0 ))
#     stop(paste("If you define a number of Pseudo Absences to be selected, you have to give a lone integer
#          or a vector of ",length(data.biomod['species'])," integer",sep="") )
  
  # Available strategy define ?
	if(strategy!='random' && strategy!='per' && strategy!='squares' && strategy!='circles' && strategy!='sre') stop("\n strategy must be one of random, per, squares, circles, sre \n") 
	
  # 2. Pseudo Absences selection =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  # find the coordinates of possible absences
    coor <- data.biomod@coord
    pres <- which(data.biomod@data.species==1)
    true.abs <- which(data.biomod@data.species==0)
  	abs <- (1:length(data.biomod@data.species))[-c(pres,true.abs)]
  
  cat('\nGiven DataSet summary :\n')
  cat('\tpres = ', length(pres), '\ttrue.abs = ', length(true.abs), '\tposs.abs = ', length(abs),'\n')
  cat('\tstategy = ', strategy, '\tdistance = ', distance, '\tnb.points = ', nb.points, '\n' )
  
  if(!is.null(nb.points)){ # if nb point is less than real absences case
    if(length(true.abs) > nb.points){
      cat('You try to select less absences that true.absenses <- original dataset returned\n')
      tab.out <- data.frame(matrix(FALSE, 
                        nrow = length(data.biomod@data.species),
                        ncol = 1))
      tab.out[c(pres,true.abs),1] <- TRUE
      colnames(tab.out) <- 'all.data'
      return(tab.out)
    }    
  }

  
    out <- rep(FALSE, length(abs))
    
    # running the different strategies
  	if(strategy=='random') abs.set <- abs

  	if(strategy=='per'){
  		abs.set <- abs[coor[abs,1] > max(coor[pres,1]) | coor[abs,1] < min(coor[pres,1]) |
  				   coor[abs,2] > max(coor[pres,2]) | coor[abs,2] < min(coor[pres,2])]
  	}
    
  	if(strategy=='squares'){
  		for(i in 1:length(pres)) {
  			out <- out + (coor[abs,1] > (coor[pres[i],1] + distance) | coor[abs,1] < (coor[pres[i],1] - distance) |
  			     	  coor[abs,2] > (coor[pres[i],2] + distance) | coor[abs,2] < (coor[pres[i],2] - distance))
  		}
  		abs.set <- abs[out==length(pres)]
  	}
    
  	if(strategy=='circles'){
  		for(i in 1:length(pres))
  			out <- out + ( sqrt((coor[abs,1]-coor[pres[i],1])^2 + (coor[abs,2]-coor[pres[i],2])^2) > distance)
  		abs.set <- abs[out==length(pres)]
  	}
    
  	if(strategy == 'sre'){
  		pred <- sre(.allAvailableAbs(data.biomod@data.species), data.biomod@data.env.var, data.biomod@data.env.var)
#       abs.set <- raster::subset(abs, pred[-(1:length(pres))] == 0)
      abs.set <- abs[abs %in% which(pred==0)]
  	}
  	
    # selecting all the possible obsences or....
    abs.selected <- matrix(c(true.abs,abs.set) ,ncol=1)
    colnames(abs.selected) <- 'PA.all.abs'
    nb.abs.pos <- length(c(true.abs,abs.set))
    
  	# ... selecting only a limited number of absences from the whole bank
  	if(!is.null(nb.points) && (nb.points > length(true.abs) )){
      if ((nb.points - length(true.abs) ) < length(abs.set)){
        abs.selected <- c() # initializing selected absences
#         meth.name <- paste(meth.name, "partial", sep=".")
        for (r in 1:nb.repet){ # loop on repetition
          abs.selected <- cbind(abs.selected,c(true.abs,sample(abs.set,(nb.points - length(true.abs)) )) )
        }
        colnames(abs.selected) <- paste('PA',1:nb.repet,sep="")
        nb.abs.pos <- nb.points
      }
    }
        
#   return(abs.selected)
#   print(dim(abs.selected))
#   print(summary(abs.selected))
  tab.out <- data.frame(matrix(FALSE, 
                        nrow = length(data.biomod@data.species),
                        ncol = ncol(abs.selected)))
  colnames(tab.out) <- colnames(abs.selected)
  for (j in 1:ncol(abs.selected)){
    tab.out[abs.selected[,j],j] <- TRUE
  }
  cat('\n', ncol(tab.out), 'Pseudo Absences dataset created (',toString(colnames(tab.out)),')\n')
  cat('\n-=-=-=- Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n')
  return(tab.out)
}

# additional hidden functions
.allAvailableAbs <- function(data.biomod.species){
  out <- data.biomod.species
  if( sum(is.na(out)>0) )
    out[is.na(out)] <- 0
  return(out)
}

.CleverCut <- function(x){
  switch(EXPR=x,
         '1' = return(c(1,1)),
         '2' = return(c(1,2)),
         '3' = return(c(2,2)),
         '4' = return(c(2,2)),
         '5' = return(c(2,3)),
         '6' = return(c(2,3)),
         '7' = return(c(3,3)),
         '8' = return(c(3,3)),
         return(c(3,3)))
}
