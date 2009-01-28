`pseudo.abs` <-
function(coor=NULL, status, env=NULL, strategy='random', distance=0, nb.points=NULL, add.pres=TRUE, 
		    species.name= 'SpNoName', create.dataset=FALSE, plot=FALSE, acol='grey80', pcol='red')

{	
	if(strategy=='sre' && is.null(env)) stop("\n you must enter some environmental data to use the sre strategy \n")
	if(strategy!='sre' && is.null(coor)) stop("\n you must enter coordinates for this strategy \n")
	if(plot && is.null(coor)) stop("\n you must enter coordinates for plotting \n")
	
  nam <- paste(species.name, strategy, sep='.')
	pres <- which(status==1)
	abs <- (1:length(status))[-pres]

	out <- rep(F, length(abs))

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
		nam <- paste(nam, distance, sep='.')
	}
	if(strategy=='circles'){
		for(i in 1:length(pres))
			out <- out + ( sqrt((coor[abs,1]-coor[pres[i],1])^2 + (coor[abs,2]-coor[pres[i],2])^2) > distance)
		abs.set <- abs[out==length(pres)]
		nam <- paste(nam, distance, sep='.')
	}
	if(strategy == 'sre'){
		assign("NbVar", ncol(env), pos=1)
		pred <- sre(status, env, env)
		abs.set <- subset(abs, pred[-(1:length(pres))] == 0)
	}
	if(!is.null(nb.points)) {
		abs.set <- sample(abs.set,nb.points)
		nam <- paste(nam, "partial", sep=".")
      }

	if(add.pres) out.set <- c(pres, abs.set)
	if(create.dataset) if(is.null(coor)) { assign(paste("Dataset",nam,sep="."), status[out.set], pos=1) 
                     } else assign(paste("Dataset",nam,sep="."), cbind(coor[out.set,],status[out.set]), pos=1)
	assign(nam, out.set, pos=1)

	if(plot){
            x11()
		plot(coor[abs.set,2]~coor[abs.set,1], xlim=c(min(coor[,1]), max(coor[,1])), ylim=c(min(coor[,2]), max(coor[,2])), col=acol, xlab="", ylab="", main=nam, xaxt='n', yaxt='n')
		par(new=T);plot(coor[pres,2]~coor[pres,1], col=pcol, ylim=c(min(coor[,2]), max(coor[,2])), xlim=c(min(coor[,1]), max(coor[,1])), xlab="", ylab="", xaxt='n', yaxt='n')
	}
}

