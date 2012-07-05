`FilteringTransformation` <-
function(ProbData, CutOffdata)
{	
	FUN2 <- function(x,y){
		x[!(x>y)]=0
		if(ncol(x)==1) return(x[,1])
		else return(x)
	}
	return(sweep(data.matrix(ProbData), 2, CutOffdata, FUN2))
}

