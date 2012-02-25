`BinaryTransformation` <-
function(ProbData, CutOffdata)
{
	FUN2 <- function(x,y){
		moa <- apply((x>y),2,as.integer)
		if(ncol(moa)==1) return(moa[,1])
		else return(moa)
	}
	return(sweep(data.matrix(ProbData), 2, CutOffdata, FUN2))
}

