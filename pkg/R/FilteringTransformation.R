`FilteringTransformation` <-
function(ProbData, CutOffdata)
{
    if(is.data.frame(ProbData)) {
        N <- dim(ProbData)[2]
        i <- 1
        while(i <= N) {
            if(sum(ProbData[,i])!=0) ProbData[ProbData[,i] < CutOffdata[i, 1],i] <- 0
            i <- i + 1
        }
    }
    else if(sum(ProbData) != 0) ProbData[ProbData < CutOffdata] <- 0
    
    return(ProbData)
}

