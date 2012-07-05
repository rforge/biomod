`.scope` <-
function(enviroTrain, Smoother, degree)
{
    XXX <- enviroTrain
    deg <- degree
    vnames <- names(XXX[])
    step.list <- as.list(vnames)
    names(step.list) <- vnames
    NbVar <- dim(enviroTrain)[2]
    i <- 1
    while(i <= NbVar) {
        vname <- names(XXX)[i]
        # loops through independent variable names
        junk <- c(paste("1 + ",vname, sep=""))
        # minimum scope
        if(is.numeric(XXX[,i])) {
            junk <- c(junk, paste(Smoother, "(", vname, ",", deg, ")", sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        else if(is.factor(XXX[,i])) {
            junk <- c(junk, paste(vname, sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        step.list[[vname]] <- junk
        i <- i + 1
    }
    return(step.list)
}

