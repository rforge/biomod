`scopeGAM` <-
function(enviroTrain)
{
    NbVar <- dim(enviroTrain)[2]
    vnames <- names(enviroTrain[])
    i <- 1
    junk2 <- c()
    while(i <= NbVar) {
        vname <- names(enviroTrain)[i]
        junk <- c()
        junk <- c(junk, paste("s(",vname,",bs='ts')"))        
        junk2 <- c(paste(junk2), paste(junk))
        i <- i + 1
    }
    junk2 <- eval(parse(text=paste("~", paste(junk2, collapse="+"))))
    return(junk2)
}

