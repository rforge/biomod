`scopeExpSyst` <-
function(enviroTrain, mod)
{
    XXX <- enviroTrain
    NbVar <- dim(enviroTrain)[2]
    vnames <- names(XXX[])
    i <- 1
    junk2 <- c()
    while(i <= NbVar) {
        vname <- names(XXX)[i]
        junk <- c()
        
        if(mod == "NNET") junk <- c(junk, paste(vname))
        if(mod == "MDA")  junk <- c(junk, paste(vname))
        if(mod == "GLMs") junk <- c(junk, paste(vname))
        if(mod == "CTA")  junk <- c(junk, paste(vname))
        if(mod == "GBM")  junk <- c(junk, paste(vname))
        
        if(mod == "GLMq") {
            if(is.numeric(XXX[,i]))      junk <- c(junk, paste(vname, "+I(", vname, "^2)+I(",vname, "^3)"))
            else if(is.factor(XXX[,i]))  junk <- c(junk, paste(vname))
        }
        if(mod == "GLMp") {
            if(is.numeric(XXX[,i]))     junk <- c(junk, paste(vname, "+I(", vname, "^2)+I(",vname, "^3)+", "poly(", vname, ",2) + poly(", vname, ",3)"))
            else if(is.factor(XXX[,i])) junk <- c(junk, paste(vname))
        }
        junk2 <- c(paste(junk2), paste(junk))
        i <- i + 1
    }
    junk2 <- eval(parse(text=paste("~", paste(junk2, collapse="+"))))
    return(junk2)
}

