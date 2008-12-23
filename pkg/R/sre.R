`sre` <-
function(Response=NULL, Explanatory=NULL, NewData=NULL, Perc025=F, Perc05=F)
{
    NbVar <- dim(Explanatory)[2]
    #assign("NbVar", NbVar, where=1)
    Var.Names <- colnames(Explanatory)

  if(is.numeric(Response)) {
    if(Perc025!=F | Perc05!=F){
      temp <- cbind(Response, Explanatory)
      temp <- temp[temp[,1]==1,]
      temp2 <- temp
      i <- 2
      while(i<=ncol(temp)){
        Q <- quantile(temp[,i],  probs=c(.025,.05,0.95, 0.975))
        if(Perc025==T) temp2 <- temp2[temp2[,i]>=Q[1] & temp2[,i]<=Q[4],]
        if(Perc05==T) temp2 <- temp2[temp2[,i]>=Q[2] & temp2[,i]<=Q[3],]
        i <- i+1
      }
      Data <- temp2[,2:ncol(temp2)]
    }
    else Data <- Explanatory[Response == 1,  ]
    if(Perc025!=F & Perc05!=F) stop("\n Select only percentile option at one time. Either Perc025=T & Perc05=F, OR Perc025=F & Perc05=T OR Perc025=F & Perc05=F \n")

        prediction <- matrix(0, nrow=dim(NewData)[1], ncol=2*NbVar)
        a <- 1
        j <- 1
        while(j <= NbVar) {
            prediction[, a] <- eval(parse(text=paste("NewData$",paste(Var.Names[j]), collapse=""))) >= min(eval(parse(text=paste("Data$", paste(Var.Names[j]), collapse=""))))
            prediction[, a + 1] <- eval(parse(text=paste("NewData$", paste(Var.Names[j]), collapse=""))) <= max(eval(parse(text=paste("Data$",paste(Var.Names[j]), collapse=""))))
            j <- j + 1
            a <- a + 2
        }
        Pred <- apply(prediction[, 1:(2 * NbVar)], 1, FUN = function(x)
        { sum(x)/(2 * NbVar)
        }
        )
        Pred <- trunc(Pred)
        return(Pred)
    }
    else {
        NbSp <- dim(Response)[2]
        Pred <- as.data.frame(matrix(0, nrow=dim(NewData)[1], ncol=NbSp, dimnames=list(seq(dim(NewData)[1]), colnames(Response))))
        z <- 1
        while(z <= NbSp) {

            if(Perc025==T | Perc05==T){
        temp <- cbind(Response[,z], Explanatory)
        temp <- temp[temp[,1]==1,]
        temp2 <- temp
        i <- 2
        while(i<=ncol(temp)){
          Q <- quantile(temp[,i], probs=c(.025,.05,0.95, 0.975))
          if(Perc025==T) temp2 <- temp2[temp2[,i]>=Q[1] & temp2[,i]<=Q[4],]
          if(Perc05==T) temp2 <- temp2[temp2[,i]>=Q[2] & temp2[,i]<=Q[3],]
          i <- i+1
        }
        Data <- temp2[,2:ncol(temp2)]
      }
      else Data <- Explanatory[Response[,z] == 1,]
            prediction <- matrix(0, nrow=dim(NewData)[1], ncol =2 * NbVar)
            a <- 1
            j <- 1
            while(j <= NbVar) {
                prediction[,a] <- eval(parse(text=paste("NewData$", paste(Var.Names[j]),collapse=""))) >= min(eval(parse(text=paste("Data$", paste(Var.Names[j]), collapse=""))))
                prediction[,a+1] <- eval(parse(text=paste("NewData$", paste(Var.Names[j]),collapse=""))) <= max(eval(parse(text=paste("Data$", paste(Var.Names[j]), collapse=""))))
                j <- j + 1
                a <- a + 2
            }
            Pred[,z] <- apply(prediction[,1:(2*NbVar)], 1,FUN = function(x)
            { sum(x)/(2 * NbVar)
            }
            )
            Pred[,z] <- trunc(Pred[,z])
            z <- z + 1
        }
        return(Pred)
    }
}

