`testnull` <-
function(object, Prev, dat){

    if(object$deviance == object$null.deviance){
        if(Prev < 0.5) pred <- rep(0, nrow(DataEvalBIOMOD))
        if(Prev >= 0.5) pred <- rep(1, nrow(DataEvalBIOMOD))
    }
    else pred <- predict(object, dat, type="response")    
    return(pred)
}

