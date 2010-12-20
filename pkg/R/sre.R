sre <- function (Response = NULL, Explanatory = NULL, NewData = NULL, Quant = 0.025) {    if (Quant >= 0.5 | Quant < 0)         stop("\n settings in Quant should be a value between 0 and 0.5 ")     quants <- c(0 + Quant, 1 - Quant)    if (class(Explanatory)[1] != "RasterStack") {        Response <- as.data.frame(Response)        NbVar <- ncol(Explanatory)
        if (class(NewData)[1] != "RasterStack") 
            Pred <- as.data.frame(matrix(0, nr = nrow(NewData), 
                nc = ncol(Response), dimnames = list(seq(nrow(NewData)), 
                  colnames(Response))))
        for (i in 1:ncol(Response)) {
            ref <- Explanatory[Response[, i] == 1, ]
            if (class(NewData)[1] == "RasterStack") {
                TF <- subset(NewData, 1)
                TF <- TF >= TF@data@min
            }
            else TF <- rep(1, nrow(NewData))
            
            for (j in 1:NbVar) {
                capQ <- quantile(ref[, j], probs = quants)
                if (class(NewData)[1] != "RasterStack") {
                  TF <- TF * (NewData[, names(ref)[j]] >= capQ[1])
                  TF <- TF * (NewData[, names(ref)[j]] <= capQ[2])
                }
                else {
                  TFmin <- NewData@layers[[which(NewData@layernames == 
                    names(ref)[j])]] >= capQ[1]
                  TFmax <- NewData@layers[[which(NewData@layernames == 
                    names(ref)[j])]] <= capQ[2]
                  TF <- TF * TFmin * TFmax
                }
            }
            if (class(TF)[1] != "RasterLayer") 
                Pred[, i] <- TF
            else Pred <- TF
        }
    }    
    if (class(NewData)[1] != "RasterStack" & ncol(Response) == 1)     	Pred <- Pred[[1]]    return(Pred)}