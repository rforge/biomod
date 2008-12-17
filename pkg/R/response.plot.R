`response.plot` <-
function(mod, Sp){

    Data <- DataBIOMOD[,1:Biomod.material[["NbVar"]]]
    Xp <- as.data.frame(matrix(sapply(Data, mean), nrow(Data), Biomod.material[["NbVar"]], byrow=TRUE, dimnames=list(NULL, colnames(Data))))
    if(!Biomod.material[["algo.choice"]][mod]) stop("The selected model has not been run")
    if(Sp>length(Biomod.material[["species.names"]])) stop(paste("Wrong species index : there was only", length(Biomod.material[["species.names"]]), "species modelled", sep=" ")) 
    
    object <- eval(parse(text=load(paste(getwd(), "/models/", Biomod.material[["species.names"]][Sp], "_", mod, sep=""))))
       
    if(mod=="ANN") if(sum(search()=="package:nnet")==0) library(nnet)
    if(mod=="CTA") if(sum(search()=="package:rpart")==0) library(rpart)
    if(mod=="MARS" | mod=="MDA") if(sum(search()=="package:mda")==0) library(mda)
    if(mod=="RF") if(sum(search()=="package:randomForest")==0) library(randomForest,  verbose=F)
    
    x11()
    if(!is.na(VarImportance)[[1]]) sqnb <- ceiling(sqrt(Biomod.material[["NbVar"]]+1)) else sqnb <- ceiling(sqrt(Biomod.material[["NbVar"]]))
    layout(matrix(c(rep(1,sqnb),2:(sqnb^2+1)),nc=sqnb, byrow=T), widths=rep(1,sqnb), heights=c(0.3,rep(1,ceiling(Biomod.material[["NbVar"]]/sqnb))))
    
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5,y=0.8,pos=1,cex=1.6,labels=paste("Response curves for species", Sp, "  by  ",mod,sep=""),col="#4c57eb")
    par(mar = c(2,2,3.5,1))

    for(i in 1:Biomod.material[["NbVar"]]) {
    
        xr <- sapply(Data, range)
        Xp1 <- Xp
        Xp1[,i] <- seq(xr[1,i], xr[2,i],  len=nrow(Data))

        if(mod=="GLM" | mod=="GAM") Xf <- predict(object, as.data.frame(Xp1), type="response")
        #if(mod=="GBM") Xf <-  predict.gbm(object, as.data.frame(Xp1), object$n.trees, type='response')
        if(mod=="GBM") Xf <-  predict.gbm(object, as.data.frame(Xp1), Models.information[[Sp]]$GBM$best.iter, type="response")
        if(mod=="CTA") Xf <- as.numeric(predict(object, Xp1, type="vector"))
        if(mod=="ANN") Xf <- Rescaler2(as.numeric(predict(object, as.data.frame(Xp1), type="raw")), OriMinMax=Models.information[[Sp]]$ANN$RawPred)
        if(mod=="MARS") Xf <- Rescaler2(as.numeric(predict(object, as.data.frame(Xp1))), OriMinMax=Models.information[[Sp]]$MARS$RawPred)
        if(mod=="MDA") Xf <- Rescaler2(predict(object, as.data.frame(Xp1), type="post")[,2], OriMinMax=Models.information[[Sp]]$MDA$RawPred)
        if(mod=="RF") Xf <- Rescaler2(predict(object, as.data.frame(Xp1), type="prob")[,2], OriMinMax=Models.information[[Sp]]$RF$RawPred)

        #if (length(unique(Xf)) !=1)    # pour ne pas afficher les variables ne présenatant aucune variations  ->  non sélectionnées par le modèle
        plot(Xp1[ ,i], Xf, ylim=c(0,1), xlab="", ylab="", type="l", main=names(Data)[i])     
    } 
    if(!is.na(VarImportance)[[1]]){
        plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE, main="Variable Importance")
        varimp <- t(VarImportance[[Sp]][mod,])
        vec <- ""
        for(i in 1:Biomod.material[["NbVar"]]) vec <- paste(vec, row.names(varimp)[i],"   ", round(varimp[i,1],digits=3), "\n", sep="") 
        text(x=0.3,y=0.4,pos=4,cex=if(Biomod.material[["NbVar"]]<8){1.2} else {1/(Biomod.material[["NbVar"]]/8)},labels=vec)
    }    
    if(mod=="ANN")  detach(package:nnet)
    if(mod=="CTA") detach(package:rpart)
    if(mod=="MARS" | mod=="MDA") detach(package:mda)
    if(mod=="RF") detach(package:randomForest)            
}

