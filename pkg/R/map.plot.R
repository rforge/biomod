`map.plot` <-
function(Sp, repnb=1, dim4=1, coor, method, format.type, wanted, color.gradient='red', Proj.name=NULL){  
    
    Th <- c('Kappa','TSS','Roc')
    mod <- c(Biomod.material[["algo"]], 'all')
    FT <- c('binary', 'filtered', 'probs')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc'")
    if(sum(FT == format.type) == 0) stop("\n : uncorrect format.type name , should be one of 'binary' 'filtered' or 'probs'")
    if(nrow(coor) != nrow(DataBIOMOD)) stop("Uncorrect mapping coordinates : coor and species data are not of the same length")
    if(Sp>length(Biomod.material[["species.names"]])) stop(paste("Wrong species index : there was only", length(Biomod.material[["species.names"]]), "species modelled", sep=" "))
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n") 
    
    if(wanted=="prediction"){
        eval(parse(text=paste("load('", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][Sp],"')", sep="")))
        sp.data <- eval(parse(text=paste("Pred_", Biomod.material[["species.names"]][Sp], sep="")))
    }
    if(wanted=="projection"){
        if(is.null(Proj.name)) stop("you haven't entered a projection name in 'Proj.name' \n") 
        eval(parse(text=paste("load('", getwd(), "/proj.", Proj.name, "/Proj_", Proj.name, "_", Biomod.material[["species.names"]][Sp],"')", sep="")))
        sp.data <- eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material[["species.names"]][Sp], sep="")))
    }

    if(color.gradient=='grey') {
        color.system <- c()
        for(i in seq(93,10,length.out=100)) color.system <- c(color.system, gray(i/100))
        color.system <- c(gray(0.93), color.system, gray(0))
    }
    if(color.gradient=='blue')
        color.system <- c(
        'grey88', 
        rainbow(45, start=0.5, end=0.65), 
        rainbow(10, start=0.65, end=0.7), 
        rainbow(45, start=0.7, end=0.85), 
        'red')
        
    if(color.gradient=='red')  
        color.system <- c(
        'grey88', c(
        rep(c(colors()[c(417,417,515)]), each=5), 
        rev(rainbow(55, start=0.13, end=0.23 )), 
        rev(rainbow(50, start=0.08, end=0.13 )[seq(1,50,length.out=15)]), 
        rev(rainbow(50, end=0.08)[seq(1,50,length.out=15)])), 
        'brown2')


    if(Biomod.material$evaluation.choice[[method]]){
    vec2 <- vec <- ""
    met <- eval(parse(text=paste("Evaluation.results.", method, sep="")))[[Sp]]
    
    for(i in 1:sum(Biomod.material$algo.choice)) {
        vec <- paste(vec, row.names(met)[i], "\n")
        vec2 <- paste(vec2, "  =  ", round(as.numeric(met[i,1]),digits=3), "\n") 
    }                                                                               
    } else {vec2 <- "" ; vec <- "not available" }
    
           
    x11()
    layout(matrix(c(rep(1,6),5:9,2,15:19,2,10:14,3,20:24,4), nr=5, byrow=T), widths=c(1,1,1,1,1,1), heights=c(1,0.5,4,0.5,4))
    par(mar = c(0.1, 0.1, 0.1, 0.1)) 
    pbox <- function(co){ plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE) ; polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col=co,border=NA) }
    pbox("#f5fcba") 
    if(wanted=='prediction') ttext <- paste("Current Prediction for", Biomod.material[["species.names"]][Sp], "in", sep=" ") else ttext <- paste("Projection for", Biomod.material[["species.names"]][Sp], "in", sep=" ") 
    if(format.type=='probs') ttext <- paste(ttext, "probabilities") else ttext <- paste(ttext, format.type, "format by", method, sep=" ") 
    text(x=0.5,y=0.8,pos=1,cex=1.6,labels=ttext,col="#4c57eb")
    pbox("grey98")  
    if(format.type=='binary') {legend(0.2,1,legend=list('0','1'),cex=2,fill=color.system[c(1,102)] ,bty='n')
    } else legend(0.2,0.92,legend=list('1','','','','','0.75','','','','','0.5','','','','','0.25','','','','','0'),cex=1, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')     
    pbox("grey98") 
    pbox("grey98")
    text(x=0.5,y=1,pos=1,cex=1.4,labels="Evaluation",col="#4c57eb")
    text(x=0.5,y=0.92,pos=1,label="by", col="#4c57eb", cex=1.4) 
    text(x=0.5,y=0.84,pos=1,label=method, col="#4c57eb", cex=1.4)      
    text(x=0.06,y=0.35,pos=4,label=vec, col="#4c57eb", cex=1.2)
    text(x=0.36,y=0.35,pos=4,label=vec2, col="#4c57eb", cex=1.2) 

    for(k in 1:10){
        pbox("grey98")
        if(k == 1) text(x=0.5, y=0.8, pos=1, cex=1.2, labels='Data input', col="#4c57eb")
        else text(x=0.5, y=0.8, pos=1, cex=1.2, labels=Biomod.material[["algo"]][Biomod.material[["algo.choice"]]][k-1], col="#4c57eb")
    }                    
                    
    level.plot(DataBIOMOD[,Biomod.material[["NbVar"]]+Sp], CoorXY, show.scale=F, title="")
    for(k in Biomod.material[["algo"]][Biomod.material[["algo.choice"]]]) { 
       if(format.type == 'probs') level.plot(sp.data[,k], XY=coor, color.gradient=color.gradient, show.scale=F, level.range=c(0,1000), title="")
       if(format.type == 'binary') level.plot(BinaryTransformation(sp.data[,k], as.numeric(met[k,4])), XY=coor, color.gradient=color.gradient, show.scale=F, title="")
       if(format.type == 'filtered') level.plot(FilteringTransformation(sp.data[,k], as.numeric(met[k,4])), XY=coor, color.gradient=color.gradient, show.scale=F, title="")
    }
    if(sum(Biomod.material[["algo.choice"]]) != 9) for(k in 1:(9-sum(Biomod.material[["algo.choice"]]))) pbox("grey98")                                 
                    
}

