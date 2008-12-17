`PredictionBestModel` <-
function(ANN=T, CTA=T, GAM=T, GBM=T,GLM=T, MARS=T, MDA=T, RF=T, SRE=T, Bin.trans=T, Filt.trans=T, method='all')
{
    Th <- c('Kappa','TSS','Roc', 'all')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc'")
    
    if(method == 'all') for(k in 1:3) PredictionBestModel(ANN, CTA, GAM, GBM, GLM, MARS, MDA, RF, SRE, Bin.trans, Filt.trans, method=Th[k])
                      
    else { if(Biomod.material[["evaluation.choice"]][method]){
    
        algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
        algo.c[names(which(!Biomod.material[["algo.choice"]]))] <- F

        g <- as.data.frame(matrix(0, nrow=Biomod.material[["NbSpecies"]], ncol=7, dimnames=list(Biomod.material[["species.names"]], c("Best.Model", "Cal", "Eval", "Tot", "CutOff", "Se", "Sp"))))
        gg <- ggg <- gggg <- matrix(0, nrow= dim(DataBIOMOD)[1], ncol=Biomod.material[["NbSpecies"]], dimnames=list(seq(dim(DataBIOMOD)[1]), Biomod.material[["species.names"]]))
    
        #GAM>GLM>GBM>RF>MARS>MDA>CTA>ANN>SRE
        G <- as.data.frame(matrix(0, ncol=9, nrow=Biomod.material[["NbSpecies"]], dimnames=list(seq(Biomod.material[["NbSpecies"]]), c("GAM", "GLM", "GBM", "RF", "MARS", "MDA", "CTA", "ANN", "SRE"))))
    
        i <- 1
        while(i <= Biomod.material[["NbSpecies"]]) {
            if(exists("DataEvalBIOMOD"))   j <- 2   else    j <- 1
    
            eval(parse(text=paste("load('", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"')", sep="")))
            sp.data <- eval(parse(text=paste("Pred_", Biomod.material[["species.names"]][i], sep="")))
    
            for(a in Biomod.material[["algo"]][algo.c]) if(a != 'SRE') eval(parse(text=paste("G$", a, "[i] <- Evaluation.results.", method, "[[i]][a,j]", sep="")))
            temp <- factor(which.max(G[i,]), levels=seq(along=colnames(G)), labels=colnames(G))   
                                                                                              
            for(a in Biomod.material[["algo"]][algo.c]){
                 if(a == temp) {
                      g[i, 2:7] <- eval(parse(text=paste("Evaluation.results.", method, sep="")))[[i]][a,1:6]
                      gg[,i] <- sp.data[,a]
                      if(Bin.trans)  ggg[,i] <- BinaryTransformation(gg[,i], g[i,5])
                      if(Filt.trans) gggg[,i] <- FilteringTransformation(gg[,i], g[i,5])
                 }
            }
            i <- i + 1
        }
        g[,1] <- factor(max.col(G), levels=seq(along=colnames(G)), labels=colnames(G))     
        
        assign(paste("BestModelBy",method,sep=""), g)
        eval(parse(text=paste("save(BestModelBy",method,", file='", getwd(), "/pred/BestModelBy",method,"')", sep="")))
            write.table(g, file=paste(getwd(),"/pred/BestModelBy",method,".txt", sep=""), row.names=F)
        
        assign(paste("PredBestModelBy",method,sep=""), as.data.frame(gg))
        eval(parse(text=paste("save(PredBestModelBy",method,", file='", getwd(), "/pred/PredBestModelBy",method,"')", sep="")))
            write.table(gg, file=paste(getwd(),"/pred/PredBestModelBy",method,".txt", sep=""), row.names=F)
            
        if(Bin.trans) {assign(paste("PredBestModelBy",method,"_Bin",sep=""), as.data.frame(ggg))
        eval(parse(text=paste("save(PredBestModelBy",method,"_Bin, file='", getwd(), "/pred/PredBestModelBy",method,"_Bin')", sep="")))
            write.table(ggg, file=paste(getwd(),"/pred/PredBestModelBy",method,"_Bin.txt", sep=""), row.names=F)}
        
        if(Filt.trans) {assign(paste("PredBestModelBy",method,"_Filt",sep=""), as.data.frame(gggg))
        eval(parse(text=paste("save(PredBestModelBy",method,"_Filt, file='", getwd(), "/pred/PredBestModelBy",method,"_Filt')", sep="")))
            write.table(gggg, file=paste(getwd(),"/pred/PredBestModelBy",method,"_Filt.txt", sep=""), row.names=F)}
    
    }}
}

