`CurrentPred` <-
function(GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, MDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc, BinKappa, BinTSS, FiltRoc, FiltKappa, FiltTSS)
{
    algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    algo.c[names(which(!Biomod.material[["algo.choice"]]))] <- F
    
    for(i in 1:Biomod.material[["NbSpecies"]]){
      
        tha <- thk <- tht <- c()    
        for(a in Biomod.material[["algo"]][algo.c]){
            if(Biomod.material[["evaluation.choice"]]["Roc"]) if(a != 'SRE') tha <- c(tha, as.numeric(Evaluation.results.Roc[[i]][a,4]))
            if(Biomod.material[["evaluation.choice"]]["Kappa"]) thk <- c(thk, Evaluation.results.Kappa[[i]][a,4])
            if(Biomod.material[["evaluation.choice"]]["TSS"]) tht <- c(tht, Evaluation.results.TSS[[i]][a,4])
        }
      
        eval(parse(text=paste("load('", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"')", sep="")))
        sp.data <- eval(parse(text=paste("Pred_", Biomod.material[["species.names"]][i], sep="")))

        if(BinRoc && Biomod.material[["evaluation.choice"]]["Roc"]){ 
            mata <- BinaryTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(c(tha,0.5)))
            assign(paste("Pred",Biomod.material[["species.names"]][i],"BinRoc", sep="_"), mata)
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_BinRoc, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_BinRoc')", sep="")))
            write.table(mata, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_BinRoc.txt", sep=""), row.names=F)
        }
        if(BinKappa && Biomod.material[["evaluation.choice"]]["Kappa"]){ 
            matk <- BinaryTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(thk))
            assign(paste("Pred",Biomod.material[["species.names"]][i],"BinKappa", sep="_"), matk)
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_BinKappa, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_BinKappa')", sep="")))
            write.table(matk, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_BinKappa.txt", sep=""), row.names=F)
        }
        if(BinTSS && Biomod.material[["evaluation.choice"]]["TSS"]){
            matt <- BinaryTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(tht))           
            assign(paste("Pred",Biomod.material[["species.names"]][i],"BinTSS", sep="_"), matt) 
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_BinTSS, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_BinTSS')", sep="")))
            write.table(matt, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_BinTSS.txt", sep=""), row.names=F)
        }
        
        if(FiltRoc && Biomod.material[["evaluation.choice"]]["Roc"]){ 
            mata <- FilteringTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(c(tha,0.5)))
            assign(paste("Pred",Biomod.material[["species.names"]][i],"FiltRoc", sep="_"), mata)
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_FiltRoc, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_FiltRoc')", sep="")))
            write.table(mata, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_FiltRoc.txt", sep=""), row.names=F)
        }
        if(FiltKappa && Biomod.material[["evaluation.choice"]]["Kappa"]){ 
            matk <- FilteringTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(thk))
            assign(paste("Pred",Biomod.material[["species.names"]][i],"FiltKappa", sep="_"), matk)
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_FiltKappa, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_FiltKappa')", sep="")))
            write.table(matk, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_FiltKappa.txt", sep=""), row.names=F)
        }
        if(FiltTSS && Biomod.material[["evaluation.choice"]]["TSS"]){
            matt <- FilteringTransformation(as.data.frame(sp.data[,algo.c]), as.data.frame(tht))           
            assign(paste("Pred",Biomod.material[["species.names"]][i],"FiltTSS", sep="_"), matt) 
            eval(parse(text=paste("save(Pred_",Biomod.material[["species.names"]][i],"_FiltTSS, file='", getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i],"_FiltTSS')", sep="")))
            write.table(matt, file=paste(getwd(),"/pred/Pred_",Biomod.material[["species.names"]][i],"_FiltTSS.txt", sep=""), row.names=F)
        }
        
    }
}

