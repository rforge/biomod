`ProjectionBestModel` <-
function(Proj.name, Bin.trans=TRUE, Filt.trans=TRUE, method='all')
{
    Th <- c('Kappa','TSS','Roc','all')
    if(sum(Th == method) == 0) stop("\n : uncorrect method name , should be one of 'Kappa' 'TSS' 'Roc' 'all'")
    
    if(method == 'all') for(k in 1:3) ProjectionBestModel(Proj.name, Bin.trans, Filt.trans, method=Th[k])
    
    else { if(Biomod.material[["evaluation.choice"]][method]) {
    
        g <- gg <- ggg <-  as.data.frame(matrix(0, nrow=Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]], ncol=Biomod.material[["NbSpecies"]], dimnames= list(seq(Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]]), Biomod.material[["species.names"]])))

        i <- 1
        while(i <= Biomod.material[["NbSpecies"]]) {
            eval(parse(text=paste("load('pred/BestModelBy", method, "')", sep="")))
            j <- as.character(eval(parse(text=paste("BestModelBy",method,"[i,1]",sep=""))))
            load(paste(getwd(), "/proj.", Proj.name, "/Proj_",Proj.name, "_", Biomod.material[["species.names"]][i],sep=''))
            g[,i] <- eval(parse(text=paste("Proj",Proj.name,Biomod.material[["species.names"]][i],sep='_')))[,j]
            if(Bin.trans)    eval(parse(text=paste("gg[,i] <- BinaryTransformation(g[,i], as.numeric(Evaluation.results.", method, "[[i]][j,4]))",sep='')))
            if(Filt.trans)   eval(parse(text=paste("ggg[,i] <- FilteringTransformation(g[,i], as.numeric(Evaluation.results.", method, "[[i]][j,4]))",sep='')))
            i <- i + 1
        }
             
        assign(paste("Proj_",Proj.name,"_BestModelBy",method,sep=""), as.data.frame(g))
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,", file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"')", sep="")))
            write.table(g, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,".txt", sep=""), row.names=F)
            
        if(Bin.trans) {assign(paste("Proj_",Proj.name,"_BestModelBy",method,"_Bin",sep=""), as.data.frame(gg))
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,"_Bin, file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Bin')", sep="")))
            write.table(gg, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Bin.txt", sep=""), row.names=F)}
        
        if(Filt.trans) {assign(paste("Proj_",Proj.name,"_BestModelBy",method,"_Filt",sep=""), as.data.frame(ggg))
        eval(parse(text=paste("save(Proj_",Proj.name,"_BestModelBy",method,"_Filt, file='", getwd(), "/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Filt')", sep="")))
            write.table(ggg, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_BestModelBy",method,"_Filt.txt", sep=""), row.names=F)}
        
    }}
}

