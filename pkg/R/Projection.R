`Projection` <-
function(Proj=NULL, Proj.name, GLM=F, GBM=F, GAM=F, CTA=F, ANN=F, SRE=F, Perc025=F, Perc05=F,
MDA=F, MARS=F, RF=F, BinRoc=F, BinKappa=F, BinTSS=F, FiltRoc=F, FiltKappa=F, FiltTSS=F)
{
    require(nnet, quietly=T)                                                                    
    require(rpart, quietly=T)
    require(Hmisc, quietly=T)
    require(Design, quietly=T)
    require(MASS, quietly=T)
    require(gbm, quietly=T)
    require(mda, quietly=T)
    require(randomForest, quietly=T)
    require(gam, quietly=T)	

    dir.create(paste(getwd(), "/proj.", Proj.name, sep=""), showWarnings=F)

    if(BinRoc && !Biomod.material$evaluation.choice["Roc"]) { BinRoc=F; cat("\n Roc cannot be used to transform probabilities into binary values, it was not selected in Models()")}
    if(FiltRoc && !Biomod.material$evaluation.choice["Roc"]) { FiltRoc=F; cat("\n Roc cannot be used to transform probabilities into filtered values, it was not selected in Models()")}   
    if(BinKappa && !Biomod.material$evaluation.choice["Kappa"]) { BinKappa=F; cat("\n Kappa cannot be used to transform probabilities into binary values, it was not selected in Models()")}
    if(FiltKappa && !Biomod.material$evaluation.choice["Kappa"]) { FiltKappa=F; cat("\n Kappa cannot be used to transform probabilities into filtered values, it was not selected in Models()")}
    if(BinTSS && !Biomod.material$evaluation.choice["TSS"]) { BinTSS=F; cat("\n TSS cannot be used to transform probabilities into binary values, it was not selected in Models()")}
    if(FiltTSS && !Biomod.material$evaluation.choice["TSS"]) { FiltTSS=F; cat("\n TSS cannot be used to transform probabilities into filtered values, it was not selected in Models()")}

    Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]] <- nrow(Proj)
    assign("Biomod.material", Biomod.material, pos=1)

    algo.c <- algo.d <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    algo.d['SRE'] <- algo.d['GLM'] <- algo.d['GAM'] <- F
 
    w <- names(which(!Biomod.material[["algo.choice"]][names(which(algo.c))]))
    if(length(w) > 0) cat(paste("\n The following model has not been used to render projections : ", paste(w, sep=" "),"\n it has not been trained \n", sep=""))     
    algo.c[names(which(!Biomod.material[["algo.choice"]]))] <- F
    assign("proj.choice", algo.c, pos=1) 
    
    algo.cc <- algo.c
    algo.cc['SRE'] <- F
    
    i <- 1
    while(i <= Biomod.material[["NbSpecies"]]){ cat(paste(Biomod.material[["species.names"]][i], " \n"))
       
        #create matrices to store projections
        g <- gg <- ggg <- gggg <- k <- kk <- kkk<- matrix(NA, nrow=nrow(Proj), ncol=9, dimnames=list(seq(nrow(Proj)), Biomod.material[["algo"]])) 
    
        for(a in Biomod.material[["algo"]][algo.c]){
            if(a != 'SRE') object <- eval(parse(text=load(paste(getwd(), "/models/", Biomod.material[["species.names"]][i], "_", a, sep=""))))      

            #établir les projections, cas particuliers de GLM et GAM qui remplissent les matrices en même temps
            if(a == 'GLM') {
                if(object$deviance == object$null.deviance) {  algo.cc["GLM"] <- F
                    if((sum(DataBIOMOD[,Biomod.material[["NbVar"]]+i])/nrow(DataBIOMOD)) < 0.5) g[,i] <- rep(0, nrow(g))
                    else g[,a] <- gg[,a] <- ggg[,a] <- gggg[,a] <- k[,a] <- kk[,a] <- kkk[,a] <- rep(1000, nrow(g))  #en fait je pense que ça sert à rien vu que c'est refait à l'étape des transformations
                } else g[,a] <- as.integer(as.numeric(predict.glm(object, Proj, type="response")) *1000)
            }

            if(a == 'GAM') {
                if(object$deviance == object$null.deviance) {  algo.cc["GAM"] <- F
                    if((sum(DataBIOMOD[,Biomod.material[["NbVar"]]+i])/nrow(DataBIOMOD)) < 0.5) g[,i] <- rep(0, nrow(g))
                    else g[,a] <- gg[,a] <- ggg[,a] <- gggg[,a] <- k[,a] <- kk[,a] <- kkk[,a] <- rep(1000, nrow(g)) 
                } else g[,a] <- as.integer(as.numeric(predict.gam(object, Proj, type="response")) *1000)
            }
            
            if(a == 'GBM') g[,a] <- as.integer(as.numeric(predict.gbm(object,Proj,Models.information[[i]]$GBM$best.iter, type='response')) *1000)
            if(a == 'CTA') g[,a] <- as.integer(as.numeric(predict(object, Proj, type="vector")) *1000)
            if(a == 'ANN') g[,a] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="raw")), type="range", OriMinMax=Models.information[[i]]$ANN$RawPred) *1000) 
            if(a == 'SRE') g[,a] <- as.integer(as.numeric(sre(DataBIOMOD[,Biomod.material[["NbVar"]]+i], DataBIOMOD[, 1:Biomod.material[["NbVar"]]], Proj, Perc025, Perc05)) *1000)
            if(a == 'MDA') g[,a] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="post")[,2]), type="range", OriMinMax=Models.information[[i]]$MDA$RawPred) *1000) 
            if(a == 'MARS') g[,a] <- as.integer(Rescaler2(as.numeric(predict(object, Proj)), type="range", OriMinMax=Models.information[[i]]$MARS$RawPred) *1000) 
            if(a == 'RF') g[,a] <- as.integer(Rescaler2(as.numeric(predict(object, Proj, type="prob")[,2]), type="range", OriMinMax=Models.information[[i]]$RF$RawPred) *1000) 
            
            #making the binary and filtered transformations if wanted
            if(algo.cc[a]){
                if(BinRoc) gg[,a] <- as.numeric(BinaryTransformation(g[,a], as.numeric(Evaluation.results.Roc[[i]][a,4])))
                if(FiltRoc) ggg[,a] <- as.numeric(FilteringTransformation(g[,a], as.numeric(Evaluation.results.Roc[[i]][a,4])))
                if(BinKappa) gggg[,a] <- as.numeric(BinaryTransformation(g[,a], Evaluation.results.Kappa[[i]][a,4]))
                if(FiltKappa) k[,a] <- as.numeric(FilteringTransformation(g[,a], Evaluation.results.Kappa[[i]][a,4]))
                if(BinTSS) kk[,a] <- as.numeric(BinaryTransformation(g[,a], Evaluation.results.TSS[[i]][a,4]))
                if(FiltTSS) kkk[,a] <- as.numeric(FilteringTransformation(g[,a], Evaluation.results.TSS[[i]][a,4]))
            }
             ggg[,'SRE']<-k[,'SRE']<-kkk[,'SRE'] <-  g[,'SRE']
             gg[,'SRE']<-gggg[,'SRE']<-kk[,'SRE'] <-  g[,'SRE']/1000
        }     

        #exportation des objets créés dans le dossier de travail            
        assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i], sep="_"), g)
        eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],", file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"')", sep="")))
        write.table(g, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],".txt", sep=""), row.names=F)
        
        if(BinRoc){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"BinRoc", sep="_"), gg)
                   eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinRoc, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinRoc')", sep="")))
                   write.table(gg, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinRoc.txt", sep=""), row.names=F)}  
                    
        if(FiltRoc){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"FiltRoc", sep="_"), ggg)
                    eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltRoc, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltRoc')", sep="")))
                    write.table(ggg, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltRoc.txt", sep=""), row.names=F)}   
                         
        if(BinKappa){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"BinKappa", sep="_"), gggg)
                     eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinKappa, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinKappa')", sep="")))
                     write.table(gggg, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinKappa.txt", sep=""), row.names=F)}
                     
        if(FiltKappa){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"FiltKappa", sep="_"), k)
                      eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltKappa, file='", getwd(),"/proj.", Proj.name, "/Proj_", Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltKappa')", sep="")))
                      write.table(k, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltKappa.txt", sep=""), row.names=F)}
                      
        if(BinTSS){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"BinTSS", sep="_"), kk)
                   eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinTSS, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinTSS')", sep="")))
                   write.table(kk, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_BinTSS.txt", sep=""), row.names=F)}
    
        if(FiltTSS){assign(paste("Proj",Proj.name,Biomod.material[["species.names"]][i],"FiltTSS", sep="_"), kkk)
                    eval(parse(text=paste("save(Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltTSS, file='", getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltTSS')", sep="")))
                    write.table(kkk, file=paste(getwd(),"/proj.", Proj.name, "/Proj_",Proj.name,"_",Biomod.material[["species.names"]][i],"_FiltTSS.txt", sep=""), row.names=F)}
                    
        i <- i+1                                       
    }
}

