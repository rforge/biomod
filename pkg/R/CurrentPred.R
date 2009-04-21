`CurrentPred` <-
function(GLM=TRUE, GBM=TRUE, GAM=TRUE, CTA=TRUE, ANN=TRUE, SRE=TRUE, MDA=TRUE, MARS=TRUE, RF=TRUE, BinRoc=FALSE, BinKappa=FALSE, BinTSS=FALSE, FiltRoc=FALSE, FiltKappa=FALSE, FiltTSS=FALSE)
{
   
    if(BinRoc && !Biomod.material$evaluation.choice["Roc"] | FiltRoc && !Biomod.material$evaluation.choice["Roc"]) { BinRoc <- FiltRoc <- F ; cat("Roc cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinKappa && !Biomod.material$evaluation.choice["Kappa"] | FiltKappa && !Biomod.material$evaluation.choice["Kappa"]) { BinKappa <- FiltKappa <- F ; cat("Kappa cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    if(BinTSS && !Biomod.material$evaluation.choice["TSS"] | FiltTSS && !Biomod.material$evaluation.choice["TSS"]) { BinTSS <- FiltTSS <- F ; cat("TSS cannot be used to transform probabilities into binary or filtered values, it was not selected in Models() \n ")}
    
    #determine the models for which the operation is wanted and possible
    algo.c <- c(ANN=ANN, CTA=CTA, GAM=GAM, GBM=GBM, GLM=GLM, MARS=MARS, MDA=MDA, RF=RF, SRE=SRE)
    algo.c[names(which(!Biomod.material$algo.choice))] <- F 
      
    transfor <- c("BinRoc", "BinKappa", "BinTSS", "FiltRoc", "FiltKappa", "FiltTSS")
    
    
    for(i in 1:Biomod.material$NbSpecies){
          
        #load the species data from hard disk
        eval(parse(text=paste("load('", getwd(), "/pred/Pred_", Biomod.material$species.names[i],"')", sep="")))
        sp.data <- eval(parse(text=paste("Pred_", Biomod.material$species.names[i], sep="")))

        for(transfo in transfor){ if(eval(parse(text=transfo))){  #go if only the transformation is wanted and possible, i.e. the evaluation method was selected
            
            sp.data.new <- array(NA, dim(sp.data), dimnames=dimnames(sp.data))
        
            for(j in 1:(dim(sp.data)[4])){  #loop by number of pseudo-absences repetitions
                for(k in 1:(dim(sp.data)[3])){ #loop by number of Datasplit calibration process repetitions
       
                    if(Biomod.material$NbRepPA == 0) nam <- paste(Biomod.material$species.names[i], "_full", sep="") else nam <- paste(Biomod.material$species.names[i], "_PA", j, sep="")
                    if(k!=1) nam <- paste(nam, "_rep", k-1, sep="")
                        
                    th <- c() 
                    if(transfo == "BinRoc" | transfo == "FiltRoc") {   for(a in Biomod.material$algo[algo.c]) if(a != 'SRE') th <- c(th, as.numeric(Evaluation.results.Roc[[nam]][a,4]))  ; th <- c(th, 0.5) }
                    if(transfo == "BinKappa" | transfo == "FiltKappa") for(a in Biomod.material$algo[algo.c]) th <- c(th, as.numeric(Evaluation.results.Kappa[[nam]][a,4]))
                    if(transfo == "BinTSS" | transfo == "FiltTSS")     for(a in Biomod.material$algo[algo.c]) th <- c(th, as.numeric(Evaluation.results.TSS[[nam]][a,4]))
                      
                    if(transfo == "BinRoc" | transfo == "BinKappa" | transfo == "BinTSS") sp.data.new[,algo.c,k,j] <- as.matrix(BinaryTransformation(as.data.frame(sp.data[,algo.c,k,j]), as.data.frame(th)))                      
                    if(transfo == "FiltRoc" | transfo == "FiltKappa" | transfo == "FiltTSS") sp.data.new[,algo.c,k,j] <- as.matrix(FilteringTransformation(as.data.frame(sp.data[,algo.c,k,j]), as.data.frame(th)))

                }
            }

            assign(paste("Pred",Biomod.material$species.names[i], transfo, sep="_"), sp.data.new)
            eval(parse(text=paste("save(Pred_",Biomod.material$species.names[i],"_", transfo, ", file='", getwd(), "/pred/Pred_", Biomod.material$species.names[i],"_", transfo,"')", sep="")))

        }}  #transfo loop        
    } #species loop
    
}
