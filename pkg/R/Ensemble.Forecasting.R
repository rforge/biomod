`Ensemble.Forecasting` <-
function(ANN=T,CTA=T,GAM=T,GBM=T,GLM=T,MARS=T,MDA=T,RF=T,SRE=T, Proj.name, weight.method, decay=1.6, PCA.median=T, binary=T, bin.method='Roc', Test=F)
{
    if(bin.method!='Roc' && bin.method!='Kappa' && bin.method!='TSS') stop("\n bin.method should be one of 'Roc', 'Kappa' or 'TSS'  \n")
    if(weight.method!='Roc' && weight.method!='Kappa' && weight.method!='TSS') stop("\n weight.method should be one of 'Roc', 'Kappa', or 'TSS' \n") 

    Th <- c('Roc', 'Kappa', 'TSS')
    list.out <- vector('list', Biomod.material[["NbSpecies"]])
    names(list.out) <- Biomod.material[["species.names"]]

    G<-gg<-ggg<-gggg<-ggggg<-gggggg <- matrix(nc=Biomod.material[["NbSpecies"]],nr=Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]],dimnames= list(1:Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]], Biomod.material[["species.names"]]))
    thm <- thpond <- c()
    
    ens.choice <- proj.choice
    for(j in Biomod.material[["algo"]]) if(!eval(parse(text=j))) ens.choice[j] <- F
    
    for(i in 1:Biomod.material[["NbSpecies"]]){
        cat(paste(Biomod.material[["species.names"]][i], " \n"))
        out <- list()   # for storing the information on weights awarded, evaluation results, pca.median model selected, etc.

        consensus.mat <- matrix(NA, nr=Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]], nc=6)
        load(paste(getwd(), "/proj.", Proj.name, "/Proj_", Proj.name, "_", Biomod.material[["species.names"]][i], sep=""))      
        sp.data <- eval(parse(text=paste("Proj_", Proj.name, "_", Biomod.material[["species.names"]][i], sep="")))
        consensus.mat[,1] <- apply(sp.data[,ens.choice], 1, mean)
        consensus.mat[,3] <- apply(sp.data[,ens.choice], 1, median)
        
        if(Test){
            consensus.test <- matrix(NA, nr=nrow(DataBIOMOD), nc=6)
            load(paste(getwd(), "/pred/Pred_", Biomod.material[["species.names"]][i], sep=""))      
            test.data <- eval(parse(text=paste("Pred_", Biomod.material[["species.names"]][i], sep="")))        
            consensus.test[,1] <- apply(test.data[,ens.choice], 1, mean)
            consensus.test[,3] <- apply(test.data[,ens.choice], 1, median)
        }
        #--------------- defining the weigths        
        k <- proj.choice
        for(a in Biomod.material[["algo"]]) if(ens.choice[a]) k[a] <- as.numeric(eval(parse(text=paste("Evaluation.results.", weight.method, sep='')))[[i]][a,1]) else k[a] <- 0 
        if(weight.method=='Roc') k['SRE'] <- 0
        
        z <-rep(1,sum(k!=0))  
      	for(j in 2:sum(k!=0)) z[j] <- z[j-1] * decay 
      	z <- c(rep(0,(9-sum(k!=0))) ,z/sum(z))
        W <- rep(0,9) 
        for(m in 1:9) W[m] <- z[sum(k[m]>k)+1]
    
        #--------------- applying weights to projections
        consensus.mat[,2] <- apply((sp.data[,ens.choice]*rep(W[ens.choice], each=Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]])), 1, sum)
        if(Test) consensus.test[,2] <- apply((test.data[,ens.choice]*rep(W[ens.choice], each=nrow(DataBIOMOD))), 1, sum)
              
        if(binary){
            thmi <- thpondi <- c()
            for(a in Biomod.material[["algo"]][Biomod.material[["algo.choice"]]]) {
                thmi    <- c(thmi,    eval(parse(text=paste("Evaluation.results.", bin.method, sep="")))[[Biomod.material[["species.names"]][i]]][a,4])
                thpondi <- c(thpondi, eval(parse(text=paste("Evaluation.results.", weight.method, sep="")))[[Biomod.material[["species.names"]][i]]][a,4])
            }
            thm <- c(thm, mean(as.numeric(thmi), na.rm=T))
            thpondi[is.na(thpondi)] <- 0
            thpond <- c(thpond, sum(as.numeric(thpondi)*W[Biomod.material[["algo.choice"]]]))   
        }

        if(PCA.median){
            if(sum(search()=="package:ade4")==0) library(ade4)  
              
            cons <- dudi.pca(sp.data[,ens.choice], scale=T, scannf = F, nf=2)
            pca.select <- colnames(sp.data[,ens.choice])[which.min(abs(cons$co[,2]))]
            #x11()
            #s.corcircle(cons$co, lab = colnames(sp.data[,ens.choice]), full = FALSE, box = F, sub=Biomod.material[["species.names"]][i])
        }
        
        #doing the methods' binary results means across models 
        for(j in 1:3){
            kdata <- rep(0, Biomod.material[[paste("proj.", Proj.name, ".length", sep="")]]) 
            for(k in Biomod.material[["algo"]][ens.choice]) if(k!='SRE') kdata <- kdata + BinaryTransformation(sp.data[,k], as.numeric(eval(parse(text=paste("Evaluation.results.", Th[j], sep="")))[[Biomod.material[["species.names"]][i]]][k,4]))
            if(SRE) kdata <- kdata + sp.data[,'SRE']/1000
            consensus.mat[,3+j] <- kdata / sum(ens.choice) *1000
            
            if(Test){
                gdata <- rep(0, nrow(DataBIOMOD))
                for(k in Biomod.material[["algo"]][ens.choice]) if(k!='SRE') gdata <- gdata + BinaryTransformation(test.data[,k], as.numeric(eval(parse(text=paste("Evaluation.results.", Th[j], sep="")))[[Biomod.material[["species.names"]][i]]][k,4])) 
                if(SRE) kdata <- kdata + sp.data[,'SRE']/1000
                consensus.test[,3+j] <- gdata / sum(ens.choice) *1000
            }
        }
        
        #test the methods with AUC
        if(Test){
            test <- matrix(nc=1,nr=6,dimnames=list(c('prob.mean','prob.mean.weighted','median','Roc.mean','Kappa.mean','TSS.mean'),"score")) 
            for(j in 1:6) test[j,1] <- somers2(consensus.test[,j], DataBIOMOD[,Biomod.material[["NbVar"]]+i])["C"]
        }  
        
        G[,i] <- consensus.mat[,1]
        gg[,i] <- consensus.mat[,2]
        ggg[,i] <- consensus.mat[,3]
        gggg[,i] <- consensus.mat[,4]
        ggggg[,i] <- consensus.mat[,5]
        gggggg[,i] <- consensus.mat[,6]
       
        out[["stats"]] <- matrix(cbind(if(Test)test else rep(NA,6), if(binary) c(thm[i],thpond[i],NA,500,500,500) else rep(NA,6)), nc=2, dimnames=list(c('prob.mean','prob.mean.weighted','median','Roc.mean','Kappa.mean','TSS.mean'), c("score","threshold"))) 
        out[["weights"]] <- as.data.frame(rbind(Biomod.material[["algo"]], round(W,digits=4)))
        if(PCA.median) out[["PCA.median"]] <- pca.select
        
        list.out[[i]] <- out
    }  
     
    assign(paste("consensus_", Proj.name, "_mean", sep=""), G)
    assign(paste("consensus_", Proj.name, "_mean_weighted", sep=""), gg)
    assign(paste("consensus_", Proj.name, "_median", sep=""), ggg)
    assign(paste("consensus_", Proj.name, "_Roc_mean", sep=""), gggg)
    assign(paste("consensus_", Proj.name, "_Kappa_mean", sep=""), ggggg)
    assign(paste("consensus_", Proj.name, "_TSS_mean", sep=""), gggggg)
    
    eval(parse(text=paste("save(consensus_", Proj.name, "_mean, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_mean')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_mean_weighted, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_mean_weighted')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_median, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_median')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_Roc_mean, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_Roc_mean')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_Kappa_mean, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_Kappa_mean')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_TSS_mean, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_TSS_mean')", sep="")))
    
    
    if(binary){                                                                              
    assign(paste("consensus_", Proj.name, "_mean_bin", sep=""), BinaryTransformation(G,as.data.frame(thm))*1000)
    assign(paste("consensus_", Proj.name, "_mean_weighted_bin", sep=""), BinaryTransformation(gg,as.data.frame(thpond))*1000)
    assign(paste("consensus_", Proj.name, "_Roc_mean_bin", sep=""), BinaryTransformation(gggg,as.data.frame(rep(500,Biomod.material[["NbSpecies"]])))*1000)
    assign(paste("consensus_", Proj.name, "_Kappa_mean_bin", sep=""), BinaryTransformation(ggggg,as.data.frame(rep(500,Biomod.material[["NbSpecies"]])))*1000)
    assign(paste("consensus_", Proj.name, "_TSS_mean_bin", sep=""), BinaryTransformation(gggggg,as.data.frame(rep(500,Biomod.material[["NbSpecies"]])))*1000)
    
    eval(parse(text=paste("save(consensus_", Proj.name, "_mean_bin, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_mean_bin')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_mean_weighted_bin, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_mean_weighted_bin')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_Roc_mean_bin, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_Roc_mean_bin')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_Kappa_mean_bin, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_Kappa_mean_bin')", sep="")))
    eval(parse(text=paste("save(consensus_", Proj.name, "_TSS_mean_bin, file='", getwd(),"/proj.", Proj.name, "/consensus_", Proj.name, "_TSS_mean_bin')", sep="")))
    }
      
    assign(paste("consensus_", Proj.name,"_results", sep=""), list.out, pos=1)
    return(list.out)
}

