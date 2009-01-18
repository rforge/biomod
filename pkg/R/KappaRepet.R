`KappaRepet` <-
function(Obs, Fit, TSS=FALSE)
{
	if(sum(Obs)==0) stop("\n The observed data only contains 0")
	tab <- as.data.frame(matrix(0, nrow=101, ncol=2))  ### il faut préciser que: "nrow = 101" sinon le data frame se remplit de "NA" au fur et à mesure plutôt que de "0" (et si il y a des NA cela pose problème plus loin).
	
	if(length(unique(Fit))==1){
    Misc<-table(as.vector(Fit) >= as.numeric(unique(Fit)), Obs) ### Robin modified here...(avant ça plantait là)
		if(TSS!=T) a <- KappaStat(Misc)
		else a <- TSS.Stat(Misc)
		TP <- Misc[4]
		TN <- Misc[1]
		ca0 <- (TN * 100)/sum(Misc[,1])
		ca1 <- (TP * 100)/sum(Misc[,2])
		if(is.na(ca0)) ca0<-0
		if(is.na(ca1)) ca1<-0
		if(TSS!=T) invisible(list(Kappa=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
		else invisible(list(TSS=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
	}
	else{
		Quant <- quantile(Fit)
		for(j in 0:100){
			Seuil <- Quant[1] + (j*((Quant[5] - Quant[1])/100))
			Misc<-table(Fit >= Seuil, Obs)
			if(TSS!=T) a <- KappaStat(Misc) else a <- TSS.Stat(Misc)
			if(a != "NA") if(a > 0) {tab[j+1, 1] <- Seuil; tab[j+1, 2] <- a}   ### Question here: Pourquoi ne pas garder les Kappa / TSS < = ??.
			rm(Misc, Seuil)
		}
		
		t <- max(tab[,2],na.rm=T)
		seuil <- tab[tab[,2]==t,1]   ### Note: Ici il se peut qu'on aie plus de 1 seuil...dans ce cas le plus bas est gardé.
		if(t > 0) {
			Misc<-table(Fit >= seuil[1], Obs)
			TP <- Misc[4]
			TN <- Misc[1]
			ca0 <- (TN * 100)/sum(Misc[,1])
			ca1 <- (TP * 100)/sum(Misc[,2])
			if(TSS!=T) invisible(list(Kappa = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
			else invisible(list(TSS = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
		}
		else {
			if(TSS!=T) invisible(list(Kappa = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
			else invisible(list(TSS = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
		}
	}
}

