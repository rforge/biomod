`KappaSRE` <-
function(Obs, Fit, TSS=F)
{
if(sum(Obs)==0) stop("\n\n The observed data only contains 0")
Misc <- table(Fit, Obs)
TP <- Misc[4]
TN <- Misc[1]
ca0 <- (TN * 100)/sum(Misc[, 1])
ca1 <- (TP * 100)/sum(Misc[, 2])
if(TSS!=T){
	t <- KappaStat(Misc)
	invisible(list(Kappa=t, se=ca1, sp=ca0))
}
else{
	t <- TSS.Stat(Misc)
	invisible(list(TSS=t, se=ca1, sp=ca0))
}
}

