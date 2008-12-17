`KappaStat` <-
function(Misc)
{
if(dim(Misc)[1]==1){
	if(row.names(Misc)[1]=="FALSE") Misc <- rbind(Misc, c(0,0))
	else{
		a <- Misc
		Misc <- c(0,0)
		Misc <- rbind(Misc, a)
	}
}
n <- sum(Misc)
n.1 <- sum(Misc[,1])
n.2 <- sum(Misc[,2])
n1. <- sum(Misc[1,])
n2. <- sum(Misc[2,])
Po <- (1/n) * (Misc[1,1] + Misc[2,2])
Pe <- ((1/n)^2) * ((n1. * n.1) + (n2. * n.2))
K <- (Po - Pe)/(1 - Pe)
#cat("\n Kappa=", K, "\n")
return(K)
}

