`TSS.Stat` <-
function(Misc)
{
	if(dim(Misc)[1]==1){
		if(row.names(Misc)[1]=="FALSE") Misc<-rbind(Misc, c(0,0))
		else {
			a<-Misc
			Misc<-c(0,0)
			Misc<-rbind(Misc, a)

		}
	}
	n <- sum(Misc)
	a <- Misc[1,1]
	b <- Misc[1,2]
	c <- Misc[2,1]
	d <- Misc[2,2]
	sens<-a/(a+c)
	spec<-d/(b+d)
	K <- (sens + spec) - 1        #TSS
	return(K)
}

