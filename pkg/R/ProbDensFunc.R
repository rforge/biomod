`ProbDensFunc` <-
function(initial, projections, plothist=TRUE, cvsn=TRUE, groups=NULL, resolution=5, save.pdf=FALSE, name="ProbDF_plot"){
    
    if(!is.null(groups) && !is.matrix(groups)) stop("\n 'groups' should be a matrix \n")
    if(!is.null(groups) && ncol(groups)!=ncol(projections)) stop("\n 'groups' and 'projections' do not have the same number of columns \n")
    
    #enabeling NAs in data (typical in BIOMOD if some models are switched off)
    groups <- groups[,which(projections[1,]!="NA")]
    projections <- projections[,which(projections[1,]!="NA")]
  	
  	#area stores the species range change calculations
    area <- (apply(projections,2,sum) / sum(initial==1) -1) * 100
  	a <- round( (min(area, na.rm=T)-(resolution+10))/10 ) *10
  	b <- round( (max(area, na.rm=T)+(resolution+10))/10 ) *10
  	p <- hist(area, breaks = seq(a,b,resolution), plot=F) 
    p$density <- p$counts / sum(p$counts)
  	
  	
    #analysis of the distribution density and calculation of the probability of events 
  	a.s <- sort(area)
		lim <- c(0.5, 0.75, 0.90, 0.95)
		nb <- round(ncol(projections)*lim)
		for(j in 1:4) {
  			g <- rep(NA,ncol(projections)-nb[j])
  
  			for(i in 1:length(g)) g[i] <- diff(range(a.s[i:(i+nb[j])]))
  			
        assign(paste("low", lim[j]*100, sep=''),a.s[which.min(g)])
  			assign(paste("high", lim[j]*100, sep=''),a.s[which.min(g)]+g[which.min(g)])
		}
		#return list
  	out <- list()
    out[["stats"]] <- matrix(c(low50,low75,low90,low95, high50,high75,high90,high95), nc=2, nr=4, dimnames=list(c("50%","75%","90%","95%"), c("lower limit", "upper limit")))    		


    if(save.pdf) pdf(paste(name, ".pdf", sep=""))

    if(!is.null(groups)){
    		if(!save.pdf) x11()
    		par(mfrow=c(1,nrow(groups)))
    		color.samp <- list()                                          
    		for(pa in 1:nrow(groups)){
    
            lv <- levels(as.factor(as.matrix(groups[pa,])))
      			color.samp[[pa]] <- colors()[sample(c(90,417,552,616,382,11,150,468,28,31,420,476,333),length(lv))]
      			g <- hist(area, breaks = seq(a,b,resolution), plot=F)
      			fac <- (max(g$counts) / sum(g$counts)) / max(g$density)
      			g$density <- g$density * fac
      
      			plot(g, freq=F, border='grey88', main=paste('groups',pa,sep=" "), xlab="Species range change (%)", ylab="Event   occurence   probability")
      			lines(density(area, width=30)$x, density(area, width=30)$y*fac)
      			for(i in 1:length(lv)){
        				div <- length(area) / length(area[groups[pa,]==lv[i]])
        				lines(density(area[groups[pa,]==lv[i]],width=30)$x, density(area[groups[pa,]==lv[i]], width=30)$y /div*fac, col=color.samp[[pa]][i])
    			  }
    			  lv <- as.factor(as.matrix(groups[pa,]))
    			  leg <- list()
            for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j]
    			  legend(( if(low95-a > b-high90){low90-40} else{high90+5}),max(p$counts)/sum(p$counts)*1.05, legend=leg, bty='n',fill=color.samp[[pa]])
  		  }
    }
    
    if(cvsn){
        #calculation of the 2 axes independently (lost vs new sites)
        area2 <- (apply(projections[which(initial==1),],2,sum) / sum(initial==1) -1) * 100
    		area3 <- area - area2
    
        if(!save.pdf) x11()
        par(mfrow=c(1,nrow(groups))) 
        
        for(i in 1:nrow(groups)){
      			lv <- as.factor(as.matrix(groups[i,]))
      			leg <- list() 
            for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j] 
      			
            levels(lv) <- 1:length(levels(lv))
      			plot(area3~area2, xlab='current', ylab='new', ylim=c(0,if(max(area3)<100){100}else{max(area3)+30}), xlim=c(-100,0), col=color.samp[[i]][lv], pch=20)
      			legend(-100, if(max(area3)<100){100}else{max(area3)+30}, legend=leg, bty='n',fill=color.samp[[i]])
      			abline(0,-1, col='grey80') 
            abline(100,-1, col='grey80')
      			text(x=-97,y=103,pos=1,label="SRC = 0", col="black", cex=0.8)
            if(max(area3)<100)text(x=-3,y=103,pos=1,label="SRC = 100", col="black", cex=0.8) 
  		      else text(x=if(max(area3)<200){-(max(area3)+33-100)}else{-96},y=if(max(area3)<200){max(area3)+33}else{203},pos=1,label="SRC = 100", col="black", cex=0.8) 
        }
    }
    
    #if plot of distribution plot wanted
    if(plothist){
    		if(!save.pdf) x11()
        par(mfrow=c(1,1))
    		
  			hist( (low95+high95)/2, breaks=c(low95,high95), col="aliceblue", xlim=c(a,b), ylim=c(0,max(p$density)*1.2), xlab="",ylab="", main="")
  			hist( (low90+high90)/2, breaks=c(low90,high90), col="slategray1", xlim=c(a,b), ylim=c(0,max(p$density)*1.2), xlab="",ylab="", main="", add=T)
  			hist( (low75+high75)/2, breaks=c(low75,high75), col="steelblue1", xlim=c(a,b), ylim=c(0,max(p$density)*1.2), xlab="",ylab="", main="", add=T)
  			hist( (low50+high50)/2, breaks=c(low50,high50), col="dodgerblue1", xlim=c(a,b), ylim=c(0,max(p$density)*1.2), xlab="",ylab="", main="", add=T)
  			abline(v=c(low95,high95, low90,high90, low75,high75, low50,high50),col= c(rep("aliceblue",2),rep("slategray1",2),rep("steelblue1",2),rep("dodgerblue1",2)), lwd=1.7)
  			legend(( if(low95-a > b-high90){low90-40} else{high90+5}),max(p$counts)/sum(p$counts)*1.05,
  			legend=list('95%','90%','75%','50%'), bty='n', fill=c("aliceblue","slategray1","steelblue1","dodgerblue1"), cex=0.8, title='distrib. of data' )
  			
        par(new=T)
        plot(p, freq=F, col="white", xlim=c(a,b), ylim=c(0,max(p$density)*1.2), main="Probability density function", xlab="Species range change (%)", ylab="Event   occurence   probability")
  	}
  	
  	if(save.pdf) dev.off()
  	return(out)
}

