`level.plot` <-
function(data.in, XY, color.gradient='red', cex=1, level.range=c(min(data.in),max(data.in)), show.scale=TRUE, title="level plot", save.file="no"){  
    
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n") 
    if(ncol(XY)!=2) stop("\n wrong coordinates given in 'XY' : there should be two columns \n")
    if(nrow(XY)!=length(data.in)) stop("\n data and coordinates should be of the same length \n")

    if(exists("multiple")) multiple.plot <- TRUE  else multiple.plot <- FALSE

    SRC <- F
    SRCvalues <- c(-2,-1,0,1)
    if(length(unique(data.in)) <= 4){
        nb <- 0
        for(i in 1:length(unique(data.in))) if(sum(unique(data.in)[i]==SRCvalues) == 1) nb <- nb+1
        if(nb==length(unique(data.in))) SRC <- T
        #if data made of 0s and 1s (binary)
        if(length(unique(data.in))==2 && sum(match(unique(data.in), SRCvalues))==7) SRC=F
    }

    if(SRC){
        color.system <- c("red", "lightgreen", "grey", "darkgreen") 
        gg <- data.in + 3
        title <- paste("SRC plot ", title, sep="")
    } else {
        if(color.gradient=='grey') {
            color.system <- c()
            for(i in seq(93,10,length.out=100)) color.system <- c(color.system, gray(i/100))
            color.system <- c(gray(0.93), color.system, gray(0))
        }
        if(color.gradient=='blue') {
            color.system <- c('grey88',
            rainbow(45, start=0.5, end=0.65),                       
            rainbow(10, start=0.65, end=0.7),
            rainbow(45, start=0.7, end=0.85),
            'red')
        }
        if(color.gradient=='red') {    
            color.system <- c(
            'grey88',
            c(rep(c(colors()[c(417,417,515)]), each=5),
            rev(rainbow(55, start=0.13, end=0.23 )),
            rev(rainbow(50, start=0.08, end=0.13 )[seq(1,50,length.out=15)]),
            rev(rainbow(50, end=0.08)[seq(1,50,length.out=15)])), 
            'brown2')
        } 
        
        #if range wanted is broader than possible, set to actual range limits
        if(level.range[1]<min(data.in)) level.range[1] <- min(data.in)  
        if(level.range[2]>max(data.in)) level.range[2] <- max(data.in)  
        
        #determine the color code to assess to each value  
        g <- gg <- data.in
        gg[gg <= level.range[1]]  <- level.range[1]
        gg[gg >= level.range[2]] <- level.range[2]
        gg <- gg-min(g) 
        gg <- gg/max(gg)*100 + 1
       
        #over and under-ranged values set to limits of color range
        gg[g < level.range[1]] <- 1
        gg[g > level.range[2]] <- 102
    }    
    
    if(save.file == "pdf") pdf(paste(title, ".pdf", sep=""))
    if(save.file == "jpeg") jpeg(paste(title, ".jpeg", sep=""))
    if(save.file == "tiff") tiff(paste(title, ".tiff", sep=""))      
    if(save.file == "postscript") postscript(paste(title, ".eps", sep=""))
  
    if(show.scale){
        
        if(!multiple.plot) layout(matrix(c(1,2),nr=1), widths=c(5,1), heights=c(1,1))
        plot(XY[,2]~XY[,1], col=color.system[gg], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
        par(mar=c(0.1,0.1,0.1,0.1))
        plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE) 
        polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
        
        if(SRC){ legend(0,0.8,legend=list(' (1) new', ' (0) stable', '(-1) kept', '(-2) lost'),cex=1, fill=rev(color.system),bty='n') 
         } else {
          if(level.range[1] == min(data.in)) lmin <- round(level.range[1], digits=2) else lmin <- paste(round(level.range[1], digits=2), " or lower", sep="")
          if(level.range[2] == max(data.in)) lmax <- round(level.range[2], digits=2) else {lmax <- paste(round(level.range[2], digits=2), " or over", sep="") ; color.system[102] <- "grey70"}

          if(!multiple.plot){
              legend(0.2,0.92,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
              '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin),cex=1, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')
          
          } else legend(0.2,1.05,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
          '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin), cex=legendcex, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')
        
        }
    }
     else plot(XY[,2]~XY[,1], col=color.system[gg], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)  
     
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
     
}

