`level.plot` <-
function(data.in, XY, color.gradient='red', cex=1, level.range=c(min(data.in),max(data.in)), show.scale=TRUE, title="level plot"){  
    
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n") 
    if(ncol(XY)!=2) stop("\n wrong coordinates given in 'XY' : there should be two columns \n")
    if(nrow(XY)!=length(data.in)) stop("\n data and coordinates should be of the same length \n")

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
      
      
    #SRCvalues <- c(-2,-1,0,1)
    #if(sum(match(SRCvalues, unique(data.in))) == length(data.in)){ SRC plot}  
      
      
    g <- gg <- data.in
    gg <- gg-level.range[1] 
    gg <- gg/level.range[2] *100+1
    gg[g<=level.range[1]] <- 1
    gg[g>=level.range[2]] <- 102
    
    if(show.scale){
        layout(matrix(c(1,2),nr=1), widths=c(5,1), heights=c(1,1))
        plot(XY[,2]~XY[,1], col=color.system[gg], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
        par(mar=c(0.1,0.1,0.1,0.1))
        plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE) 
        polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
        legend(0.2,0.92,legend=list(round(level.range[2], digits=2),'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
         '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',round(level.range[1], digits=2)),cex=1, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')
    }
     else plot(XY[,2]~XY[,1], col=color.system[gg], pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)  
}

