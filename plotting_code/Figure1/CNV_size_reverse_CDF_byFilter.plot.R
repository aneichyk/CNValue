#!/usr/bin/env R

#rCNV Map Project
#Spring 2017
#Talkowski Lab & Collaborators

#Copyright (c) 2017 Ryan Collins
#Distributed under terms of the MIT License

#Code to generate reverse CDFs for CNV size

#####Set parameters
WRKDIR <- "/Users/collins/Desktop/RCollins/Talkowski_Local/CNV_DB/rCNV_map/"
options(scipen=1000,stringsAsFactors=F)
logvect.all <- log10(as.vector(sapply(4:7,function(x){return((1:9)*10^x)})))
logvect.mids <- log10(as.vector(sapply(4:6,function(x){return(c(5,10)*10^x)})))
logvect.mains <- log10(as.vector(sapply(4:6,function(x){return(10^x)})))
pctvect.all <- log10(as.vector(sapply(-1:2,function(x){return((1:9)*10^x)})))
pctvect.mids <- log10(as.vector(sapply(-1:1,function(x){return(c(5,10)*10^x)})))
pctvect.mains <- log10(as.vector(sapply(-1:2,function(x){return(10^x)})))
cols.groups <- c("#A5A6A7","#00BFF4","#EC008D","#FFCB00")

#####Master function
plotsizes <- function(filt,cols,xaxis=T,yaxis=T,
                      mar=c(3.5,3.5,0.5,0.5)){
  #Read data & convert to log-scaled CDF
  sizes <- lapply(list("CTRL","NEURO","SOMA","CNCR"),function(group){
    dat <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.",group,
                            ".noMaxSize.E2.",filt,".txt",sep=""))[,1]
    dat <- log10(dat)
    cdf <- sapply(logvect.all,function(min){
      length(which(dat>=min))/length(dat)
    })
    cdf <- log10(100*cdf)
    cdf[which(!is.finite(cdf))] <- -100
    return(cdf)
  })

  #Read data & compute medians
  medians <- lapply(list("CTRL","NEURO","SOMA","CNCR"),function(group){
    dat <- read.table(paste(WRKDIR,"plot_data/figure1/CNV_size.",group,
                            ".noMaxSize.E2.",filt,".txt",sep=""))[,1]
    dat <- log10(dat)
    return(median(dat))
  })
  medians <- unlist(medians)

  #Set parameters
  par(mar=mar,bty="n")

  #Prepare plot
  plot(x=log10(c(50000,5000000)),y=log10(c(0.03,100)),
       type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")

  #Gridlines - vertical & horizontal interleaved
  abline(v=logvect.all,col="gray95")
  abline(h=pctvect.all,col="gray95")
  abline(v=logvect.mids,col="gray85")
  abline(h=pctvect.mids,col="gray85")
  abline(v=logvect.mains,col="gray75")
  abline(h=pctvect.mains,col="gray75")

  #X-axis
  axis(1,logvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
  axis(1,logvect.mids,labels=NA,tck=-0.02,col="gray20",lwd=1)
  axis(1,logvect.mains,labels=NA,tck=-0.025,col="black",lwd=1.5)
  if(xaxis==T){
    axis(1,at=logvect.mids,line=-0.6,tick=F,las=2,cex.axis=0.7,
         labels=c("50kb","100kb","500kb","1Mb","5Mb","10Mb"))
  }

  #Y-axis
  axis(2,pctvect.all,labels=NA,col="gray50",tck=-0.015,lwd=0.5)
  axis(2,pctvect.mids,labels=NA,tck=-0.02,col="gray20",lwd=1)
  axis(2,pctvect.mains,labels=NA,tck=-0.025,col="black",lwd=1.5)
  if(yaxis==T){
    axis(2,at=c(log10(0.1),pctvect.mids),line=-0.5,tick=F,las=2,
         labels=paste(c(0.1,0.5,1,5,10,50,100),"%",sep=""))
  }

  #Plot CDFs
  sapply(1:length(sizes),function(i){
    points(logvect.all,sizes[[i]],type="l",col=cols[i],lwd=3)
  })

  #Medians
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1),
       col="white",border="white")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1),ytop=par("usr")[4])
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1)-0.1,col="white")
  segments(x0=logvect.all,x1=logvect.all,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray95")
  segments(x0=logvect.mids,x1=logvect.mids,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray85")
  segments(x0=logvect.mains,x1=logvect.mains,
           y0=par("usr")[3],y1=log10(0.1)-0.1,col="gray75")
  points(x=medians,y=rep(mean(c(par("usr")[3],log10(0.1)-0.1)),4),
         pch=18,col=cols,cex=2)
  points(x=medians,y=rep(mean(c(par("usr")[3],log10(0.1)-0.1)),4),
         pch=18,col="black",cex=0.5)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1)-0.1,ytop=log10(0.1),
       col="white",border="white")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=log10(0.1)-0.1,col=NA)
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=log10(0.1),ytop=par("usr")[4],col=NA)
}

#Run plots
pdf(paste(WRKDIR,"rCNV_map_paper/Figures/Figure1/CNVsize_by_filter_by_group.reverse_CDF.pdf",sep=""),
    width=7,height=1.8)
par(mfrow=c(1,4))
plotsizes(filt="coding",cols=cols.groups)
plotsizes(filt="haplosufficient",cols=cols.groups)
plotsizes(filt="noncoding",cols=cols.groups)
plotsizes(filt="intergenic",cols=cols.groups)
dev.off()


