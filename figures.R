#---------------------------------------------------------------------------------
# Levels and color palettes definitions
#---------------------------------------------------------------------------------
# Levels for p1 (in percentile levels) (the 2nd one is zoomed at the end)
BRKS0=c(0,50,70,80,90,93,95,97,98,99,99.3,99.5,99.7,99.8,99.9,99.93,99.95,100)
BRKS1=c(0,50,80,90,95,98,99,99.5,99.8,99.9,99.95,99.98,99.99,99.995,99.998,99.999,99.9995,100)
BRKS2=c(BRKS1[-c(3:5,length(BRKS1))],99.9996,99.9997,99.9998,100)

# Levels for FAR (idem)
BRKS3=c(-1000,0,10,20,33,50,66,75,80,82.5,85,87.5,90,92.5,95,97.5,99,100)
BRKS4=c(BRKS3[-c(3:5,length(BRKS3))],99.5,99.75,99.9,100)
BRKS5=c(-1000,seq(0,30,le=16),100)

# Corresponding colors
COLS0=c("lightblue",rev(heat.colors(length(BRKS0)-7)),"red3","red4","magenta4","magenta","pink")
COLS1=c("lightblue",rev(heat.colors(length(BRKS1)-6)),"red3","red4","magenta4","violet")
COLS2=c(COLS1[-c(3:5,length(COLS1))],"magenta3","magenta","violet","white")

#---------------------------------------------------------------------------------
# A few things for the legends
#---------------------------------------------------------------------------------
# Colors / names for the 3 methods
method.color=function(m){
  if (m=="calend") out=1
  if (m=="annmax") out=4
  if (m=="locmax") out=2
  out
}

method.name=function(m){
  if (m=="calend") out="Calendar"
  if (m=="annmax") out="Annual-maxima"
  if (m=="locmax") out="Local-maxima"
  out
}

# Months names
MONTHS=c("January","February","March","April",
	 "May","June","July","August",
	 "September","October","November","December")

# Functions that convert percentile levels into probabilities (with scientific writing)
perc2prob=function(lev){
  lev=1-lev/100
  out=c(); for (l in lev){
    if (l>=0.099) out=c(out,round(l,1))
    else if (l>=0.0099) out=c(out,round(l,2))
    else if (l>=0.00099) out=c(out,round(l,3))
    else out=c(out,format(l,sci=T))
  }; return(out)
}

prob2leg=function(probs){
  out=c(); for (p in probs){
    dum=round(p,digits=2); if (dum==0) dum=signif(p,digits=1)
    out=c(out,dum)
  }; return(out)
}

# Dates in titles/legends
day2date=function(id){
  date=DATES[iDAYS0[id],]
  return(paste(MONTHS[date$m],date$d))
}

#---------------------------------------------------------------------------------
# Figure to illustrate the methodology
#---------------------------------------------------------------------------------
myfig.meth=function(s=ALL.STAT.1D,id1=1,id2=1,title=""){
  # Indices d,w of the requested time window (id1:id2)
  id0=floor((id1+id2)/2)
  iw0=id2-id1+1
  # Value of the year of interest
  dum0=s[[1]]$value[id0,iw0]
  # Get data, fit, pdf and p1 for all methods
  meth=names(s)
  dum=list(); pdf=list(); p1=list(); for (i in 1:length(s)){
    if (meth[i]=="annmax") dum[[i]]=s[[i]]$data1[1,iw0,]
    else dum[[i]]=s[[i]]$data1[id0,iw0,]
    pdf[[i]]=mypdffun(s[[i]]$distrib,s[[i]]$fit1[id0,iw0,])
    p1[[i]]=myperclevel(dum0,s[[i]]$distrib,s[[i]]$fit1[id0,iw0,])
  }
  # Plot params
  par(mar=c(2.5,0.5,2.5,0.5),mgp=c(1.1,0.1,0),las=1,tcl=0.3)
  xlim=range(unlist(dum)); xlim=xlim+c(-2,1)*0.1*(xlim[2]-xlim[1])
  ylim=c(); for (i in 1:length(pdf)) ylim=c(ylim,pdf[[i]](mean(dum[[i]]))); ylim=c(-1,1.5)*max(ylim)
  leg1=c(); for (m in names(s)) leg1=c(leg1,method.name(m))
  leg2=paste("p1 =",prob2leg(1-unlist(p1)))
  leg=paste(leg1," (",leg2,")",sep="")
  col="darkviolet"; for (m in names(s)) col=c(col,method.color(m))
  # Empty plot
  plot(pdf[[1]],xlim[1],xlim[2],ylim=ylim,axes=F,xlab="",ylab="",type="n")
  # Value dum0
  abline(v=dum0,lwd=2,col=col[1]); text(dum0,ylim[2],Y0,pos=4,font=2,col=col[1])
  # Pdfs
  for (i in 1:length(pdf)) plot(pdf[[i]],xlim[1],xlim[2],lwd=3,col=col[i+1],add=T)
  # Barcodes
  ycb=seq(2,20,le=length(s)+1); ycb=rbind(ycb[1:length(s)],ycb[-1]-1); ycb=ycb*ylim[1]/20
  for (i in 1:length(dum)){ for (d in dum[[i]]) segments(d,ycb[1,i],d,ycb[2,i],lwd=0.5,col=col[i+1]) }
  # Axes, title and legends
  axis(1,pretty(xlim,5))
  titday=day2date(id1); if (id2!=id1) titday=paste(titday,"-",day2date(id2))
  title(main=paste(title,VAR,CITY,"pdf for",titday),adj=0); title(xlab=VARLEG)
  legend("topleft",bty="n",leg=leg,col=col[-1],text.col=col[-1],lwd=3)
#  legend("topleft",bty="n",leg=leg1,col=col[-1],text.col=col[-1],lwd=3)
#  legend("topright",bty="n",leg=leg2,lwd=NA,text.col=col[-1])
}

#---------------------------------------------------------------------------------
# Figure local p1 (1d)
#---------------------------------------------------------------------------------
myfig.1d=function(s=ALL.STAT.1D,lev=BRKS1,col=COLS1,city=CITY,main="p1",title="",unique.colorbar=T){

# Panel plot
npanel=length(s)
if (unique.colorbar)  layout(matrix(1:(npanel+1),1,npanel+1),width=c(rep(0.31,npanel),0.07))
if (!unique.colorbar) layout(matrix(1:(2*npanel),1,2*npanel),width=rep(c(0.31,0.07),npanel))

# Function for each panel plot
dumplot=function(mat,lev,col,ltitle,rtitle,graph.col,colorbar=F){
  # Plot matrix
  par(mar=c(3.5,2.5,1.5,1),las=1,tcl=0.3)
  image(1:nrow(mat),1:ncol(mat),mat,axes=F,main="",xlab="",ylab="",breaks=lev,col=col)
  for (i in 1:4) axis(i,at=c(-999,999),lab=F)
  # X axis : calendar days
  d0=DATES[iDAYS0,]; m0=unique(d0$m)
  mv=which(d0$d==1); mv[mv==1]=NA; abline(v=mv-0.5,lty=2)
  ml=which(d0$d %in% c(1,10,20)) ; axis(1,at=ml,lab=d0$d[ml],mgp=c(0,0.1,0))
  mt=which(d0$d %in% c(5,15,25,30)) ; axis(1,at=mt,lab=F)
  mm=c(); for (m in m0) mm=c(mm,mean(which(d0$m==m))); axis(1,at=mm,lab=MONTHS[m0],tcl=0,mgp=c(0,0.9,0))
  # Y axis : duration of time windows
  axis(2,at=1:nWIN,lab=WINDOWS,mgp=c(0,0.2,0))
  # Titles
  title(ltitle,adj=0); title(rtitle,adj=1,col.main=graph.col)
  title(xlab="Central day",mgp=c(2,0,0))
  title(ylab="Duration in days",mgp=c(1.5,0,0))
  # Adding the global max
  ijmax=which(mat==max(mat,na.rm=T),arr.ind=T)
  points(ijmax,pch="X",font=2,cex=1.5)
  # Adding max for each duration
  imax=as.numeric(apply(mat,2,which.max))
  points(imax,1:nWIN,pch=16,cex=0.6)
  # Colorbar if requested
  if (colorbar){
    par(mar=c(3.5,0,1.5,1))
    levleg=lev[2:(length(lev)-1)]; if (main %in% c("p1","p0")) levleg=perc2prob(levleg)
    mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
  }
}

# Plotting the n panels
llev=lev; if (!typeof(lev)=="list") { llev=list(); for (m in names(s)) llev[[paste(m)]]=lev }
lcol=col; if (!typeof(col)=="list") { lcol=list(); for (m in names(s)) lcol[[paste(m)]]=col }
ltit=c(paste(main,city,EVE),rep("",npanel-1)); ltit=paste(title,ltit)
rtit=c(); for (m in names(s)) rtit=c(rtit,method.name(m)); #rtit=paste(title,rtit)
cols=c(); for (m in names(s)) cols=c(cols,method.color(m))
istat=which(c("p1","p0","FAR")==main)
if (unique.colorbar) cb=c(rep(FALSE,length(s)-1),TRUE) else cb=rep(TRUE,length(s))
for (i in 1:length(s))
  dumplot(100*as.matrix(s[[i]][[istat]]),llev[[i]],lcol[[i]],ltit[i],rtit[i],cols[i],cb[i])

}

#---------------------------------------------------------------------------------
# Figure p1max all domains (2d)
#---------------------------------------------------------------------------------
myfig.2d=function(S=ALL.STAT.2D,lev=BRKS2,col=COLS2,main="p1",title="",unique.colorbar=T){

# Panel plot
npanel=length(S)
if (unique.colorbar)  layout(matrix(1:(npanel+1),1,npanel+1),width=c(rep(0.31,npanel),0.07))
if (!unique.colorbar) layout(matrix(1:(2*npanel),1,2*npanel),width=rep(c(0.31,0.07),npanel))

# Function for each panel plot
dumplot=function(mat,lev,col,ltitle,rtitle,graph.col,colorbar=FALSE){
  # Plot matrix
  par(mar=c(3.5,2.5,1.5,1),las=2,tcl=0.3)
  image(1:nrow(mat),1:ncol(mat),mat,axes=F,main="",xlab="",ylab="",breaks=lev,col=col)
  for (i in 1:4) axis(i,at=c(-999,999),lab=F)
  # Axes
  axis(1,at=1:(nSIZ+1),lab=c("0x0",names(DOMAINS)),mgp=c(0,0.1,0))
  axis(2,at=1:length(WINDOWS),lab=WINDOWS,mgp=c(0,0.2,0))
  # Titles
  title(ltitle,adj=0); title(rtitle,adj=1,col.main=graph.col)
  title(xlab=paste("Domain size in",SPATIAL),mgp=c(2.5,0,0))
  title(ylab="Duration in days",mgp=c(1.5,0,0))
  # Adding the global max
  ijmax=which(mat==max(mat,na.rm=T),arr.ind=T)
  points(ijmax,pch="X",font=2,cex=1.5)
  # Colorbar if requested
  if (colorbar){
    par(mar=c(3.5,0,1.5,1))
    levleg=lev[2:(length(lev)-1)]; if (main %in% c("p1","p0")) levleg=perc2prob(levleg)
    mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
  }
}

# Plotting the n panels with a unique colorbar or not
llev=lev; if (!typeof(lev)=="list") { llev=list(); for (m in names(S)) llev[[paste(m)]]=lev }
lcol=col; if (!typeof(col)=="list") { lcol=list(); for (m in names(S)) lcol[[paste(m)]]=col }
ltit=c(paste(main,EVE),rep("",npanel-1)); ltit=paste(title,ltit)
rtit=c(); for (m in names(S)) rtit=c(rtit,method.name(m)); #rtit=paste(title,rtit)
cols=c(); for (m in names(S)) cols=c(cols,method.color(m))
istat=which(c("p1","p0","FAR")==main)
if (unique.colorbar) cb=c(rep(FALSE,length(S)-1),TRUE) else cb=rep(TRUE,length(S))
for (i in 1:length(S))
  dumplot(100*as.matrix(S[[i]][,,istat]),llev[[i]],lcol[[i]],ltit[i],rtit[i],cols[i],cb[i])

}

