#---------------------------------------------------------------------------------
# A few useful functions
#---------------------------------------------------------------------------------
# Running mean accounting for NAs and centering
myrunavg=function(x,n) filter(x,rep(1/n,n))

# Maximum accouting for NAs
mymax=function(x) {out=NA; if (any(!is.na(x))) out=max(x,na.rm=T); out}

# Various fits of parametric laws, output is systematically 3 parameters (possibly NA)
myfitparams=function(dum,distrib,ksi=NULL){
#dum[iY0]=NA ### to exclude the value of the year of interest (we don't want to do that)
  dum=na.omit(dum); out=rep(NA,3); if (length(dum)>0){
  if (distrib=="gauss")  out=c(mean(dum),sd(dum),NA)
  if (distrib=="gev")    out=unname(fevd(dum,method="Lmoments")$results)
  if (distrib=="gevfix") out=unname(fgev(dum,shape=ksi)$param)
  if (distrib=="gumbel") out=c(unname(fevd(dum,type="Gumbel")$results$par),0)
  if (distrib=="snorm")  out=unname(snormFit(dum)$par)
  }
  out
}

# PDF function for various parametric laws
mypdffun=function(distrib,params){
  if (distrib=="gauss")  out=function(x) dnorm(x,params[1],params[2])
  if (distrib %in% c("gev","gevfix","gumbel")) {
    evfit=params; out=function(x) devd(x,evfit[1],evfit[2],evfit[3])}
  if (distrib=="snorm")  out=function(x) dsnorm(x,params[1],params[2],params[3])
  out
}

# Percentile level for various parametric laws
myperclevel=function(dum0,distrib,params){
  out=NA; if (!is.na(dum0) & any(!is.na(params))){
  if (distrib=="gauss")  out=pnorm(dum0,params[1],params[2])
  if (distrib %in% c("gev","gevfix","gumbel")) out=pevd(dum0,params[1],params[2],params[3])
  if (distrib=="snorm")  out=psnorm(dum0,params[1],params[2],params[3])
  }
  out
}

#---------------------------------------------------------------------------------
# Function that field-averages a myno() object over a spatial domain.
#   > data (list) and domain (vector c(lon0,lat0), c(lon1,lon2,lat1,lat2), or c("country1",...))
#   < list with var, time and timeu
# Warning:
#   - data$lon must be in ascending order (e.g. without break 360/0 or 180/-180)
#   - the "country" option requires the CN object created by country_analysis.R
#---------------------------------------------------------------------------------
myfldmean=function(data,domain,coslat=TRUE){
  if (length(dim(data$var))==2) data$var=array(data$var,dim=c(dim(data$var),1))
  if (class(domain)=="numeric"){
    if (length(domain)==2){
      ilon=which(data$lon==domain[1]); ilat=which(data$lat==domain[2])
      out=data$var[ilon,ilat,]
    }
    if (length(domain)==4){
      ilon=which(data$lon>=domain[1] & data$lon<=domain[2])
      ilat=which(data$lat>=domain[3] & data$lat<=domain[4])
      mat=matrix(c(data$var[ilon,ilat,]),ncol=length(data$time))
      wghts=1 ; if (coslat) wghts=rep(cos(data$lat[ilat]*pi/180),each=length(ilon))
      out=apply(mat,2,mywmean,wghts)
    }
  }
  if (class(domain)=="character"){
    dmat=matrix(c(data$var),ncol=length(data$time))
    wghts=rep(1,nrow(dmat)) ; if (coslat) wghts=rep(cos(data$lat*pi/180),each=length(data$lon))
    ii=which(c(CN$map)%in%domain)
    out=apply(dmat[ii,],2,mywmean,wghts[ii])
  }
  list(var=out,time=data$time,timeu=data$timeu)
}

#---------------------------------------------------------------------------------
# Function that returns a list of stats (p1, p0, far) for a given location (station,
# lon-lat domain, country) and for all time windows requested in main.R (WINDOWS).
#   > ts (vector or myno() object), domain (also used for detrending), annual.maxima
#     for fitting annual max rather than calendar days, calendar.flex for possible
#     calendar local window, distrib (gauss, gev, etc.), and fixed.ksi (for EVT distribs)
#   < stats in a data.frame object with colnames = WINDOWS
# Warning: the version provided below only includes a very simple detrending (the 
#          correction for climate change does not depend on the time window), which
#          can only be used with temperature variables (see C&R2018 for details)
#---------------------------------------------------------------------------------
compute.stat.1d=function(ts,domain,annual.maxima=FALSE,calendar.flex=0,
			 distrib="gauss",fixed.ksi=NULL){
  #---- If ts is a myno() object, it gets averaged over domain
  if (typeof(ts)=="list") ts=myfldmean(ts,domain)$var
  #---- Building arrays dum (samples) and dumY0 (values of the year of interest)
  dum=array(NA,dim=c(nD0,nWIN,nYRS)); dumY0=matrix(NA,nD0,nWIN)
  # Loop on time windows
  for (iw in 1:nWIN){ w=WINDOWS[iw]
    # Running mean of ts over w
    tsw=myrunavg(ts,w)
    # dumY0 = value year of interest
    for (id in 1:nD0){ d=iDAYS0[id] ; dumY0[id,iw]=tsw[d] }
    # dum = sample of the nYRS values :
    # - if annual.maxima, same sample for all calendar days (only id=1 is filled)
    if (annual.maxima){
      dumw=c(); for (y in PERIOD) dumw=c(dumw,mymax(tsw[which(DATES$y==y)]))
      dum[1,iw,]=dumw
    }
    # - else, loop on calendar days
    if (!annual.maxima){ for (id in 1:nD0){ d=iDAYS0[id]
      idays=which(DATES$m==DATES$m[d] & DATES$d==DATES$d[d])
      if (calendar.flex==0) dum[id,iw,]=tsw[idays]
      else {
        flex=-calendar.flex:calendar.flex
        idays=rep(idays,each=length(flex))+rep(flex,nYRS)
        idays[idays<1]=1 ; idays[idays>nrow(DATES)]=nrow(DATES)
        dum[id,iw,]=apply(matrix(tsw[idays],nrow=length(flex)),2,mymax)
      }
    }}
  }
  #---- Computing data detrended wrt. t=Y0 (dum1, for p1) and t=1 (dum0, for p0)
  tr=myfldmean(TREND,domain)$var
  mydetrend=function(ts,year) ts-(tr-tr[which(PERIOD==year)])
  dum1=dum; dum0=dum; for (iw in 1:nWIN){ for (id in 1:nD0){ if (any(!is.na(dum[id,iw,]))){
      dum1[id,iw,]=mydetrend(dum[id,iw,],Y0)
      dum0[id,iw,]=mydetrend(dum[id,iw,],PERIOD[1])
  }}}
  #---- Fitting parameters of the requested distribution
  if (distrib=="gev" & !is.null(fixed.ksi)){
    fit1=dum1[,,1:3]; fit0=dum0[,,1:3]
    ksi=rep(fixed.ksi,length(WINDOWS))[1:length(WINDOWS)]
    for (iw in 1:nWIN){ for (id in 1:nD0){ if (any(!is.na(dum[id,iw,]))){
      fit1[id,iw,]=myfitparams(dum1[id,iw,],"gevfix",ksi=ksi[iw])
      fit0[id,iw,]=myfitparams(dum0[id,iw,],"gevfix",ksi=ksi[iw])
    }}}
  }
  else {
    fit1=aperm(apply(dum1,1:2,myfitparams,distrib),c(2,3,1))
    fit0=aperm(apply(dum0,1:2,myfitparams,distrib),c(2,3,1))
  }
  #---- If annual.maxima, fits are repeated over all days
  if (annual.maxima){ for (id in 2:nD0){ fit1[id,,]=fit1[1,,] ; fit0[id,,]=fit0[1,,] }}  
  #---- Computing p1 (percentile level of dumY0 within fit1)
  out.p1=abind(dumY0,fit1,along=3)
  out.p1=apply(out.p1,1:2,function(x) myperclevel(x[1],distrib,x[-1]))
  out.p1=as.data.frame(out.p1); names(out.p1)=WINDOWS
  #---- Computing p0 (percentile level of dumY0 within fit0)
  out.p0=abind(dumY0,fit0,along=3)
  out.p0=apply(out.p0,1:2,function(x) myperclevel(x[1],distrib,x[-1]))
  out.p0=as.data.frame(out.p0); names(out.p0)=WINDOWS
  #---- Computing FAR with default FAR=0 if p0=p1=1 (warning: p0 and p1 are percentile levels here!)
  out.far=1-(1-out.p0)/(1-out.p1)
  out.far[which(out.p0==1 & out.p1==1,arr.ind=T)]=0
  #---- Output
  list(p1=out.p1,p0=out.p0,far=out.far,
       data=dum,data1=dum1,data0=dum0,value=dumY0,
       fit1=fit1,fit0=fit0,distrib=distrib)
}

#---------------------------------------------------------------------------------
# Functions that gets/prints information on the max of compute.stat.1d()
#   > output of compute.stat.1d() and stat to optimize (by default: p1)
#   < matrix nwin x 6 (p1, p0, far, iday1, iday2, 0 [domain size])
#---------------------------------------------------------------------------------
get.statmax.1d=function(stat.1d,stat.max="p1"){
  out=c(); for (iw in 1:nWIN){
    imax=which.max(stat.1d[[paste(stat.max)]][,iw])
    outw=rep(NA,6); if (length(imax)>0){
      day1=imax-trunc((WINDOWS[iw]-1)/2); day2=day1+WINDOWS[iw]-1
      outw=c(stat.1d$p1[imax,iw],stat.1d$p0[imax,iw],stat.1d$far[imax,iw],day1,day2,0)
    }
    out=rbind(out,outw)
  }
  unname(out)
}

print.max.eve.1d=function(stat.1d,stat.max="p1"){
  dum=get.statmax.1d(stat.1d,stat.max)
  imax=which(c("p1","p0","far")==stat.max)
  dum=unlist(dum[which.max(dum[,imax]),])
  print(paste("Period =",paste(dates[iDAYS0][c(dum[4],dum[5])],collapse="-")))
  print(paste("--- p1 =",dum[1],"(percentile level)"))
  print(paste("--- p0 =",dum[2],"(percentile level)"))
  print(paste("-- far =",dum[3]))
}

#---------------------------------------------------------------------------------
# Function that returns an array of p1 max (by default) for all DOMAINS/WINDOWS
# requested in main.R (including the local domain, see argument stat.1d)
# + the associated p0 and FAR values.
#   > data (from myno()), stat.1d (from compute.stat.1d()), stat.max (p1, p0 or far),
#     + same arguments as compute.stat.1d.
#   < array nSIZ x nWIN x 6 (p1, p0, far, iday1, iday2, domain [index of DOMAINS])
#---------------------------------------------------------------------------------
compute.statmax.2d=function(dat,stat.1d,stat.max="p1",annual.maxima=FALSE,calendar.flex=0,
			    distrib="gauss",fixed.ksi=NULL){
  # Initialization
  out=array(NA,dim=c(nSIZ,nWIN,6))
  # Loop on SIZES
  for (is in 1:nSIZ){
    print(paste("Size:",is,"/",nSIZ))
    # Loop on DOMAINS corresponding to SIZES[is]
    ndom=nrow(DOMAINS[[is]])
    dums=c(); for (id in 1:ndom){ d=DOMAINS[[is]][id,]
      print(paste("Domain:",id,"/",ndom)); print(d)
      # Computing local stats on the domain (with compute.stat.1d)
      statsd=compute.stat.1d(dat,d,annual.maxima,calendar.flex,distrib,fixed.ksi)
      # For each time duration, optimizing the time window
      # and saving associated p1 p0 far iday1 iday2 (with get.statmax.1d)
      dumsd=get.statmax.1d(statsd,stat.max)[,-6]
      # Saving in a suitable array format
      dums=abind(dums,t(dumsd),along=3)
    }
    # For each domain size, optimizing the spatial domain
    # and saving associated p1 p0 far iday1 iday2 idom
    dum=c(); for (iw in 1:nWIN){
      imax=which.max(c(dums[which(c("p1","p0","far")==stat.max),iw,]))
      dum=rbind(dum,c(dums[,iw,imax],imax))
    }
    out[is,,]=dum
  }
  # Concatenating with local stats
  out=abind(get.statmax.1d(stat.1d,stat.max),out,along=1)
  # Output
  out
}

#---------------------------------------------------------------------------------
# Function that prints information on the max of compute.stat.2d()
#   > output of compute.stat.1d() and indexes of domain size et time duration
#---------------------------------------------------------------------------------
print.max.eve.2d=function(stat.2d,idom,iwin){
  dum=stat.2d[idom,iwin,]
  print(paste("Period =",paste(dates[iDAYS0][c(dum[4],dum[5])],collapse="-")))
  print(paste("Domain =",paste(DOMAINS[[idom]][dum[6],],collapse=";")))
  print(paste("--- p1 =",dum[1],"(percentile level)"))
  print(paste("--- p0 =",dum[2],"(percentile level)"))
  print(paste("-- far =",dum[3]))
}


