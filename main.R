#---------------------------------------------------------------------------------
# Load rbase.R (packages and base functions) + a few more packages
#---------------------------------------------------------------------------------
source("rbase.R")
library(extRemes)   # GEV
library(evd)        # GEV with fixed params
library(geosphere)  # geometry on Earth (e.g. distance)

#---------------------------------------------------------------------------------
# Namelist of the event
#---------------------------------------------------------------------------------
# Name of the event
EVE="EHW03"

# Variable of interest
VAR="T"
VARLEG="Temperature (K)"

# Period to consider for searching the event
DAY1=20030501     # start date
DAY2=20030930     # end date

# Full period
PERIOD=1950:2018

# Spatial analysis: by lon-lat domains or groups of countries
SPATIAL="lon-lat"
#SPATIAL="country"

#---------------------------------------------------------------------------------
# Get data
#---------------------------------------------------------------------------------
# Get raw data
if (SPATIAL=="lon-lat") XRAW=myno("../data/EOBS_tas_day_raw_Europe_2.5d.nc")
if (SPATIAL=="country") XRAW=myno("../data/EOBS_tas_day_raw_Europe_0.5d.nc")
dts=time2mdy(XRAW$time,XRAW$timeu)
XRAW$var=XRAW$var[,,which(dts$y %in% PERIOD)]
XRAW$time=XRAW$time[which(dts$y %in% PERIOD)]

# Get long-term T trend (HIST+RCP CMIP5 smooth multi-model mean)
if (SPATIAL=="lon-lat") TREND=myno("../data/ENS_tasld_historical_rcp85_JJA_raw_smooth_Europe_lsm.nc")
if (SPATIAL=="country") TREND=myno("../data/ENS_tas_historical_rcp85_JJA_raw_smooth_Europe_0.5d.nc")
dts=time2mdy(TREND$time,TREND$timeu)
TREND$var=TREND$var[,,which(dts$y %in% PERIOD)]
TREND$time=TREND$time[which(dts$y %in% PERIOD)]

#---------------------------------------------------------------------------------
# A few helpful assigns
#---------------------------------------------------------------------------------
nYRS=length(PERIOD)                         # number of years
DATES=time2mdy(XRAW$time,XRAW$timeu)        # dates as data.frame $m $d $y
dates=mdy2yyyymmdd(DATES)                   # dates as number yyyymmdd2mdy
iDAYS0=which(dates>=DAY1 & dates<=DAY2)     # indices of the event period within all dates
nD0=length(iDAYS0)                          # number of days of the event period
Y0=DATES$y[iDAYS0][1]                       # year of the event (to be used for detrending)
iY0=which(PERIOD==Y0)                       # index of the year of the event within all years

#---------------------------------------------------------------------------------
# Load subroutines for functions and figures
#---------------------------------------------------------------------------------
source("functions.R")
source("figures.R")

#---------------------------------------------------------------------------------
# Time windows and spatial domains to consider for searching the event
#---------------------------------------------------------------------------------
# Time windows (vector of durations in days)
WINDOWS=c(1:15,20,25,30,35,40,45,60,75,90,120,150,180)
#WINDOWS=c(1:30,seq(35,90,5),seq(100,150,10))
nWIN=length(WINDOWS)

# Spatial domains (list of sizes with several possible domains per size)
# >> Computed in the lon-lat or country subroutine which provides
#    - DOMAINS (list of domains),
#    - SIZES   (matrix of sizes),
#    - nSIZ    (nb of sizes)
source(paste(SPATIAL,"_analysis.R",sep=""))

#---------------------------------------------------------------------------------
# Time-only analysis at one particular location (1D)
#---------------------------------------------------------------------------------
# Location of interest
LOC=c(2.5,47.5)       ; CITY="Paris"   # one grid point
#LOC=c(-5,10,40,52.5) ; CITY="France"  # one lon-lat domain
#LOC="France"         ; CITY=dom       # one country (only works with country_analysis.R)

# Compute 1d stats for the 3 methods
# See details in compute.stat.1d(), subroutine functions.R
ALL.STAT.1D=list()
ALL.STAT.1D$calend=compute.stat.1d(XRAW,LOC,annual.max=F,calendar.flex=0,dis="gauss")
ALL.STAT.1D$annmax=compute.stat.1d(XRAW,LOC,annual.max=T,calendar.flex=0,dis="gev",fixed.ksi=0)
ALL.STAT.1D$locmax=compute.stat.1d(XRAW,LOC,annual.max=F,calendar.flex=7,dis="gauss")

# Date vs. duration plot of p1 for the 3 methods
myfig.1d(ALL.STAT.1D,lev=BRKS0,col=COLS0,city=CITY,main="p1",title=c("a)","b)","c)"))

# Associated FAR
myfig.1d(ALL.STAT.1D,lev=BRKS4,col=COLS0,city=CITY,main="FAR",title=c("d)","e)","f)"))

# Print optimal event definition for each method (by default: optimize p1)
for (meth in names(ALL.STAT.1D)){
  print(paste("*** Event definition for method",meth,":"))
  print.max.eve.1d(ALL.STAT.1D[[paste(meth)]])
  print("***")
}

# Illustration of the methodology for the event definition of calendar method
dum=get.statmax.1d(ALL.STAT.1D$calend)
idays=dum[which.max(dum[,1]),4:5]
myfig.meth(ALL.STAT.1D,idays[1],idays[2])

#---------------------------------------------------------------------------------
# Full space-time analysis (2D)
#---------------------------------------------------------------------------------
# Compute 2d stats for the 3 methods
# See details in compute.stat.2d(), subroutine functions.R
# (Warning: depending on the number of domains, it can take some time...)
ALL.STAT.2D=list()
ALL.STAT.2D$calend=compute.statmax.2d(XRAW,ALL.STAT.1D$calend,ann=F,cal=0,dis="gauss")
ALL.STAT.2D$annmax=compute.statmax.2d(XRAW,ALL.STAT.1D$annmax,ann=T,cal=0,dis="gev",fixed.ksi=0)
ALL.STAT.2D$locmax=compute.statmax.2d(XRAW,ALL.STAT.1D$locmax,ann=F,cal=7,dis="gauss")

# Domain size vs. time duration plot of minimum p1 for the 3 methods
myfig.2d(ALL.STAT.2D,lev=BRKS0,col=COLS0,main="p1",title=c("a)","b)","c)"))

# Associated FAR
myfig.2d(ALL.STAT.2D,lev=BRKS4,col=COLS0,main="FAR",title=c("d)","e)","f)"))

# Print optimal event definition for each method (by default: optimize p1)
for (meth in names(ALL.STAT.2D)){
  print(paste("*** Event definition for method",meth,":"))
  imax=which(ALL.STAT.2D[[paste(meth)]][,,1]==max(ALL.STAT.2D[[paste(meth)]][,,1]),arr.ind=T)
  print.max.eve.2d(ALL.STAT.2D[[paste(meth)]],imax[1],imax[2])
  print("***")
}

