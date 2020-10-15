library(geosphere)

#--------------------------------------------
# Data grid
#--------------------------------------------
PTS=data.frame(lon=rep(XRAW$lon,length(XRAW$lat)),lat=rep(XRAW$lat,each=length(XRAW$lon)))

#--------------------------------------------
# Processing of countries belonging to PTS
#--------------------------------------------
# Original countries
cn=map.where("world",PTS)

# Remove sub-regions (e.g. France:Corsica > France)
for (i in 1:length(cn)) cn[i]=strsplit(cn[i],":")[[1]][1]

# Print list of countries belonging to the domain
#sort(unique(na.omit(cn)))

# Selection of countries for the analysis
retained=c("Albania","Andorra","Austria","Belarus","Belgium","Bosnia and Herzegovina","Bulgaria",
	   "Croatia","Czech Republic","Denmark","Estonia","Finland","France","Germany","Greece",
	   "Hungary","Ireland","Italy","Kosovo","Latvia","Lithuania","Luxembourg","Macedonia",
	   "Moldova","Montenegro","Netherlands","Norway","Poland","Portugal","Romania","Serbia",
	   "Slovakia","Slovenia","Spain","Sweden","Switzerland","UK","Ukraine")
cn[which(! cn %in% retained)]=NA
##
cn.names=sort(unique(na.omit(cn)))

# Print size of countries (in gridpoints... only to have a rough idea)
dum=c(); for (c in cn.names) dum=c(dum,length(which(cn==c)))
#cbind(cn.names,dum)[order(dum),]

# Grouping of (small) countries
cn[which(cn=="Andorra")]="Spain"
#cn[which(cn %in% c("Belgium","Luxembourg","Netherlands"))]="Benelux"
cn[which(cn %in% c("Albania","Bosnia and Herzegovina","Croatia","Kosovo","Macedonia","Montenegro","Serbia","Slovenia"))]="Ex-Yugoslavia"
cn[which(cn=="Moldova")]="Ukraine"
##
cn.names=sort(unique(na.omit(cn)))

# Country numbers
cnn=rep(NA,length(cn))
for (i in 1:length(cnn)) { if (!is.na(cn[i])) cnn[i]=which(cn.names==cn[i]) }

# List with countries information
CN=list(names=cn.names,map=matrix(cn,nrow=length(XRAW$lon)),map2=matrix(cnn,nrow=length(XRAW$lon)))
nCN=length(CN$names)

#--------------------------------------------
# Make domains from countries
#--------------------------------------------

# Vector of country numbers (sizes of domains)
SIZES=c(1,2,3,4,5,6,8,10,12,15,20,nCN)
SIZES=c(1:12,15,18,21,24,nCN)
nSIZ=length(SIZES)

# Central coordinates of countries (no latitude weighting...)
cn.centers=c(); for (c in CN$names)
  cn.centers=rbind(cn.centers,apply(PTS[which(c(CN$map)==c),],2,mean))
#
#mymap(XRAW$lon,XRAW$lat,CN$map2); points(cn.centers)

# Nearest neighbours
cn.nn=c(); for (i in 1:nCN)
  cn.nn=rbind(cn.nn,CN$names[order(distm(cn.centers[i,],cn.centers))]) 

# List of possible domains (with duplicates removed)
DOMAINS=list(); for (i in SIZES){
  domi=c(); for (j in 1:nCN) domi=rbind(domi,sort(cn.nn[j,1:i]))
print(paste(i,":",length(which(duplicated(domi)))))
  if (length(which(duplicated(domi)))>0) domi=domi[-which(duplicated(domi)),]
  if (length(domi)==i) domi=matrix(domi,ncol=i)
  DOMAINS[[paste("c",i,sep="")]]=domi
}

# Plot random domain
plot.random.domain=function(nc=1){
  mymap(XRAW$lon,XRAW$lat,CN$map2)
  dum=DOMAINS[[paste("c",nc,sep="")]]
  points(PTS[which(c(CN$map) %in% dum[sample(1:nrow(dum),1),]),])
}

#pdf("../figures/random_domain.pdf")
#plot.random.domain(5)
#dev.off()
