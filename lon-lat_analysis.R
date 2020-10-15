###############################################################################################
# SUBROUTINES
###############################################################################################

#---------------------------------------------------------------------------------
# Funnction that computes all possible domains of various sizes within lonlim/latlim.
#   > sizes (matrix 2 x nsiz, sizes in nb of grid points, with row1>=row2 [see below]),
#     lonlim/latlim (vectors 2 values), dlon/dlat (resolution in deg).
#   < list with names = sizes in nb of grid points (1x1, 2x1, 2x2, etc..)
#---------------------------------------------------------------------------------
make.all.possible.domains=function(sizes,lonlim,latlim,dlon,dlat){
  # Maximum size of possible domains
  sizlonmax=1+(lonlim[2]-lonlim[1])/dlon ; sizlatmax=1+(latlim[2]-latlim[1])/dlat
  # Subroutine that generates domains of a given size
  dumfun=function(s){
    if (s[1]>sizlonmax | s[2]>sizlatmax) return(NULL)
    else {
      x1=seq(lonlim[1],lonlim[2]-(s[1]-1)*dlon,dlon)
      y1=seq(latlim[1],latlim[2]-(s[2]-1)*dlat,dlat)
      x2=x1+(s[1]-1)*dlon; y2=y1+(s[2]-1)*dlat; nx=length(x1); ny=length(y1)
      return(cbind(rep(x1,each=ny),rep(x2,each=ny),rep(y1,nx),rep(y2,nx)))
    }
  }
  # Loop on sizes: for size si, generates rectangular domains si[1]xsi[2] AND si[2]xsi[1]
  doms=list(); for (i in 1:ncol(sizes)){ si=sizes[,i]
    doms.i=dumfun(si)
    if (si[1]!=si[2]) doms.i=rbind(doms.i,dumfun(si[2:1]))
    doms[[paste(si,collapse="x")]]=doms.i
  }
  # Ouput
  doms
}
  
###############################################################################################
# MAIN
###############################################################################################

# Limits of the region considered for searching the event
LONLIM=c(-10,40) 
LATLIM=c(35,70)

# Area that must be encompassed within possible domains (at least 1 grid point)
AREA0=c(2.5,2.5,47.5,47.5)

# Resolution of gridded data
GRIDRES=XRAW$lon[2]-XRAW$lon[1]
  
# Domain sizes considered for searching the event, in number of grid points
# (matrix with 2 rows and row1>=row2 [the function then generates portrait or landscape])
SIZES=rbind(c(1,rep(c(2,seq(3,15,2)),each=2),17,19,21),c(rep(c(1:2,seq(3,15,2)),each=2),15,15))
SIZES=rbind(c(1,3,5,7,11,15,21),c(1,3,5,7,11,15,15))
nSIZ=ncol(SIZES)
print(paste("Initial number of sizes:",nSIZ))

# Generating all possible domains for each SIZE
DOMAINS=make.all.possible.domains(SIZES,LONLIM,LATLIM,GRIDRES,GRIDRES)
nSIZ=length(DOMAINS)
print(paste("Initial number of domains:",sum(sapply(DOMAINS,nrow)),"(",nSIZ,"sizes)"))

# Removing domains that do not contain AREA0
killdom=function(x) (x[1]>AREA0[2] | x[2]<AREA0[1] | x[3]>AREA0[4] | x[4]<AREA0[3])
for (i in 1:nSIZ){ ai=DOMAINS[[i]]; bi=c(); for (j in 1:nrow(ai)){ aij=ai[j,]
  if (!killdom(aij)) bi=rbind(bi,aij) } ; DOMAINS[[i]]=unname(bi)
}
nSIZ=length(DOMAINS)
print(paste("After area filter:",sum(sapply(DOMAINS,nrow)),"(",nSIZ,"sizes)"))

# Removing intermediate domains (2.5°) for sizes >= 3 grid points
killdom=function(x) { tx=trunc(x); ((tx[1]!=x[1] & tx[2]!=x[2]) | (tx[3]!=x[3] & tx[4]!=x[4])) }
i0=min(which(SIZES[1,]<=3 & SIZES[2,]>=3))
for (i in 4:nSIZ){ ai=DOMAINS[[i]]; bi=c(); for (j in 1:nrow(ai)){ aij=ai[j,]
  if (!killdom(aij)) bi=rbind(bi,aij) } ; DOMAINS[[i]]=unname(bi)
}
nSIZ=length(DOMAINS)
print(paste("After 5deg filter:",sum(sapply(DOMAINS,nrow)),"(",nSIZ,"sizes)"))

# Removing domains that do not have data (e.g. sea)
killdom=function(x) is.na(myfldmean(TREND,x)$var[1])
for (i in 1:nSIZ){ ai=DOMAINS[[i]]; bi=c(); for (j in 1:nrow(ai)){ aij=ai[j,]
  if (!killdom(aij)) bi=rbind(bi,aij) } ; DOMAINS[[i]]=unname(bi)
}
nSIZ=length(DOMAINS)
print(paste("After data filter:",sum(sapply(DOMAINS,nrow)),"(",nSIZ,"sizes)"))

# Restricting SIZES to those with generated domains
SIZES=matrix(unlist(strsplit(names(DOMAINS),"x")),nrow=2)

