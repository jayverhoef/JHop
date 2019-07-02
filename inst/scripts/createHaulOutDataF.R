# install haul-out data from Josh's github site
library(devtools)
devtools::install_github('jmlondon/akpvhaulout')
library(akpvhaulout)

data(tbl_percent_glacial)
dt = tbl_percent_glacial
dt = as.data.frame(dt[,1:12])


dt = dt[!is.na(dt$date_time),]
dt$yr = as.POSIXlt(dt[,"date_time"])$year + 1900
dt$dy = as.POSIXlt(dt[,"date_time"])$yday
dt$hr = (as.POSIXlt(dt[,"date_time"])$hour + 16)%%24

dfal = dt[dt$dy > as.POSIXlt("2004-07-15")$yday & 
  dt$dy < as.POSIXlt("2004-09-30")$yday,]
dspr = dt[dt$dy <= as.POSIXlt("2004-07-15")$yday & 
  dt$dy > as.POSIXlt("2004-05-01")$yday,]

# Do spring first
spenoList = unique(dspr$speno)
#order by datadatetime within speno, and remove duplicate times
  dupTF = NULL
  d2 = NULL
  for(i in 1:length(spenoList)) {
    tmp = dspr[dspr$speno == spenoList[i],]
    ntmp = length(tmp[,1])
    oindx = order(tmp$date_time)
    tmp = tmp[oindx,]
    dupTF = c(dupTF,
      any(duplicated(tmp$date_time)) )
    tmp = tmp[!duplicated(tmp$date_time),] 
    d2 = rbind(d2,tmp)
  }
dspr = d2
dspr$speno = as.factor(as.character(dspr$speno))

dspr$yrhr = (dspr$yr - min(dspr$yr))*365*24 + dspr$dy*24 + as.POSIXlt(dspr[,"date_time"])$hour
dspr$yrhr0 = dspr$yrhr
for(i in 1:length(levels(dspr$speno)))
	dspr[dspr$speno == levels(dspr$speno)[i],'yrhr0'] = 
		dspr[dspr$speno == levels(dspr$speno)[i],'yrhr'] -
    min(dspr[dspr$speno == levels(dspr$speno)[i],'yrhr'])


dspr = dspr[order(dspr$speno,dspr$yrhr0),]
dspr$dry = dspr$percent_dry/100
dspr$y = (dspr$dry*(length(dspr$dry) - 1) + 0.5)/
  length(dspr$dry)
dspr$dystd = (dspr$dy-as.POSIXlt("2004-06-15")$yday)/30
dspr$speno = as.factor(as.character(dspr$speno))

# Do fall next
spenoList = unique(dfal$speno)
#order by datadatetime within speno, and remove duplicate times
  dupTF = NULL
  d2 = NULL
  for(i in 1:length(spenoList)) {
    tmp = dfal[dfal$speno == spenoList[i],]
    ntmp = length(tmp[,1])
    oindx = order(tmp$date_time)
    tmp = tmp[oindx,]
    dupTF = c(dupTF,
      any(duplicated(tmp$date_time)) )
    tmp = tmp[!duplicated(tmp$date_time),] 
    d2 = rbind(d2,tmp)
  }
dfal = d2
dfal$speno = as.factor(as.character(dfal$speno))

dfal$yrhr = (dfal$yr - min(dfal$yr))*365*24 + dfal$dy*24 + as.POSIXlt(dfal[,"date_time"])$hour
dfal$yrhr0 = dfal$yrhr
for(i in 1:length(levels(dfal$speno)))
	dfal[dfal$speno == levels(dfal$speno)[i],'yrhr0'] = 
		dfal[dfal$speno == levels(dfal$speno)[i],'yrhr'] -
    min(dfal[dfal$speno == levels(dfal$speno)[i],'yrhr'])


dfal = dfal[order(dfal$speno,dfal$yrhr0),]
dfal$dry = dfal$percent_dry/100
dfal$y = (dfal$dry*(length(dfal$dry) - 1) + 0.5)/
  length(dfal$dry)
dfal$dystd = (dfal$dy-as.POSIXlt("2004-08-15")$yday)/30
dfal$speno = as.factor(as.character(dfal$speno))


path1='/media/jay/Hitachi2GB/00NMML/ActiveRPack/Harb_Seal_Surv/'
dHOspr = dspr
save(dHOspr,file=paste0(path1,'JHop/data/dHOspr.rda'))
dHOfal = dfal
save(dHOfal,file=paste0(path1,'JHop/data/dHOfal.rda'))





