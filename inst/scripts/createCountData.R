path1='/media/jay/Hitachi2GB/00NMML/ActiveRPack/Harb_Seal_Surv/'
path2='JHop/inst/raw data/'
fl = 'JHI_Seal_GroundCounts_Combined_92-02and07-08_jnw.csv'
d1 = read.csv(paste0(path1,path2,fl))

cutDate = strptime('2000/7/15','%Y/%m/%d')$yday
yday = as.POSIXlt(as.POSIXct(d1$Date))$yday
pup = d1[,'Pups_Out']
pup[yday > cutDate] = NA

dgrnd = data.frame(
  date = as.POSIXct(d1$Date),
  year = as.POSIXlt(as.POSIXct(d1$Date))$year + 1900,
  yday =  yday,
  hour = strptime(d1$Time, format = '%H:%M')$hour,
  minu = strptime(d1$Time, format = '%H:%M')$min,
  adu = d1[,'Non.Pups_Out'],
  pup = pup,
  sky = d1$Sky,
  prec = d1$Precipitation
)

save(dgrnd,file=paste0(path1,'JHop/data/dgrnd.rda'))

path1='/media/jay/Hitachi2GB/00NMML/ActiveRPack/Harb_Seal_Surv/'
path2='JHop/inst/raw data/'
fl = 'dair.csv'
GlacEsts = read.csv(paste0(path1,path2,fl))

yday = strptime(paste0(GlacEsts$year,'/',substr(GlacEsts$moday,1,1),
  '/',substr(GlacEsts$moday,2,3)), '%Y/%m/%d')$yday
dair = data.frame(GlacEsts, yday = yday)

save(dair,file=paste0(path1,'JHop/data/dair.rda'))
