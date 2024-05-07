
library(lubridate)
library(raster)
library(dplyr)
library(adehabitatHR)
library(adehabitatLT)
library(sp)
library(trip)
library(ggplot2)
library(patchwork)
library(terra)
library(reshape2)


### Load raster mask and generate spatial coordinate system objects
### load raster mask and generate spatial coordinate system objects
mask <- raster("C:/D1MPA trimestral/mask.tif") # mask
mask<-mask*0
plot(mask)

EPSG6932<-crs(mask) # extract the coordinate spatial system from the object

EPSG4326<-CRS("+proj=longlat +datum=WGS84 +no_defs")


D1<-shapefile("C:/D1MPA_Review/D1_6932.shp")




###------Humpback whale----------------

hump<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/humpback_whales/fastOutTracks_humpback.csv")
hump2<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/humpback_whales/tracks.csv")

hump2df<-data.frame(id=as.factor(hump2$individual_id),
                    date=hump2$date,
                    lc=hump2$location_quality,
                    lon=hump2$decimal_longitude,
                    lat=hump2$decimal_latitude)

humps<-rbind(hump,hump2df)

humps$timestamp<-as.POSIXct(strptime(humps$date, 
                                     format="%Y-%m-%d %H:%M:%S", tz="GMT"))

humps <- humps[!duplicated(humps[, c("lon", "lat","timestamp","id")]), ]

humpd<-subset(humps,lon<0 & lon>(-150))

humpd$month<-month(humpd$timestamp)

humpd<-na.omit(humpd)

humpsp<-SpatialPointsDataFrame(humpd[4:5],humpd,proj4string = EPSG4326)

hbwspW<-sp::spTransform(humpsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

hbwspW<-raster::intersect(hbwspW,D1)
#plot(hbwspW)

hbw<-data.frame(hbwspW)

hbw$id<-factor(hbw$id)

hbw<- hbw[!duplicated(hbw[, c("timestamp","id")]), ]

hbw.lt<-as.ltraj(xy=hbw[15:16],date=hbw$timestamp,id=hbw$id,proj4string = EPSG6932)

hbw.ld<-ld(hbw.lt)

ggplot(hbw.ld,aes(dist))+geom_histogram()

hbw.ld2<-subset(hbw.ld,dist<50000) #eliminate positions based on consecutive points distance 

length(unique(hbw.ld2$id))
length((hbw.ld2$id))
summary(as.factor(year(hbw.ld2$date)))
summary(as.factor(month(hbw.ld2$date)))

#write.csv(hbw.ld2,"C:/D1MPA trimestral/Data/Whale_Humpback.csv")

hbw.ld2$quarter<-quarter(hbw.ld2$date)

hbw.ld2$idy<-paste(hbw.ld2$id,year(hbw.ld2$date),month(hbw.ld2$date),week(hbw.ld2$date))

hbw.q1<-subset(hbw.ld2,quarter=="1")
hbw.q2<-subset(hbw.ld2,quarter=="2")


# JFM
d <- data.frame(x=hbw.q1$x,y=hbw.q1$y, tms=hbw.q1$date, id=hbw.q1$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/31.45

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Whale_Humpback_JFM.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Whale_Humpback_JFM.tif", overwrite=TRUE)


# AMJ 

d <- data.frame(x=hbw.q2$x,y=hbw.q2$y, tms=hbw.q2$date, id=hbw.q2$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/35.7

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Whale_Humpback_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Whale_Humpback_AMJ.tif", overwrite=TRUE)


### -------------- minke whale----------------

minke<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/minke_whales/fastOutTracks_minke.csv")

minke$timestamp<-as.POSIXct(strptime(minke$date, 
                                     format="%Y-%m-%d %H:%M:%S", tz="GMT"))

minkes <- minke[!duplicated(minke[, c("lon", "lat","timestamp","id")]), ]

minked<-subset(minkes,lon<0 & lon>(-150))

minked$month<-month(minked$timestamp)

minked<-na.omit(minked)

minkem<-data.frame(minked$month)

minkesp<-SpatialPointsDataFrame(minked[4:5],minked,proj4string = EPSG4326)

minkespW<-sp::spTransform(minkesp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
minkespW<-raster::intersect(minkespW,D1)

amw<-data.frame(minkespW)

amw$id<-factor(amw$id)

amw<- amw[!duplicated(amw[, c("timestamp","id")]), ]

amw.lt<-as.ltraj(xy=amw[15:16],date=amw$timestamp,id=amw$id,proj4string = EPSG6932)

amw.ld<-ld(amw.lt)

ggplot(amw.ld,aes(dist))+geom_histogram()

amw.ld2<-subset(amw.ld,dist<50000) #eliminate positions based on consecutive points distance 

amw.ld2$quarter<-quarter(amw.ld2$date)

amw.ld2$idy<-paste(amw.ld2$id,year(amw.ld2$date),month(amw.ld2$date),week(amw.ld2$date))

#write.csv(amw.ld2,"C:/D1MPA trimestral/Data/Whale_AntarcticMinke.csv")

amw.q1<-subset(amw.ld2,quarter=="1")
amw.q2<-subset(amw.ld2,quarter=="2")


# JFM
d <- data.frame(x=amw.q1$x,y=amw.q1$y, tms=amw.q1$date, id=amw.q1$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/0.49

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Whale_AntMinke_JFM.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Whale_AntMinke_JFM.tif", overwrite=TRUE)


# AMJ 

d <- data.frame(x=amw.q2$x,y=amw.q2$y, tms=amw.q2$date, id=amw.q2$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/0.20

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Whale_AntMinke_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Whale_AntMinke_AMJ.tif", overwrite=TRUE)


### ---------- Crabeater Seal----------------

#RAATD

crab<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/SCAR_EGBAMM_RAATD_2018_Standardised/RAATD_CRAS_standardized.csv")


crab$timestamp<-as.POSIXct(strptime(paste(paste(crab$year,crab$month,crab$day,sep="-"),
                                          crab$time), 
                                    format="%Y-%m-%d %H:%M:%S", tz="GMT"))

crab<-subset(crab,location_quality!="Z")

crabd<-data.frame(Lon=crab$decimal_longitude,
                  Lat=crab$decimal_latitude,
                  id=as.factor(crab$individual_id),Spp=crab$abbreviated_name,
                  timestamp=crab$timestamp)

crabd<-subset(crabd,Lon<0 & Lon>(-150))

crabd$month<-month(crabd$timestamp)

crabd<-na.omit(crabd)

crabm<-data.frame(crabd$month)

crabsp<-SpatialPointsDataFrame(crabd[1:2],crabd,proj4string = EPSG4326)

crabspW<-sp::spTransform(crabsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

crabspW<-raster::intersect(crabspW,D1)

#plot(crabspW)

crab<-data.frame(crabspW)

crab$id<-factor(crab$id)

crab<- crab[!duplicated(crab[, c("timestamp","id")]), ]

crab.lt<-as.ltraj(xy=crab[14:15],date=crab$timestamp,id=crab$id,proj4string = EPSG6932)

crab.ld<-ld(crab.lt)

ggplot(crab.ld,aes(dist))+geom_histogram()

crab.ld2<-subset(crab.ld,dist<50000) #eliminate positions based on consecutive points distance 

crab.ld2$quarter<-quarter(crab.ld2$date)
crab.ld2$idy<-paste(crab.ld2$id,year(crab.ld2$date),month=month(crab.ld2$date))


#write.csv(crab.ld2,"C:/D1MPA trimestral/Data/Seal_Crabeater.csv")

crab.q1<-subset(crab.ld2,quarter=="1")
crab.q2<-subset(crab.ld2,quarter=="2")
crab.q3<-subset(crab.ld2,quarter=="3")
crab.q4<-subset(crab.ld2,quarter=="4")


# amj

d <- data.frame(x=crab.q2$x,y=crab.q2$y, tms=crab.q2$date, id=crab.q2$id)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/71.221

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Crabeater_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Crabeater_AMJ.tif", overwrite=TRUE)


# jas


d <- data.frame(x=crab.q3$x,y=crab.q3$y, tms=crab.q3$date, id=crab.q3$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/77.44

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Crabeater_JAS.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Crabeater_JAS.tif", overwrite=TRUE)

#ond

d <- data.frame(x=crab.q4$x,y=crab.q4$y, tms=crab.q4$date, id=crab.q4$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/31.11

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Crabeater_OND.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Crabeater_OND.tif", overwrite=TRUE)



### ---------Fur Seal---------------
# Hinke et al. 2017

jh17<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/pone.0170132.s001/Data/satellite telemetry.csv")

head(jh17)

summary(as.factor(jh17$Loc.Qual))

jh17<-subset(jh17,Loc.Qual!="Z")
jh17<-subset(jh17,Loc.Qual!="0")
#jh17<-subset(jh17,Loc.Qual!="")

jh17$timestamp<-as.POSIXct(strptime(paste(jh17$Date,jh17$Time), format="%m/%d/%Y %H:%M:%S", tz="GMT"))

jh17<-data.frame(Lon=jh17$Longitude,Lat=jh17$Latitude,id=as.factor(jh17$Deployment),
                 Spp=jh17$Spp,timestamp=jh17$timestamp,Col=jh17$Site)
head(jh17)

summary(as.factor(jh17$Spp))

jhf<-subset(jh17,Spp=="AFS")

head(jhf)
summary(jhf$timestamp)

jhd<-jhf


# RAATD

ra<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/SCAR_EGBAMM_RAATD_2018_Standardised/RAATD_ANFS_standardized.csv")


ra$timestamp<-as.POSIXct(strptime(paste(paste(ra$year,ra$month,ra$day,sep="-"),
                                        ra$time), 
                                  format="%Y-%m-%d %H:%M:%S", tz="GMT"))
head(ra)

summary(as.factor(ra$location_quality))

ra<-subset(ra,location_quality!="Z")


raat<-data.frame(Lon=ra$decimal_longitude,
                 Lat=ra$decimal_latitude,
                 id=as.factor(ra$individual_id),Spp=ra$abbreviated_name,
                 timestamp=ra$timestamp,Col=c("SOI"))

raatd<-subset(raat,Lon<0 & Lon>(-150))

head(raatd)

raatd$Spp<-c("AFS")

###Krause et al. 2022

kr<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/FurSealShirreff/AFSPOW_postapril.csv")

kr$timestamp<-as.POSIXct(strptime(kr$Date, 
                                  format="%m/%d/%Y %H:%M", tz="GMT"))

krd<-data.frame(Lon=kr$Longitude,
                Lat=kr$Latitude,
                id=as.factor(kr$Tag),Spp=kr$Spp,
                timestamp=kr$timestamp,Col=c("CS"))

afs<-rbind(jhd,raatd,krd)

afs<-na.omit(afs)

afs$id<-as.factor(afs$id)

afs$month<-month(afs$timestamp)

afs$pres<-c(1)

summary(as.factor(afs$Col))

afsp<-SpatialPointsDataFrame(afs[1:2],afs,proj4string = EPSG4326)

afspW<-sp::spTransform(afsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

afspW<-raster::intersect(afspW,D1)

#plot(afspW)

afs<-data.frame(afspW)

afs$id<-factor(afs$id)

afs<- afs[!duplicated(afs[, c("timestamp","id")]), ]

afs.lt<-as.ltraj(xy=afs[16:17],date=afs$timestamp,id=afs$id,proj4string = EPSG6932)

afs.ld<-ld(afs.lt)

afs.ld2<-subset(afs.ld,dist<50000) #eliminate positions based on consecutive points distance 

summary(as.factor(year(afs.ld2$date)))

afs.ld2$year<-year(afs.ld2$date)

afs.ld2$IDy<-paste(afs.ld2$id,year(afs.ld2$date),month(afs.ld2$date))

#write.csv(afs.ld2,"C:/D1MPA trimestral/Data/Seal_Fur.csv")

afs.ld2$quarter<-quarter(afs.ld2$date)
afs.q1<-subset(afs.ld2,quarter=="1")
afs.q2<-na.omit(subset(afs.ld2,quarter=="2"))
afs.q3<-subset(afs.ld2,quarter=="3")
afs.q4<-subset(afs.ld2,quarter=="4")

#jfm

d <- data.frame(x=afs.q1$x,y=afs.q1$y, tms=afs.q1$date, id=afs.q1$IDy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/17.7

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Fur_JFM.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Fur_JFM.tif", overwrite=TRUE)

# amj

d <- data.frame(x=afs.q2$x,y=afs.q2$y, tms=afs.q2$date, id=afs.q2$IDy)
tr <- trip(d,correct_all = T)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/20.94

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Fur_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Fur_AMJ.tif", overwrite=TRUE)


# jas
afs.q3$IDy<-paste(afs.q3$id,year(afs.q3$date),week(afs.q3$date))
d <- data.frame(x=afs.q3$x,y=afs.q3$y, tms=afs.q3$date, id=afs.q3$IDy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/7.24

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Fur_JAS.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Fur_JAS.tif", overwrite=TRUE)

#ond

afs.q4$IDy<-paste(afs.q4$id,year(afs.q4$date),week(afs.q4$date))
d <- data.frame(x=afs.q4$x,y=afs.q4$y, tms=afs.q4$date, id=afs.q4$IDy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/4.73

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Fur_OND.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Fur_OND.tif", overwrite=TRUE)


###----------Elephant Seal------------


soes<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/SCAR_EGBAMM_RAATD_2018_Standardised/RAATD_SOES_standardized.csv")


soes$timestamp<-as.POSIXct(strptime(paste(paste(soes$year,soes$month,soes$day,sep="-"),
                                          soes$time), 
                                    format="%Y-%m-%d %H:%M:%S", tz="GMT"))
soes<-subset(soes,location_quality!="Z")

soesd<-data.frame(Lon=soes$decimal_longitude,
                  Lat=soes$decimal_latitude,
                  id=as.factor(soes$individual_id),Spp=soes$abbreviated_name,
                  timestamp=soes$timestamp)

soesd<-subset(soesd,Lon<0 & Lon>(-150))

soesd$month<-month(soesd$timestamp)

soesd<-na.omit(soesd)

soessp<-SpatialPointsDataFrame(soesd[1:2],soesd,proj4string = EPSG4326)

soesspW<-sp::spTransform(soessp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

soesspW<-raster::intersect(soesspW,D1)

soes<-data.frame(soesspW)

soes$id<-factor(soes$id)

soes<- soes[!duplicated(soes[, c("timestamp","id")]), ]

soes.lt<-as.ltraj(xy=soes[15:16],date=soes$timestamp,id=soes$id,proj4string = EPSG6932)

soes.ld<-ld(soes.lt)

#ggplot(soes.ld2,aes(dist))+geom_histogram()

soes.ld2<-subset(soes.ld,dist<30000) #eliminate positions based on consecutive points distance 

soes.ld2$idy<-paste(soes.ld2$id,year(soes.ld2$date),week(soes.ld2$date))

#write.csv(soes.ld2,"C:/D1MPA trimestral/Data/Seal_Elephant.csv")

soes.ld2$quarter<-quarter(soes.ld2$date)
soes.q1<-subset(soes.ld2,quarter=="1")
soes.q2<-na.omit(subset(soes.ld2,quarter=="2"))
soes.q3<-subset(soes.ld2,quarter=="3")
soes.q4<-subset(soes.ld2,quarter=="4")


#jfm

d <- data.frame(x=soes.q1$x,y=soes.q1$y, tms=soes.q1$date, id=soes.q1$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/156.1

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Elephant_JFM.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Elephant_JFM.tif", overwrite=TRUE)

# amj

d <- data.frame(x=soes.q2$x,y=soes.q2$y, tms=soes.q2$date, id=soes.q2$idy)
tr <- trip(d,correct_all = T)
tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/64.1

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Elephant_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Elephant_AMJ.tif", overwrite=TRUE)

# jas

d <- data.frame(x=soes.q3$x,y=soes.q3$y, tms=soes.q3$date, id=soes.q3$idy)
tr <- trip(d,correct_all = T)
tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/46.2

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Elephant_JAS.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Elephant_JAS.tif", overwrite=TRUE)

#ond
d <- data.frame(x=soes.q4$x,y=soes.q4$y, tms=soes.q4$date, id=soes.q4$idy)
tr <- trip(d)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/59

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Seal_Elephant_OND.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Seal_Elephant_OND.tif", overwrite=TRUE)


####-----------Penguins---------------


# Chinstrap, Adelie and Gentoo penguin data from two sources: 

#Hinke et al. 2017 https://doi.org/10.1371/journal.pone.0170132
# Hinke et al. 2020 https://doi.org/10.1098/rsbl.2020.0645

# Adelie Penguin: RAATD https://zenodo.org/record/3722948#.YqG-DpBBw-Q

#Gentoo Penguins: Korczak-Abshire et al. 2021 https://doi.org/10.1098/rsbl.2020.0708

#Lload and process data from each paper

# ----------Hinke et al. 2017-----------

jh17<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/pone.0170132.s001/Data/satellite telemetry.csv")

head(jh17)

#summary(as.factor(jh17$Loc.Qual))

jh17<-subset(jh17,Loc.Qual!="Z")
jh17<-subset(jh17,Loc.Qual!="0")
#jh17<-subset(jh17,Loc.Qual!="")

jh17$timestamp<-as.POSIXct(strptime(paste(jh17$Date,jh17$Time), format="%m/%d/%Y %H:%M:%S", tz="GMT"))

jh17<-data.frame(Lon=jh17$Longitude,Lat=jh17$Latitude,id=as.factor(jh17$Deployment),
                 Spp=jh17$Spp,timestamp=jh17$timestamp,Col=jh17$Site)

# ---------Hinke et al. 2020-----------

jh20<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/FledglingBottleneck/tracks.csv")

jh20$timestamp<-as.POSIXct(strptime(jh20$LOC_DATE, 
                                    format="%m/%d/%Y %H:%M:%S", tz="GMT"))

jh20<-data.frame(Lon=jh20$LONGITUDE,Lat=jh20$LATITUDE,id=as.factor(jh20$PTT),
                 Spp=jh20$SPP,timestamp=jh20$timestamp,Col=c("COPA"))

jh20$Spp[jh20$Spp=="PYD"]<-"ADPE"
jh20$Spp[jh20$Spp=="PYN"]<-"CHPE"
jh20$Spp[jh20$Spp=="PYP"]<-"GEPE"

###---- Korczak-Abshire et al. 2021--------

ka21<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/Gentoo_Korczak-abshire/manuscript_data.csv")

ka21$timestamp<-as.POSIXct(strptime(ka21$LOC_DATE, 
                                    format="%m/%d/%Y %H:%M", tz="GMT"))

ka21<-data.frame(Lon=ka21$LONGITUDE,Lat=ka21$LATITUDE,id=as.factor(ka21$PTT),Spp=c("GEPE"),
                 timestamp=ka21$timestamp,Col=ka21$COLONY)

### -----RAATD--------

ra<-read.csv("C:/D1MPA_model_2024/Predators/Tracking/SCAR_EGBAMM_RAATD_2018_Standardised/RAATD_ADPE_standardized.csv")


ra$timestamp<-as.POSIXct(strptime(paste(paste(ra$year,ra$month,ra$day,sep="-"),
                                        ra$time), 
                                  format="%Y-%m-%d %H:%M:%S", tz="GMT"))

ra<-subset(ra,location_quality!="Z")

raat<-data.frame(Lon=ra$decimal_longitude,
                 Lat=ra$decimal_latitude,
                 id=as.factor(ra$individual_id),Spp=ra$abbreviated_name,
                 timestamp=ra$timestamp,Col=c("SOI"))

raatd<-subset(raat,Lon<0 & Lon>(-150))

df5<-rbind(jh17,jh20,ka21,raatd)

pengs<-subset(df5,Spp=="ADPE"|Spp=="CHPE"|Spp=="GEPE")

head(pengs)

peng<- pengs[!duplicated(pengs[, c("Lat", "Lon","timestamp","id")]), ]

pengm<-plyr::ddply(pengs, c("id","Spp","timestamp","Col"), summarise,
                   lon=mean(Lon),
                   lat=mean(Lat)) 

pengN<-plyr::ddply(pengs, c("id","Spp"), summarise,
                   N=length(Lon)) 

pengm<-merge(pengm,pengN)

pengm<-subset(pengm,lon<0)
pengm<-subset(pengm,lat<(-50))
pengm<-subset(pengm,N>5)


#write.csv(pengm,"C:/D1MPA_model_2024/Predators/Tracking/penguins.csv") 

pengm$timestamp<-as.POSIXct(strptime(pengm$timestamp,
                                     format="%Y-%m-%d %H:%M:%S", tz="GMT"))
pengm$month<-month(pengm$timestamp)
pengm$id<-paste(pengm$Spp,pengm$id,pengm$month)

#pengm<-subset(pengm,month<9)

pengm2<-plyr::ddply(pengm, c("id"), summarise,
                    N=length(lon)) 

pengms<-merge(pengm,pengm2,by="id")

pengms<-subset(pengms,N.y>5)

adp<-subset(pengms,Spp=="ADPE")
chp<-subset(pengms,Spp=="CHPE")
gep<-subset(pengms,Spp=="GEPE")

adp<-subset(adp,lon>(-75))

gep<-subset(gep,lon<(-50) & lon>(-67))
gep<-subset(gep,lat<(-61.5))

adsp<-SpatialPointsDataFrame(adp[5:6],adp,proj4string = EPSG4326)
chsp<-SpatialPointsDataFrame(chp[5:6],chp,proj4string = EPSG4326)
gesp<-SpatialPointsDataFrame(gep[5:6],gep,proj4string = EPSG4326)

adspW<-sp::spTransform(adsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
chspW<-sp::spTransform(chsp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
gespW<-sp::spTransform(gesp,CRS("+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

plot(adspW)
plot(chspW)
plot(gespW)

### --------gentoo penguin----------------


gepW<-raster::intersect(gespW,D1)

geps<-data.frame(gepW)

head(geps)
geps<- geps[!duplicated(geps[, c("timestamp","id")]), ]

gep.lt<-as.ltraj(xy=geps[16:17],date=geps$timestamp,id=geps$id,proj4string = EPSG6932)

gep.ld<-ld(gep.lt)

ggplot(gep.ld,aes(dist))+geom_histogram()

gep.ld2<-subset(gep.ld,dist<50000) #eliminate positions based on consecutive points distance 

gep.ld2$year<-year(gep.ld2$date)

#write.csv(gep.ld2,"C:/D1MPA trimestral/Data/Penguin_Gentoo.csv")


gep.ld2$quarter<-quarter(gep.ld2$date)
gep.q1<-subset(gep.ld2,quarter=="1")
gep.q2<-na.omit(subset(gep.ld2,quarter=="2"))
gep.q3<-subset(gep.ld2,quarter=="3")
gep.q4<-subset(gep.ld2,quarter=="4")

plot(gep.q2$x,gep.q2$y)

# amj

d <- data.frame(x=gep.q2$x,y=gep.q2$y, tms=gep.q2$date, id=gep.q2$id)
tr <- trip(d,correct_all = T)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/8.71

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Penguin_Gentoo_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Penguin_Gentoo_AMJ.tif", overwrite=TRUE)

# jas

d <- data.frame(x=gep.q3$x,y=gep.q3$y, tms=gep.q3$date, id=gep.q3$id)
tr <- trip(d,correct_all = T)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/2.1

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Penguin_Gentoo_JAS.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Penguin_Gentoo_JAS.tif", overwrite=TRUE)

### -----chinstrap-------

summary(as.factor(chp$Col))

chspW<-raster::intersect(chspW,D1)

chps<-data.frame(chspW)

chps<- chps[!duplicated(chps[, c("timestamp","id")]), ]

head(chps)

chp.lt<-as.ltraj(xy=chps[16:17],date=chps$timestamp,id=chps$id,proj4string = EPSG6932)

chp.ld<-ld(chp.lt)

ggplot(chp.ld,aes(dist))+geom_histogram()


chp.ld2<-subset(chp.ld,dist<50000) #eliminate positions based on consecutive points distance 

#write.csv(chp.ld2,"C:/D1MPA trimestral/Data/Penguin_Chinstrap.csv")

chp.ld2$quarter<-quarter(chp.ld2$date)
chp.q1<-subset(chp.ld2,quarter=="1")
chp.q2<-na.omit(subset(chp.ld2,quarter=="2"))
chp.q3<-subset(chp.ld2,quarter=="3")
chp.q4<-subset(chp.ld2,quarter=="4")

# amj

d <- data.frame(x=chp.q2$x,y=chp.q2$y, tms=chp.q2$date, id=chp.q2$id)

tr <- trip(d,correct_all = T)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/1.15

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Penguin_Chinstrap_AMJ.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Penguin_Chinstrap_AMJ.tif", overwrite=TRUE)

# jas

d <- data.frame(x=chp.q3$x,y=chp.q3$y, tms=chp.q3$date, id=chp.q3$id)
tr <- trip(d,correct_all = T)

tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/0.15

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Penguin_Chinstrap_JAS.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Penguin_Chinstrap_JAS.tif", overwrite=TRUE)

##----------adelie penguin--------------

summary(as.factor(adp$Col))

adspW<-raster::intersect(adspW,D1)

adps<-data.frame(adspW)


adps<- adps[!duplicated(adps[, c("timestamp","id")]), ]

head(adps)

adp.lt<-as.ltraj(xy=adps[16:17],date=adps$timestamp,id=adps$id,proj4string = EPSG6932)

adp.ld<-ld(adp.lt)

ggplot(adp.ld,aes(x,y))+geom_point()

ggplot(adp.ld,aes(dist))+geom_histogram()

adp.ld2<-subset(adp.ld,dist<50000) #eliminate positions based on consecutive points distance 

adp.ld2$month<-month(adp.ld2$date)

#write.csv(adp.ld2,"C:/D1MPA trimestral/Data/Penguin_Adelie.csv")

adspM<-subset(adp.ld2,month=="2"|month=="3"|month=="4"|month=="1")



d <- data.frame(x=adspM$x,y=adspM$y, tms=adspM$date, id=adspM$id)
tr <- trip(d,correct_all = T)
tg <- rasterize(tr, mask, field = "tms") # in seconds

y<-cover(tg,mask)/(3600)/24

max(y)

ym<-y/10.1

plot(ym)

#writeRaster(y,"C:/D1MPA trimestral/Trip/Penguin_Adelie_Mig.tif", overwrite=TRUE)
#writeRaster(ym,"C:/D1MPA trimestral/Trip/Rescaled/Penguin_Adelie_Mig.tif", overwrite=TRUE)
