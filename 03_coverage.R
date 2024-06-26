


rm(list=ls())

gc()


library(lubridate)
library(raster)
library(dplyr)
library(sp)
library(trip)
library(ggplot2)
library(patchwork)
library(terra)
library(reshape2)




###------------- rescale data to binary ----------------

# load files

# Set the directory where TIF files are located
tif_directory <- "C:/D1MPA trimestral/Trip/Rescaled"

# List all TIF files in the directory
tif_files <- list.files(tif_directory, pattern = ".tif$", full.names = T)

# Load all TIF files into a raster stack
raster_stack <- terra::rast(tif_files)

# Define the number of intervals
num_intervals <- 2


for (i in 1:nlyr(raster_stack)) {
  # Extract the ith layer
  raster_layer <- raster_stack[[i]]
  
  # Calculate the 0.05 quantile
  quantile_005 <- quantile(values(raster_layer), probs = 0.95,na.rm=T)
  
  # Reclassify the raster layer
  reclassified_layer <- terra::ifel(raster_layer > quantile_005, 1, 0)
  
  # Replace the original layer with the reclassified layer in the raster stack
  raster_stack[[i]] <- reclassified_layer
}


# Plot the reclassified raster stack
plot(raster_stack)

# load additional binary files 

# Set the directory where TIF files are located
tif_directory2 <- "C:/D1MPA 2024 EMM/Layers2024/Non-tracking_data"

# List all TIF files in the directory
tif_files2 <- list.files(tif_directory2, pattern = ".tif$", full.names = T)


# Load all TIF files into a raster stack
raster_stack2 <- terra::rast(tif_files2)
plot(raster_stack2)


#stack both stacks together
rasters0<-c(raster_stack,raster_stack2)

rsum<-sum(rasters0) # sum all stacks

plot(rsum)

max(rsum)

#writeRaster(rsum,"C:/D1MPA trimestral/D1MPA Fishable Area/UpDated_Layers_Sum.tif",overwrite=T)


### -------2018 layers------------

# load additional binary files 

# Set the directory where TIF files are located
tif_directory3 <- "C:/D1MPA trimestral/D1MPA Fishable Area/Layers 2018/geotiff"

# List all TIF files in the directory
tif_files3 <- list.files(tif_directory3, pattern = ".tif$", full.names = T)

tif_files3
# Load all TIF files into a raster stack
raster_stack3 <- terra::rast(tif_files3)
plot(raster_stack3)

sum2018<-sum(raster_stack3)

plot(sum2018)


#writeRaster(sum2018,"C:/D1MPA trimestral/D1MPA Fishable Area/Layers_2018_Sum.tif",overwrite=T)


#stack both stacks together
rasters<-c(raster_stack,raster_stack2,raster_stack3)

rsum<-sum(rasters) # sum all stacks

plot(rsum)

max(rsum)


#writeRaster(rsum,"C:/D1MPA trimestral/D1MPA Fishable Area/All_Layers_Sum.tif",overwrite=T)



###--- coverage --------

### load subareas and strata files

strata<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Subareas_Strata.shp")
plot(strata)

summary(as.factor(strata$SubArea))

gpz<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/GPZ_2024.shp")
head(gpz)
plot(gpz)

# preferred fishable area

fa<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Fishable_Area.shp")
plot(fa)

### spatial points

spt<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Grid_Points_D1.shp")

head(spt)

stfa<- (over(x = spt, y = fa))

summary(as.factor(stfa$Fishable))

stfa<-SpatialPointsDataFrame(coordinates(spt),stfa[4],proj4string = crs(spt))

#-----strata coverage -----------

st.points<- (over(x = stfa, y = strata))

summary(as.factor(st.points$SubArea))

stratas<-data.frame(st.points[4],stfa[1], coordinates(spt))

stgrid<-SpatialPointsDataFrame(coordinates(spt),stratas)

#plot(stgrid)

stvect<-vect(stgrid)

head(stvect)

summary(as.factor(stvect$SubArea))

ext.st<-terra::extract(rasters, stvect)

head(ext.st)

st.df <- data.frame(FA=stvect$Fishable,SA=stvect$SubArea,ext.st[2:55])

head(st.df)

summary(as.factor(st.df$SA))
summary(as.factor(st.df$FA))


###--------- gpz covergae-----

gpz.points<- (over(x = stfa, y = gpz))

gpzs<-data.frame(gpz.points)
head(gpzs)
summary(as.factor(gpzs$SubAreas))

summary(as.factor(gpzs$Strata))

sfdf<-data.frame(stfa,gpzs)

head(sfdf)

gpzgrid<-SpatialPointsDataFrame(coordinates(stfa),sfdf)

gpzvect<-vect(gpzgrid)

names(gpzvect)[names(gpzvect) == "SubAreas"] <- "SA"

ext.gpz<-terra::extract(rasters, gpzvect)

gpz.df <- data.frame(FA=sfdf$Fishable,SA=gpzvect$SA,ext.gpz[2:55])

#gpz.df$SA[is.na(gpz.df$SA)]<-"Outside"

summary(as.factor(gpz.df$FA))
summary(as.factor(gpz.df$SA))


###-----coverage calcualtion ---------

total<- data.frame(st.df[3:56]) %>%
  #group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<- data.frame(st.df[1],st.df[3:56]) %>%
  group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<-subset(totalF,FA=="Yes")


SAO <- st.df[2:56] %>%
  group_by(SA) %>%
  summarise_all(sum, na.rm = TRUE)

SAF <- st.df %>%
  group_by(FA,SA) %>%
  summarise_all(sum, na.rm = TRUE)

SAF<-subset(SAF,FA=="Yes")
SAO
SAF

GPZ <- gpz.df[2:56] %>%
  group_by(SA) %>%
  summarise_all(sum, na.rm = TRUE)

GPF <- gpz.df %>%
  group_by(FA,SA) %>%
  summarise_all(sum, na.rm = TRUE)

GPF<-subset(GPF,FA=="Yes")

totm<-melt(total)
tofm<-melt(totalF)

saom<-melt(SAO,id.vars=c("SA"))
safm<-melt(SAF,id.vars=c("FA","SA"))

gpzm<-na.omit(melt(GPZ,id.vars=c("SA")))
gpfm<-na.omit(melt(GPF,id.vars=c("FA","SA")))

names(totm)[names(totm) == "value"] <- "range"
names(saom)[names(saom) == "value"] <- "SA.cover"
names(gpzm)[names(gpzm) == "value"] <- "GPZ.cover"

names(tofm)[names(tofm) == "value"] <- "rangeF"
names(safm)[names(safm) == "value"] <- "SA.coverF"
names(gpfm)[names(gpfm) == "value"] <- "GPZ.coverF"

head(gpfm)

head(safm)

head(tofm)

gpz.st<-merge(gpzm,saom,all=T)

alldf<-merge(gpz.st,totm,all=T)
head(alldf)

alldff1<-merge(alldf,tofm)
alldff2<-merge(alldff1,safm)
alldf<-merge(alldff2,gpfm)

head(alldf)

names(alldf)[names(alldf) == "variable"] <- "layers"

alldf$sa.fa.prop<-alldf$SA.coverF/alldf$SA.cover

alldf$sa.fa.prop[is.na(alldf$sa.fa.prop)]<-0


### if a layer has more than 75% of its occurrence within fishable area, it is a priority


alldf$priority<-ifelse(alldf$sa.fa.prop>0.75,"Priority","Non-priority")

alldf$GPZ.prop<-alldf$GPZ.cover/alldf$SA.cover

alldf$GPZ.prop[is.na(alldf$GPZ.prop)]<-0



head(alldf)


#summer layers


# specify layers. If it starts with X they are krill for the updated layers
# if O, they are layers from 2018
# if other, they are updated predator layers

alldf$type<-substring(alldf$layers,first=1,last=1)

summary(as.factor(alldf$type))

alldf$types[alldf$type=="K"]<-"Krill updated"
alldf$types[alldf$type=="O"]<-"layers 2018"
alldf$types[is.na(alldf$types)]<-"predators updated"

alldf$types[alldf$layers=="obj_5b_Esuperba"]<-"Krill updated"

summary(as.factor(alldf$types))

layers_2018<-subset(alldf,types=="layers 2018")
predators<-subset(alldf,types=="predators updated")

krill<-subset(alldf,types=="Krill updated" )


gpzS<-subset(predators,layers!="Penguin_Gentoo_AMJ" & layers!="Penguin_Gentoo_JAS" &
               layers!="Penguin_Chinstrap_AMJ" & layers!="Penguin_Chinstrap_JAS"&
               layers!="Seal_Crabeater_AMJ"&layers!="Seal_Crabeater_JAS"&
               layers!="Seal_Elephant_AMJ"&layers!="Seal_Elephant_JAS"&
               layers!="Seal_Fur_AMJ"&layers!="Seal_Fur_JAS"&
               layers!="Whale_Humpback_AMJ"&layers!="Whale_AntMinke_AMJ"&
               layers!="Penguin_Emperor_Cols" & layers!="obj_5b_Esuperba" & 
               layers!="AFS_BS_SOI" )
#winter layers


gpzW<-subset(predators,layers=="Penguin_Gentoo_AMJ" | layers=="Penguin_Gentoo_JAS" |
               layers=="Penguin_Chinstrap_AMJ" | layers=="Penguin_Chinstrap_JAS"|
               layers=="Seal_Crabeater_AMJ"|layers=="Seal_Crabeater_JAS"|
               layers=="Seal_Elephant_AMJ"|layers=="Seal_Elephant_JAS"|
               layers=="Seal_Fur_AMJ"|layers=="Seal_Fur_JAS"|
               layers=="Whale_Humpback_AMJ"|layers=="Whale_AntMinke_AMJ"|
               layers=="Penguin_Emperor_Cols")



#winter layers

# 2018 layers 

(ggplot(layers_2018,aes(layers,GPZ.prop,colour=priority))+
    geom_segment( aes(x=layers, xend=layers, y=0.5, yend=GPZ.prop),
                  linetype="dotted") +
    geom_point(size=2)+
    geom_hline(yintercept=c(0.5))+
    facet_wrap(SA~.)+coord_flip()+
    scale_colour_manual(values=c("blue","red"))+
    ggtitle(label="")+theme_bw()+
    xlab("layers")+ylab("Proportion of GPZ coverage")+
    ggtitle(label="Layers 2019"))

### layers 2024 

(ggplot((krill),aes(layers,GPZ.prop,colour=priority))+
    geom_segment( aes(x=layers, xend=layers, y=0.2, yend=GPZ.prop),
                  linetype="dotted") +
    geom_point(size=2)+
    geom_hline(yintercept=c(0.2))+
    facet_wrap(SA~.)+coord_flip()+
    scale_colour_manual(values=c("blue","red"))+
    ggtitle(label="")+theme_bw()+
    xlab("Summer layers")+ylab("Proportion of GPZ coverage")+
    ggtitle(label="a. Layers 2024 Krill"))/

# summer layers

(ggplot((gpzS),aes(layers,GPZ.prop,colour=priority))+
    geom_segment( aes(x=layers, xend=layers, y=0.5, yend=GPZ.prop),
                  linetype="dotted") +
    geom_point(size=2)+
    geom_hline(yintercept=c(0.5))+
    facet_wrap(SA~.)+coord_flip()+
    scale_colour_manual(values=c("blue","red"))+
    ggtitle(label="")+theme_bw()+
    xlab("Summer layers")+ylab("Proportion of GPZ coverage")+
   ggtitle(label="b. Layers 2024 October to March"))/

# a part of fish habitat is inside SOISS MPA. That might increase coverage



#Winter layers

(ggplot((gpzW),aes(layers,GPZ.prop,colour=priority))+
    geom_segment( aes(x=layers, xend=layers, y=0.5, yend=GPZ.prop),
                  linetype="dotted") +
    geom_point(size=2)+
    geom_hline(yintercept=c(0.5))+
    facet_wrap(SA~.)+coord_flip()+
    scale_colour_manual(values=c("blue","red"))+
    ggtitle(label="")+theme_bw()+
    xlab("Winter layers")+ylab("Proportion of GPZ coverage")+
    ggtitle(label="Layers 2024 April to September"))




### now, do the same with the SPZ considered



### load subareas and strata files

strata<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Subareas_Strata.shp")
plot(strata)

summary(as.factor(strata$SubArea))

gpz<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/GPZ_SPZ_cut.shp")
head(gpz)
plot(gpz)
summary(as.factor(gpz$SubArea))

# preferred fishable area

fa<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Fishable_Area.shp")
plot(fa)

### spatial points

spt<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Grid_Points_D1.shp")

head(spt)

stfa<- (over(x = spt, y = fa))

summary(as.factor(stfa$Fishable))

stfa<-SpatialPointsDataFrame(coordinates(spt),stfa[4],proj4string = crs(spt))

#-----strata coverage -----------

st.points<- (over(x = stfa, y = strata))

summary(as.factor(st.points$SubArea))

stratas<-data.frame(st.points[4],stfa[1], coordinates(spt))

stgrid<-SpatialPointsDataFrame(coordinates(spt),stratas)

#plot(stgrid)

stvect<-vect(stgrid)

head(stvect)

summary(as.factor(stvect$SubArea))

ext.st<-terra::extract(rasters, stvect)

head(ext.st)

st.df <- data.frame(FA=stvect$Fishable,SA=stvect$SubArea,ext.st[2:61])

head(st.df)

summary(as.factor(st.df$SA))
summary(as.factor(st.df$FA))


###--------- gpz covergae-----

gpz.points<- (over(x = stfa, y = gpz))

gpzs<-data.frame(gpz.points)
head(gpzs)

summary(as.factor(gpzs$SubArea))
head(stfa)
sfdf<-data.frame(stfa,gpzs)

head(sfdf)
summary(as.factor(sfdf$SubArea))

gpzgrid<-SpatialPointsDataFrame(coordinates(stfa),sfdf)

gpzvect<-vect(gpzgrid)

names(gpzvect)[names(gpzvect) == "SubArea"] <- "SA"

summary(as.factor(gpzvect$SA))

ext.gpz<-terra::extract(rasters, gpzvect)

gpz.df <- data.frame(FA=sfdf$Fishable,SA=gpzvect$SA,ext.gpz[2:61])

#gpz.df$SA[is.na(gpz.df$SA)]<-"Outside"

summary(as.factor(gpz.df$FA))
summary(as.factor(gpz.df$SA))


###-----coverage calcualtion ---------

total<- data.frame(st.df[3:62]) %>%
  #group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<- data.frame(st.df[1],st.df[3:62]) %>%
  group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<-subset(totalF,FA=="Yes")


SAO <- st.df[2:62] %>%
  group_by(SA) %>%
  summarise_all(sum, na.rm = TRUE)

SAF <- st.df %>%
  group_by(FA,SA) %>%
  summarise_all(sum, na.rm = TRUE)

SAF<-subset(SAF,FA=="Yes")
SAO
SAF

GPZ <- gpz.df[2:62] %>%
  group_by(SA) %>%
  summarise_all(sum, na.rm = TRUE)

GPF <- gpz.df %>%
  group_by(FA,SA) %>%
  summarise_all(sum, na.rm = TRUE)

GPF<-subset(GPF,FA=="Yes")


GPZ
GPF

totm<-melt(total)
tofm<-melt(totalF)

saom<-melt(SAO,id.vars=c("SA"))
safm<-melt(SAF,id.vars=c("FA","SA"))

gpzm<-na.omit(melt(GPZ,id.vars=c("SA")))
gpfm<-na.omit(melt(GPF,id.vars=c("FA","SA")))
summary(as.factor(gpzm$SA))
summary(as.factor(gpfm$SA))

names(totm)[names(totm) == "value"] <- "range"
names(saom)[names(saom) == "value"] <- "SA.cover"
names(gpzm)[names(gpzm) == "value"] <- "GPZ.cover"

names(tofm)[names(tofm) == "value"] <- "rangeF"
names(safm)[names(safm) == "value"] <- "SA.coverF"
names(gpfm)[names(gpfm) == "value"] <- "GPZ.coverF"

head(gpfm)

head(safm)

head(tofm)

gpz.st<-merge(gpzm,saom,all=T)

alldf<-merge(gpz.st,totm,all=T)
head(alldf)

alldff1<-merge(alldf,tofm)
alldff2<-merge(alldff1,safm)
alldf<-merge(alldff2,gpfm)

head(alldf)

names(alldf)[names(alldf) == "variable"] <- "layers"

alldf$sa.fa.prop<-alldf$SA.coverF/alldf$SA.cover

alldf$sa.fa.prop[is.na(alldf$sa.fa.prop)]<-0


### if a layer has more than 75% of its occurrence within fishable area, it is a priority


alldf$priority<-ifelse(alldf$sa.fa.prop>0.75,"Priority","Non-priority")

alldf$GPZ.prop<-alldf$GPZ.cover/alldf$SA.cover

alldf$GPZ.prop[is.na(alldf$GPZ.prop)]<-0



head(alldf)


#summer layers


# specify layers. If it starts with X they are krill for the updated layers
# if O, they are layers from 2018
# if other, they are updated predator layers

alldf$type<-substring(alldf$layers,first=1,last=1)

summary(as.factor(alldf$type))

alldf$types[alldf$type=="X"]<-"Krill updated"
alldf$types[alldf$type=="O"]<-"layers 2018"
alldf$types[is.na(alldf$types)]<-"predators updated"

alldf$types[alldf$layers=="obj_5b_Esuperba"]<-"Krill updated"

summary(as.factor(alldf$types))



summary(as.factor(alldf$priority))

prior<-subset(alldf,priority=="Priority")

krillSPZ<-subset(prior,types=="Krill updated" )
predators<-subset(prior,types=="predators updated")

gpzSPZ<-subset(predators,layers!="Penguin_Gentoo_AMJ" & layers!="Penguin_Gentoo_JAS" &
               layers!="Penguin_Chinstrap_AMJ" & layers!="Penguin_Chinstrap_JAS"&
               layers!="Seal_Crabeater_AMJ"&layers!="Seal_Crabeater_JAS"&
               layers!="Seal_Elephant_AMJ"&layers!="Seal_Elephant_JAS"&
               layers!="Seal_Fur_AMJ"&layers!="Seal_Fur_JAS"&
               layers!="Whale_Humpback_AMJ"&layers!="Whale_AntMinke_AMJ"&
               layers!="Penguin_Emperor_Cols" & layers!="obj_5b_Esuperba" & 
                 layers!="AFS_BS_SOI" & layers!="Seal_Crabeater_OND" )

gpzS1<-subset(gpzS,priority=="Priority")

krillS1<-subset(krill,priority=="Priority")

krillS1$unit<-c("GPZ")
krillSPZ$unit<-c("GPZ + SPZ")

gpzS1$unit<-c("GPZ")
gpzSPZ$unit<-c("GPZ + SPZ")


krilldf<-rbind(krillS1,krillSPZ)


gpzspz2<-rbind(gpzS1,gpzSPZ,subset(krilldf,layers=="X01_Krill_JR_Weddell"))
gpzspz2<-subset(gpzspz2,layers!="Seal_Crabeater_OND")


### layers 2024 GPZ and GPZ+SPZ
summary(as.factor(gpzspz2$layers))

gpzspz2$layers<-factor(gpzspz2$layers,levels=c("bathy_150_500","bathy_0_150",
                                               "Whale_Humpback_JFM",
                                               "Whale_AntMinke_JFM","Seal_Elephant_JFM",
                                              "IBA_CHP", "IBA_GEP","IBA_ADP",
                                               "X01_Krill_JR_Weddell"))
  
  (ggplot(subset(gpzspz2,SA=="SA 48.1"),aes(layers,GPZ.prop,colour=unit))+
     geom_segment( aes(x=layers, xend=layers, y=0.5, yend=GPZ.prop),
                   linetype="dotted") +
     geom_point(size=2)+
     geom_hline(yintercept=c(0.5),linetype="dashed")+
     facet_wrap(SA~.)+coord_flip()+
     scale_colour_manual(values=c("blue","red"))+
     ggtitle(label="")+theme_bw()+
     xlab("Summer layers")+ylab("Proportion")+
     ggtitle(label=" layers GPZ + SPZ coverage"))
  
  

### ------- GPZ units coverage----------



###--- coverage --------

### load subareas and strata files

gpz<-shapefile("C:/D1MPA 2024 EMM/Vertices/D1MPA_Vertices_2024.shp")
head(gpz)
plot(gpz)

# preferred fishable area

fa<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Fishable_Area.shp")
plot(fa)

### spatial points

spt<-shapefile("C:/D1MPA trimestral/D1MPA Fishable Area/Shapefiles/Grid_Points_D1.shp")

head(spt)

stfa<- (over(x = spt, y = fa))

summary(as.factor(stfa$Fishable))

stfa<-SpatialPointsDataFrame(coordinates(spt),stfa[4],proj4string = crs(spt))

#-----strata coverage -----------

st.points<- (over(x = stfa, y = strata))

summary(as.factor(st.points$SubArea))

stratas<-data.frame(st.points[4],stfa[1], coordinates(spt))

stgrid<-SpatialPointsDataFrame(coordinates(spt),stratas)

#plot(stgrid)

stvect<-vect(stgrid)

head(stvect)

summary(as.factor(stvect$SubArea))

ext.st<-terra::extract(rasters0, stvect)

head(ext.st)

st.df <- data.frame(FA=stvect$Fishable,SA=stvect$SubArea,ext.st[2:32])

head(st.df)

summary(as.factor(st.df$SA))
summary(as.factor(st.df$FA))


###--------- gpz covergae-----

head(gpz)



gpz.points<- (over(x = stfa, y = gpz))

gpzs<-data.frame(gpz.points)
head(gpzs)

summary(as.factor(gpzs$Name))

#summary(as.factor(gpzs$Strata))

sfdf<-data.frame(stfa,gpzs)

head(sfdf)

gpzgrid<-SpatialPointsDataFrame(coordinates(stfa),sfdf)

gpzvect<-vect(gpzgrid)

#names(gpzvect)[names(gpzvect) == "SubAreas"] <- "SA"

ext.gpz<-terra::extract(rasters0, gpzvect)

gpz.df <- data.frame(FA=sfdf$Fishable,SA=gpzvect$Name,Zone=gpzvect$Zone,ext.gpz[2:32])

#gpz.df$SA[is.na(gpz.df$SA)]<-"Outside"

summary(as.factor(gpz.df$FA))
summary(as.factor(gpz.df$SA))



###-----coverage calcualtion ---------

total<- data.frame(st.df[3:56]) %>%
  #group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<- data.frame(st.df[1],st.df[3:56]) %>%
  group_by(FA) %>%
  summarise_all(sum, na.rm = TRUE)

totalF<-subset(totalF,FA=="Yes")

head(gpz.df)


GPZ <- gpz.df[2:34] %>%
  group_by(SA,Zone) %>%
  summarise_all(sum, na.rm = TRUE)

GPF <- gpz.df %>%
  group_by(FA,SA,Zone) %>%
  summarise_all(sum, na.rm = TRUE)

GPF<-subset(GPF,FA=="Yes")



totm<-melt(total)
tofm<-melt(totalF)

gpzm<-na.omit(melt(GPZ,id.vars=c("SA","Zone")))
gpfm<-na.omit(melt(GPF,id.vars=c("FA","SA","Zone")))




names(totm)[names(totm) == "value"] <- "range"
names(gpzm)[names(gpzm) == "value"] <- "cover"

names(tofm)[names(tofm) == "value"] <- "rangeF"
names(gpfm)[names(gpfm) == "value"] <- "coverF"

head(gpfm)

head(tofm)


df<-merge(gpzm,totm,all=T)
head(df)
dff<-merge(df,tofm)
head(dff)

alldf<-merge(dff,gpfm)


head(alldf)

names(alldf)[names(alldf) == "variable"] <- "layers"

alldf$propcov<-alldf$cover/alldf$range
alldf$propcovF<-alldf$coverF/alldf$rangeF
alldf$pro.rangeF<-alldf$rangeF/alldf$range

summary(as.factor(alldf$FA))

gpz.cov<-subset(alldf,Zone=="GPZ")

gpz.cov<-subset(gpz.cov,propcov>0)

write.csv(gpz.cov,"C:/D1MPA 2024 EMM/gpz_coverage.csv")


ggplot(subset(gpz.cov,SA=="SOI"),aes(reorder(layers,+propcovF), propcovF))+
  geom_point(size=3)+coord_flip()+facet_wrap(SA~.)+theme_bw()+
  xlab("Layers")+ylab("Proportion of coverage")+
  
  ggplot(subset(gpz.cov,SA=="EI"),aes(reorder(layers,+propcovF), propcovF))+
  geom_point(size=3)+coord_flip()+facet_wrap(SA~.)+theme_bw()+
  xlab("Layers")+ylab("Proportion of coverage")+
  
  ggplot(subset(gpz.cov,SA=="SSIW"),aes(reorder(layers,+propcovF), propcovF))+
  geom_point(size=3)+coord_flip()+facet_wrap(SA~.)+theme_bw()+
  xlab("Layers")+ylab("Proportion of coverage")+
  
  ggplot(subset(gpz.cov,SA=="NWAP"),aes(reorder(layers,+propcovF), propcovF))+
  geom_point(size=3)+coord_flip()+facet_wrap(SA~.)+theme_bw()+
  xlab("Layers")+ylab("Proportion of coverage")+
  
  
  ggplot(subset(gpz.cov,SA=="SWAP"),aes(reorder(layers,+propcovF), propcovF))+
  geom_point(size=3)+coord_flip()+facet_wrap(SA~.)+theme_bw()+
  xlab("Layers")+ylab("Proportion of coverage")

  
