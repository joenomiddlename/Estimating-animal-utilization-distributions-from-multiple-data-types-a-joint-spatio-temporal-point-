# July and August 2014 Analysis

# Mac version 
setwd("~/Documents/Whale Project")

# PC version
setwd("~/ownCloud/Whale Project")

library(dplyr)
library(lubridate)
library(maptools)
library(rgeos)
library(rgdal)

##
longform.points <- NULL
for (year in c("09","10","11","12","13","14","15","16")) {
# Mac edition  
#set.working.dir <- paste("/Users/joe/Documents/Whale Project/Data for Joe/CRP-BG EFFORT DATA/20",year,"/DAYLIGHT_trackpoints",sep="")
# PC verison
set.working.dir <- paste("~/ownCloud/Whale Project/Data for Joe/CRP-BG EFFORT DATA/20",year,"/DAYLIGHT_trackpoints",sep="")
setwd(set.working.dir)
temp <- list.files(pattern="*.dbf", full.names=TRUE)
filenames <- unlist(lapply(strsplit(unlist(strsplit(temp,".dbf")),"./"),"[[",2))
if (year=="09") ind=(1:length(filenames))[-c(36,101)] # take out the run up the east side of Vanc.Isl (36)
if (year!="09") ind=(1:length(filenames))
ind=(1:length(filenames))
for (i in ind) {
  temp.files=readOGR(".",layer=filenames[i])
  #assign(paste(substr(filenames[i],7,28),"_", i, sep=""), temp.files)
  
    longform.points <- rbind(longform.points , 
                             data.frame(x= temp.files@coords[,1],y= temp.files@coords[,2], 
                                        Date=temp.files@data$Day,
                                        Month = as.integer(substr(as.character(temp.files@data$Day),6,7)),
                                        Time=substr(as.character(temp.files@data$Date_PDT),12,19),
                                        Day = substr(as.character(temp.files@data$Date_PDT),9,10)) )
  
  if ((i%%5)==0)	print(paste("year",year,"trackline",i,sep=" - "))
  print(i)
  print('out of')
  print(length(ind))
}
}

#saveRDS(longform.points, file = 'all_trackpoints.rds')

#load( 'JulAug_dat.RData' )
setwd("~/ownCloud/Whale Project/Final Project Scripts")
longform.points = readRDS('all_trackpoints.rds') # no need to run above

# MAC
#setwd("/Users/joe/Documents/Whale Project/Data for Joe/Coast shapefile")
# PC
setwd("~/ownCloud/Whale Project/Data for Joe/Coast shapefile")
COAST=readOGR(".",layer="coast_ALBERS")

# plot the coastline for the area of interest
dev.off()
quartz(height=6.3,width=9.5) 
par(mfrow=c(1,1),oma=c(5.5,5.5,3,3),mar=rep(0,4))
plot(NA, 
     xlim=c(800000,1300000),
     ylim=c(350000,650000),axes=F)	 # right! # SRKW.s@extent

mtext(side=1,"BC Albers",outer=F,line=2.4,cex=.95);box()

for (i in c(18001:length(COAST@polygons))) { # coastline at < 18000 is north coast and not of interest
  for (j in 1:length(unlist(COAST@polygons[[i]]@Polygons))) {
    if (coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,2]>229000 &
        coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,2]<659000 &
        coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,1]>822391 &
        coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,1]<1300000 &
        length(coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,2])>75) {
      polygon(coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,1],coordinates(unlist(COAST@polygons[[i]]@Polygons)[[j]])[,2],
              border="black",lwd=1,col=rgb(0,0,0,.2))
    }}}# BC mainland
polygon(coordinates(unlist(COAST@polygons[[1]]@Polygons)[[1]])[,1],coordinates(unlist(COAST@polygons[[1]]@Polygons)[[1]])[,2],
        border="black",lwd=1,col=rgb(0,0,0,.2))
axis(1)
axis(2, las=2)
box()

# coordinates for the 'area of interest'
polygon(c(970,1270,1270,970)*1000,c(350,350,550,550)*1000, lwd=2, border='green3', col=NA)

setwd("~/ownCloud/Whale Project/Final Project Scripts")

##### USES LATER CREATED FILE
# Are the BCCSN WW sightings covering a new area?
WW_BCCSN = readRDS('WW_BCCSN_labelled_2001-2011.rds')
proj4string(WW_BCCSN) = '+proj=longlat +ellps=WGS84'
WW_BCCSN = spTransform(WW_BCCSN, COAST@proj4string)
points(WW_BCCSN@coords[,1],WW_BCCSN@coords[,2],  
       pch=19, cex=1, col='red')
points(WW_BCCSN[WW_BCCSN$YR>=2009,]@coords[,1],WW_BCCSN[WW_BCCSN$YR>=2009,]@coords[,2],  
       pch=19, cex=1, col='orange')
# Only make up ~3% of all sightings in 2009 onwards...
#####

# bring in the "focal follow" with-whale transects
# MAC
#BGdata <- read.csv("/Users/joe/Documents/Whale Project/Data for Joe/DFO CRP-BG LOCATIONS_BCalbers.csv")
# PC
BG <- read.csv("~/ownCloud/Whale Project/Data for Joe/DFO CRP-BG LOCATIONS_BCalbers.csv")
head(BG);tail(BG);dim(BG); "22893  x  23"

# Plot BGs effort points 
points(	longform.points[,c("x","y")],col='blue',cex=0.2,pch=16)
points(	longform.points[,c("x","y")],col='green',cex=.05,pch=16)

points(	BG[,c("Alb_East","Alb_North")],col=rgb(.8,0,0,.01),cex=1.2,pch=16)
points(	BG[,c("Alb_East","Alb_North")],col='orangered',cex=.2,pch=16)

# Check for duplicates in the tracks and BG observations #
sum(paste(longform.points$Day,paste(longform.points$Month, as.character(longform.points$Time))) %in% paste(BG$DAY, paste(BG$MONTH, as.character(BG$TIME_PDT))) )
# 17292 / 1341523 are duplicated in longform.points so remove them.

longform.points = longform.points[!paste(longform.points$Day,paste(longform.points$Month, as.character(longform.points$Time))) %in% paste(BG$DAY, paste(BG$MONTH, as.character(BG$TIME_PDT))),]
# Check for duplicates in the tracks and BG observations #
sum(paste(longform.points$Day,paste(longform.points$Month, as.character(longform.points$Time))) %in% paste(BG$DAY, paste(BG$MONTH, as.character(BG$TIME_PDT))) )
#

# Load whale watch data
# MAC
#BCCSNAlbers <- read.csv("~/Documents/Whale Project/SmrCombinedduplicatesRemoved_BCCSNAlbers.csv")
# PC
BCCSNAlbers <- read.csv("~/ownCloud/Whale Project/SmrCombinedduplicatesRemoved_BCCSNAlbers.csv")

# Now we wish to subset the data into only July 2009-2016 and include only WW presence-only data
WW = BCCSNAlbers[BCCSNAlbers$Source ==  'TWM-SA-WW' & BCCSNAlbers$Year >= 2009,]
WW = SpatialPointsDataFrame(coords = cbind(WW$lon,WW$lat),
                                      data = WW,
                                      proj4string = CRS(proj4string(COAST)))
# Plot WW sightings 

points(WW$lon,WW$lat,  
       pch=16, cex=1.2, col=rgb(1,1,0,.02))
points(coordinates(WW)[,1],coordinates(WW)[,2],  
       pch=16, cex=.2, col=rgb(1,.8,0,.2))

# Add a month and day variable to the tracks data
BG$lat = BG$Alb_North
BG$lon = BG$Alb_East
WW$MONTH = WW$Month
WW$YEAR = WW$Year
WW$DAY = WW$Day
longform.points$DAY = as.integer(longform.points$Day)
longform.points = longform.points[,-which(names(longform.points)== 'Day')]
BG$TIME = hms(as.character(BG$TIME_PDT))
day(BG$TIME) = as.numeric(BG$DAY)
month(BG$TIME) = as.numeric(BG$MONTH)
year(BG$TIME) = as.numeric(BG$YEAR)
longform.points$TIME = hms(as.character(longform.points$Time))
day(longform.points$TIME) = as.numeric(longform.points$DAY)
month(longform.points$TIME) = as.numeric(longform.points$Month)
year(longform.points$TIME) = year(longform.points$Date)

table(as.numeric(seconds(BG$TIME[1:(dim(BG)[1]-1)] - BG$TIME[2:(dim(BG)[1])])))
# Not all readings are 16 seconds apart! Continuous time Random Walk needed

#save.image(file = 'final_2009_2016.RData') # save all of the above data wrangling

WW_sightings = spTransform(WW, CRS = proj4string(COAST))

# Load raster variables (e.g. sea-surface temperatures)
#install.packages("httr", dependencies = TRUE)
#install.packages("ncdf4",dependencies = TRUE) 
#install.packages("readr", dependencies = TRUE)

area_of_interest = SpatialPoints(coords = cbind(c(970,1270,1270,970)*1000,c(350,350,550,550)*1000),
                                    proj4string = CRS(proj4string(COAST)))
area_of_interest = spTransform(area_of_interest, CRSobj = CRS("+proj=longlat"))

polygon_of_interest = Polygon(coords = cbind(c(970,1270,1270,970)*1000,c(350,350,550,550)*1000))
polygon_of_interest = Polygons(list(polygon_of_interest),1)
polygon_of_interest = SpatialPolygons(list(polygon_of_interest))
plot(polygon_of_interest)
proj4string(polygon_of_interest) = CRS(proj4string(COAST))

require("xtractomatic")
# First extract the topological features (e.g. sea depth)
tpos <- c('2006-04-01', '2016-11-30') 
ypos <- area_of_interest@bbox[2,]
xpos <- area_of_interest@bbox[1,]
topo <- xtracto_3D("ETOPO180", xpos, ypos, tpos)
hist(topo$data)
image(topo$data[,,1], col = heat.colors(100))
topo2 = topo
topo2$data[topo$data[,,1] >=0] = NA
image(topo2$data[,,1], col = heat.colors(100))
#getInfo('ETOPO18')

# extract chlorophyll
chloro <- xtracto_3D("mhchlamday", xpos, ypos, tpos)
getInfo("mhchlamday")
hist(chloro$data)
image(log(chloro$data[,,1]), col = heat.colors(100))
# get the more recent high res version 
tpos3 = c('2012-01-15','2016-11-30') 
chloro2 <- xtracto_3D("erdVH3chlamday", xpos, ypos, tpos3)
getInfo("erdVH3chlamday")
hist(chloro2$data)
image(log(chloro2$data[,,1]), col = heat.colors(100))

######## Important code for converting xtractomatic
######## data into SpatialPixelsDataFrame (needed for)
######## INLAbru
topo = SpatialPixelsDataFrame(points = expand.grid(topo$longitude, topo$latitude), 
                                data = as.data.frame(array(topo$data, dim = c(dim(topo$data)[1]*dim(topo$data)[2],dim(topo$data)[3]))), 
                                proj4string = CRS("+init=epsg:3857"), 
                                tolerance = 0.000183117)
chloro = SpatialPixelsDataFrame(points = expand.grid(chloro$longitude, chloro$latitude), 
                       data = as.data.frame(array(chloro$data, dim = c(dim(chloro$data)[1]*dim(chloro$data)[2],dim(chloro$data)[3]))), 
                       proj4string = CRS("+init=epsg:3857"), 
                       tolerance = 0.000183117)
plot(chloro['V3'])
chloro2 = SpatialPixelsDataFrame(points = expand.grid(chloro2$longitude, chloro2$latitude), 
                                data = as.data.frame(array(chloro2$data, dim = c(dim(chloro2$data)[1]*dim(chloro2$data)[2],dim(chloro2$data)[3]))), 
                                proj4string = CRS("+init=epsg:3857"), 
                                tolerance = 0.000183117)

# Next extract sea surface temperatures
searchData(list('varname:sst', 'datasetname:amday'))
tpos1 <- c('2006-04-01', '2007-04-15')
tpos2 <- c('2007-04-16', '2016-11-30')
SST_1 <- xtracto_3D("agsstamday", xpos, ypos, tpos1)
SST_2 <- xtracto_3D("atsstamday", xpos, ypos, tpos2)
par(mfrow = c(1,1))
hist(SST_1$data)
hist(SST_2$data)
range(SST_2$data)
range(SST_2$data, na.rm=TRUE)
sum(SST_2$data > 0 & SST_2$data < 25 & !is.na(SST_2$data)) / sum(!is.na(SST_2$data)) # 99.8% plausible readings
# Do the implausible values lie on land?
SST_2_2 = SST_2
SST_2_2$data[SST_2$data > 0 & SST_2$data < 30 & !is.na(SST_2$data)] = 1000 # give massive value to unsensible areas
image(apply(SST_2_2$data, c(1,2), sum, na.rm=T)) # red denotes missing (i.e. land)
# white denotes 100% data quality 
# yellow - orange denotes some quality issues
# It is clear that the violating areas lie on land and close to shore

# Plot all plausible SST values
SST_2_3 = SST_2
SST_2_3$data[(SST_2_3$data < 0 | SST_2_3$data > 30) & !is.na(SST_2_3$data)] = NA
image(SST_2_3$data[,,1])
image(SST_2_3$data[,,2])
image(SST_2_3$data[,,3])
image(SST_2_3$data[,,4])

## The SST measurements around the San Juan islands are very noisy 

# Missings exist in SST2 and outliers exist (SST can't be negative or 123 degrees...)
image(SST_1$data[,,1], col = heat.colors(100))
image(SST_2$data[,,1], col = heat.colors(100))
image(SST_2$data[,,2], col = heat.colors(100))
image(SST_2$data[,,3], col = heat.colors(100))

SST_1 = SpatialPixelsDataFrame(points = expand.grid(SST_1$longitude, SST_1$latitude), 
                              data = as.data.frame(array(SST_1$data, dim = c(dim(SST_1$data)[1]*dim(SST_1$data)[2],dim(SST_1$data)[3]))), 
                              proj4string = CRS("+init=epsg:3857"), 
                              tolerance = 0.000183117)
SST_2 = SpatialPixelsDataFrame(points = expand.grid(SST_2$longitude, SST_2$latitude), 
                               data = as.data.frame(array(SST_2$data, dim = c(dim(SST_2$data)[1]*dim(SST_2$data)[2],dim(SST_2$data)[3]))), 
                               proj4string = CRS("+init=epsg:3857"), 
                               tolerance = 0.000183117)
plot(SST_1['V1'])
plot(SST_2['V1'])


# calculate the sampling effort in the bounded box per month by Brian
# first compute delta t
table(as.numeric(seconds(longform.points$TIME[1:(dim(longform.points)[1]-1)] - longform.points$TIME[2:(dim(longform.points)[1])])))
difftime(longform.points$TIME)
sum(as.numeric(seconds(longform.points$TIME[1:(dim(longform.points)[1]-1)] - longform.points$TIME[2:(dim(longform.points)[1])])) < -100)
# loop through dates and calculate delta t's for ALL times that fall within bounding box

longform.points$dt = 0
longform.points$interior = 0
count = 0
for(temp in levels(longform.points$Date))
{
  ind = longform.points$Date == temp # subset by Date
  ind2 = gContains(polygon_of_interest,
                   SpatialPoints(coords = cbind(longform.points[ind,]$x,longform.points[ind,]$y),
                                 proj4string = CRS(proj4string(COAST))),
                   byid = T) # Which points lie within the area of interest
  #ind3 = ind*ind2 # refine the indices to include only those lying within the box
  tmp = longform.points[ind,]
  dt = as.numeric(seconds(tmp$TIME[2:(dim(tmp)[1])] - tmp$TIME[1:(dim(tmp)[1]-1)]))
  longform.points[ind,'dt'] = c(0,dt)
  longform.points[ind,'interior'] = as.numeric(ind2)
  count = count+1
  print(816-count)
}

# calculate approx total search/effort time spent in area of interest by adding up the dt's in the interior
# calculate the monthly search effort
BG_effort_monthly = matrix(0, nrow = 8, ncol = 12)
for(i in 1:12){
  for(j in 1:8){
    if(sum(longform.points$TIME@year == 2008+j & longform.points$TIME@month == i) > 0)
    {
      tmp = longform.points[longform.points$TIME@year == 2008+j & longform.points$TIME@month == i,]
      BG_effort_monthly[j,i] = sum(tmp[tmp$interior==1,]$dt)
    }
  }
}
BG_effort_monthly # in seconds
# convert to hours
BG_effort_monthly = BG_effort_monthly/3600
BG_effort_monthly

# save all the above calculations
save.image('final_2009_2016_V2.RData')

## How do the composition of pods change throughout the region?
rowSums(xtabs(~ PODS + DET_ID, data = BG)>0) # 94% of Brian's sightings contain L pod ((104+18+7+7)/(104+18+7+7+4+3+2))
table(WW@data$Pod) # 43% of sightings made by WW companies contain L pod (473+287+169+707)/(sum(table(WW$Pod))-1443)

# Get vessel density
library(raster)
setwd("~/ownCloud/Whale Project/test_shapefile/Smoothed")
rlist=list.files(getwd(), pattern="tif$", full.names=FALSE)
#vessel_density = raster(rlist)
#plot(vessel_density)
#crs(vessel_density) = CRS("+init=epsg:3857")
vessel_density2 = readGDAL(rlist[1])
plot(vessel_density2)
crs(vessel_density2) = CRS("+init=epsg:3857")

# save all the above calculations
save.image('final_2009_2016_V3.RData')

######### How to convert vessel density into SpatialPixelsDataFrame
######### Needed for INLAbru

# vessel_density3 = SpatialPixelsDataFrame(data = vessel_density2@data,
#                                          grid = vessel_density2@grid,
#                                          points = coordinates(vessel_density2@grid),
#                                          proj4string = vessel_density2@proj4string)

####### Interpolator function needed for INLAbru -- change elevation to the desired dataset
# f.elev <- function(x,y) {
#   # turn coordinates into SpatialPoints object:
#   spp <- SpatialPoints(data.frame(x=x,y=y)) 
#   # attach the appropriate coordinate reference system (CRS)
#   proj4string(spp) <- CRS(proj4string(elev))
#   # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
#   v <- over(spp, elev) 
#   v[is.na(v)] <- 0 # NAs are a problem! Remove them
#   return(v$elevation)
# } 

# Look at locations of BCCSN observations
BCCSN = BCCSNAlbers[BCCSNAlbers$Source == 'BCCSN',]
setwd("~/ownCloud/Whale Project/Final Project Scripts")
COAST_simp = readRDS('COAST_simp.rds')
plot(COAST_simp)
points(x = BCCSN$lon, y = BCCSN$lat,col='blue',cex=0.2,pch=16)
points(x = BCCSN$lon, y = BCCSN$lat,col='green',cex=.05,pch=16)

# Look at locations of BCCSN Whale watch observations from 2009 - 2011
library(readxl)
BCCSN2 <- read_excel("~/ownCloud/Whale Project/Summer-ClippedandChecked.xls")
unique(BCCSN2$OBSERVER_1)
# get all Whale watch boats
WW_BCCSN = BCCSN2[BCCSN2$OBSERVER_1 %in% c('Eco-Tourism -','Eco-Tou ism -','Other'),]
table(WW_BCCSN$YR) # 2001 - 2011 only - do we need to request more? does this provide additional spatial coverage?
WW_BCCSN = WW_BCCSN[!(WW_BCCSN$ORGANIZATI %in% c('No organizat','DFO - Scienc')),] # remove scientists
unique(WW_BCCSN$ORGANIZATI)

WW_BCCSN = SpatialPointsDataFrame(coords = cbind(WW_BCCSN$LONDEC, WW_BCCSN$LATDEC), data = WW_BCCSN, proj4string = CRS("+init=epsg:3857"))
#saveRDS(WW_BCCSN, file = 'WW_BCCSN_labelled_2001-2011.rds')

# Download weather data for the time period 2009 - 2016 #
#install.packages('weathercan')
#library(weathercan)
# find weather stations within AOI
#ids = stations_search(coords = cbind(48.43,-123.37), dist = 100, interval = "day")$station_id
# download weather from these stations
#historic_weather <- weather_dl(station_ids = as.numeric(as.character(ids[1])), stn = as.data.frame(weathercan::stations), start = "2009-02-01", end = "2009-02-20", interval = 'day', verbose = T)
