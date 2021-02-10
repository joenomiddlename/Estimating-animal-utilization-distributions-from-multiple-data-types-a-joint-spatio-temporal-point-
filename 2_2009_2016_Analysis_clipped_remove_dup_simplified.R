## Steps 
# 1 Build med res mesh for stationary models
# 2 Build high res mesh for nonstationary barrier model
# 3 Map covariates onto meshes
# 4 Fit models with AND without sampling intensity layers - use DIC for model selection and use med res mesh
# 5 Fit Barrier model to 'best' model using high-res mesh and compare (with high res mesh) to stationary model with LOOCV.

# Mac version 
setwd("~/Documents/Whale Project")

# PC version
setwd("~/ownCloud/Whale Project")

library(dplyr)
library(lubridate)
library(maptools)
library(rgeos)
library(rgdal)
library(INLA)
library(inlabru)
library(ggmap)
library(plyr)

# read in 
setwd("~/ownCloud/Whale Project/Final Project Scripts")
load('final_2009_2016_V3.RData')

# load the COAST shapefiles
# filled in and clipped shapefile
COAST_simp = readRDS('COAST_mesh.rds')
COAST_plotting = readRDS('COAST_simp.rds')
COAST_simp = spTransform(COAST_simp, CRS(proj4string(COAST)))

BG$TIME2 = strptime(as.character(BG$TIME_PDT), format='%H:%M:%S')
day(BG$TIME2) = as.numeric(BG$DAY)
month(BG$TIME2) = as.numeric(BG$MONTH)
year(BG$TIME2) = as.numeric(BG$YEAR)
head(BG$TIME2)


# Early corrections and strong assumptions
# Strong assumption 1 - For Brian's data - take the central location during each follow
 BG_DF = SpatialPointsDataFrame(coords = cbind(BG$Alb_East, BG$Alb_North),
                                data = BG, proj4string = CRS(proj4string(COAST)))
 BG_DF2 = BG_DF
 BG_DF2$DET_ID = factor(BG_DF2$DET_ID)
#  center_id = numeric()
#  count = 1
#  for(i in unique(BG_DF2$DET_ID))
#  {
#    ids = which(BG_DF2$DET_ID==i)
#    num_follow = length(ids)
#    center_id[count] = median(ids)
#    count = count+1
#  }
#  BG_subset = BG_DF2[center_id,]
# # check that each follow is only of one pod
#  by(BG_DF2$PODS, BG_DF2$DET_ID, FUN = function(x){length(unique(x))}) 
#  BG_subset = BG_DF2[center_id,]
 
 # Change to first sighting for consistency with WW data. Also reduces risk of sightings being
 # found in areas with predicted zero search effort.
 first_valid_id = numeric()
 count = 1
 test = 0
 for(i in unique(BG_DF2$DET_ID))
 {
   test=0
   count2 = 1
   while(test == 0)
   {
     ids = which(BG_DF2$DET_ID==i)
     #num_follow = length(ids)
     first_valid_id[count] = ids[count2]
     # test if point is in AOI
     if(gContains(COAST_simp, BG_DF2[ids[count2],])){
       test=1
       }
     count2 = count2 + 1
   }

   count = count+1
 }
 
 
 BG_subset = BG_DF2[first_valid_id,]
 
 BG_subset$Pod = revalue(BG_subset$PODS, c("J01"="J", "J01 K01"="JK", 
                                     'J01 K01 L01' = 'JKL', 'J01 L01' = 'JL', 
                                     'K01' = 'K', 'K01 L01' = 'KL','L01' = 'L'))

########### IMPORTANT DATA SUBSETTING POINT ############

## REMOVE ALL REPEATED SIGHTINGS OF THE SAME POD IN THE WW DATA ##
## THIS IS SINCE IT IS ALMOST GUARANTEED THAT THE BOATS REMAINED WITH WHALE ##
## SINCE WE REMOVED ALL REPEATED SIGHTINGS BY BRIAN, WE DO SO WITH WW ##
# Find indices that we want 
ind_dup = duplicated(WW@data[,c('Year','Month','Day','Pod')], fromLast = FALSE) # take the first such sighting

# add TIME2 variable to WW data
 ### Verfied time conversion
 WW$DateTimes = as.POSIXct(3600*24 * (WW$Matdate - min(WW$Matdate)), 
                           tz = 'UTC', origin = '2009-05-04 14:15', 
                           format = "%Y-%m-%d %H:%M")
 WW$TIME2 = WW$DateTimes
 WW$Date = as.Date(WW$DateTimes)
 WW_sightings$DateTimes = as.POSIXct(3600*24 * (WW_sightings$Matdate - min(WW_sightings$Matdate)), 
                           tz = 'UTC', origin = '2009-05-04 14:15', 
                           format = "%Y-%m-%d %H:%M")
 WW_sightings$TIME2 = WW_sightings$DateTimes
 WW_sightings$Date = as.Date(WW_sightings$DateTimes)
 
 SRKW_days = unique(WW@data$Date[WW@data$Pod == 'SRKW'])
 
# merge BG with WW data
BG_subset = BG_subset[,c('YEAR','MONTH','DAY','Pod','TIME2')]
WW_sightings = WW_sightings[!ind_dup,c('YEAR','MONTH','DAY','Pod','TIME2')]
#colnames(WW_sightings@data)=c("YEAR","MONTH","PODS" )

total_sightings = rbind.SpatialPointsDataFrame(BG_subset, WW_sightings)

# coordinates for the 'area of interest'
AOI = Polygon(cbind(c(970,1270,1270,970)*1000,c(350,350,500,500)*1000))
AOIs = Polygons(list(AOI), 'AOI') 
#COAST=readOGR("~/ownCloud/Whale Project/Data for Joe/Modified shapefiles/",layer="Clipped_filled2")

# Load the simplified/smoothed shapefile and the high res origninal used for plotting
#COAST= readOGR("~/ownCloud/Whale Project/Data for Joe/Modified shapefiles/",layer="Max_simp")
#COAST_plotting=readOGR("~/ownCloud/Whale Project/Data for Joe/Coast shapefile/",layer="coast_ALBERS")
# 
AOI = SpatialPolygons(list(AOIs), as.integer(1), proj4string = CRS(proj4string(COAST)))
# 
# COAST_AOI = gIntersection(gSimplify(COAST, tol=2000), AOI)
# plot(COAST_AOI)
# COAST_simp = gDifference(AOI, COAST_AOI)

# plot(COAST_simp)
# gIsValid(COAST_simp)
# 
#COAST_AOI_plotting = gIntersection(gSimplify(COAST_plotting, tol=1000), AOI)
#plot(COAST_AOI_plotting)
#COAST_plotting = gDifference(AOI, COAST_AOI_plotting)
# 
#plot(COAST_plotting)
# gIsValid(COAST_plotting)

######################

# All the above have been applied to the following RData file
# Points lying on land defined by INLA have been jittered to ensure they lie in water.

#######################

# load utility functions and plotting scripts
source('~/ownCloud/Whale Project/Final Project Scripts/plotting_scripts.R')
source('~/ownCloud/Whale Project/Final Project Scripts/utility_functions.R')

#load("~/ownCloud/Whale Project/Final Project Scripts/New_mesh_WS.RData")

#tmp = gContains(poly.barrier,total_sightings, byid = T)
#sum(tmp) # 0 point lies outside polygons... 

tmp = gContains(COAST_simp,total_sightings, byid = T)
sum(!tmp) # 19 points lie on land
tmp = gContains(COAST_plotting,total_sightings, byid = T)
sum(!tmp) # 21 lie on high res shapefile - not of concern as we have removed these polygons in INLA

# points(	BG[,c("Alb_East","Alb_North")],col=rgb(.8,0,0,.01),cex=1.2,pch=16)
# points(	BG[,c("Alb_East","Alb_North")],col='orangered',cex=.2,pch=16)
# points(WW$lon,WW$lat,  
#        pch=16, cex=1.2, col=rgb(1,1,0,.02))
# points(coordinates(WW)[,1],coordinates(WW)[,2],  
#        pch=16, cex=.2, col=rgb(1,.8,0,.2))

# Check all data lies within COAST_simp
check = gContains(COAST_simp,
                  SpatialPoints(coords = cbind(WW_sightings@coords[,1], WW_sightings@coords[,2]),
                                proj4string = CRS(proj4string(COAST))),
                  byid = T)
sum(check)
sum(!check)
check2 = gContains(COAST_simp,
                   SpatialPoints(coords = cbind(BG_subset@coords[,1], BG_subset@coords[,2]),
                                 proj4string = CRS(proj4string(COAST))),
                   byid = T)
sum(check2)
sum(!check2)
# 19 lie on land in WW data, 0 is out of AOI in BG - due to resolution issues. 
plot(COAST_simp)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col=rgb(.8,0,0,.01),cex=1.2,pch=16)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col='orangered',cex=.2,pch=16) #lies on land - move south
points(coordinates(WW_sightings)[!check,],  
       pch=16, cex=1.2, col='red')
# lies on land - move west

# move BG south
# gContains(COAST_simp,
#           SpatialPoints(coords = cbind(BG$Alb_East[!check2], BG$Alb_North[!check2]-1000),
#                         proj4string = CRS(proj4string(COAST))),
#           byid = T) # lies in water now
# plot(COAST_simp)
# points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,1000),col=rgb(.8,0,0,.01),cex=3,pch=16)
# BG$Alb_North[!check2] = BG$Alb_North[!check2] - 1000

# move WW west
gContains(COAST_simp,
          SpatialPoints(coords = cbind(WW_sightings@coords[!check,1]-2000, WW_sightings@coords[!check,2]),
                        proj4string = CRS(proj4string(COAST))),
          byid = T) # a few still lie in water now
check3 = gContains(COAST_simp,
                   SpatialPoints(coords = cbind(WW_sightings@coords[!check,1]-2000, WW_sightings@coords[!check,2]),
                                 proj4string = CRS(proj4string(COAST))),
                   byid = T)
plot(COAST_simp)
points(	WW_sightings@coords[!check,][check3,]-cbind(rep(2000,sum(check3)) ,rep(0,sum(check3))),col='red',cex=0.5,pch=16)
WW_sightings@coords[!check,1][check3] = WW_sightings@coords[!check,1][check3] - 2000
#WW_sightings@coords[!check,1][check3] = WW_sightings@coords[!check,1][check3] - 2000
# find the remaining 2 points 
ggplot() + gg(COAST_simp, alpha = 0.5) +
  gg(WW_sightings[which(!check),][which(!check3),][1,], colour = 'red') + # move northeast
  gg(WW_sightings[which(!check),][which(!check3),][2,], colour = 'orange') #+ # move northeast
  #gg(WW_sightings[which(!check),][which(!check3),][3,], colour = 'yellow') + # move east
  #gg(WW_sightings[which(!check),][which(!check3),][4,], colour = 'green') # move east
WW_sightings@coords[!check,1][!check3][1] = WW_sightings@coords[!check,1][!check3][1] + 3000
WW_sightings@coords[!check,2][!check3][1] = WW_sightings@coords[!check,2][!check3][1] + 2000
#WW_sightings$lon[!check][!check3][1] = WW_sightings$lon[!check][!check3][1] + 2000
WW_sightings@coords[!check,2][!check3][2] = WW_sightings@coords[!check,2][!check3][2] +2000
WW_sightings@coords[!check,1][!check3][2] = WW_sightings@coords[!check,1][!check3][2] +2000

# quick double check
check = gContains(COAST_simp,
                  SpatialPoints(coords = WW_sightings@coords,
                                proj4string = CRS(proj4string(COAST))),
                  byid = T)
sum(check)
sum(!check)

check = gContains(COAST_simp,
                  SpatialPoints(coords = WW@coords,
                                proj4string = CRS(proj4string(COAST))),
                  byid = T)
sum(check)
sum(!check)
# check2 = gContains(COAST_simp,
#                    SpatialPoints(coords = cbind(BG$Alb_East, BG$Alb_North),
#                                  proj4string = CRS(proj4string(COAST))),
#                    byid = T)
# sum(check2)
# sum(!check2)
# success!
ggplot() + gg(COAST_simp, alpha = 0.5) +
  gg(WW_sightings)

ggplot() + gg(COAST_simp, alpha = 0.5) +
  gg(WW[which(!check),])

#load('./Final Project Scripts/final_2009_2016_V3.RData')





# mesh = inla.mesh.2d(#loc = cbind(Dat_boat_red$Alb_East, Dat_boat_red$Alb_North),
#                     boundary = COAST_simp,
#                     #max.n = 1000,
#                     cutoff = 1000,
#                     offset = c(7000,40000),
#                     max.edge = c(7000, 30000),
#                     min.angle = 21,
#                     plot.delay = 1,
#                     crs = CRS(proj4string(COAST)))
#                     #cutoff = 500,
                    #max.edge=c(1000,5000))$n


######## NEW 

boundary <- list(
  as.inla.mesh.segment(COAST_simp),
  NULL)

## Build the mesh:
mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(7000, 25000),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=1000, ## Filter away adjacent points.
                     offset=c(7000, 25000)) ## Offset for extra boundaries, if needed.
# meshbuilder() shows this mesh is good

## Build the barrier mesh - not meshbuilder() suggests 10km is smallest range estimable by barrier:

mesh_barrier <- inla.mesh.2d(boundary=boundary,
                                     max.edge=c(3000, 25000),
                                     min.angle=c(21, 21),
                                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                                     cutoff=1000, ## Filter away adjacent points.
                                     offset=c(14000, 18000)) ## Offset for extra boundaries, if needed.

## Plot the mesh:
plot(mesh)
plot(mesh_barrier)
#####

mesh$n
mesh_barrier$n

# Due to different R-INLA versions, load the mesh's pre-made for use using the above code.
# This is since different versions create different mesh$n. Code will not run without.
mesh = readRDS('mesh.rds')
mesh_barrier = readRDS('mesh_barrier.rds')

# check how well the mesh boundary covers the COAST_simp shapefile
ggplot() + gg(mesh, colour = 'blue') #+ gg(COAST_simp, colour = 'red')
 ggplot() + gg(mesh, colour = 'blue') + gg(COAST_simp, colour = 'red')
 ggplot() + gg(mesh, colour = 'blue') #+ gg(COAST_simp, colour = 'red')
 ggplot() + gg(mesh_barrier, colour = 'blue') + gg(COAST_simp, colour = 'red')


gmap(spTransform(merge(BG_DF,WW_sightings), CRS("+proj=longlat")),zoom=7) + gg(spTransform(BG_DF, CRS("+proj=longlat")),alpha = 0.2) + gg(spTransform(WW_sightings, CRS("+proj=longlat")),alpha =0.2) + gg(spTransform(COAST_simp, CRS("+proj=longlat"))) #+ gg(inla.spTransform(mesh, CRS("+proj=longlat")))
gmap(spTransform(merge(BG_DF,WW_sightings), CRS("+proj=longlat")),zoom=7) +  gg(spTransform(COAST_simp, CRS("+proj=longlat"))) #+ gg(inla.spTransform(mesh, CRS("+proj=longlat")))
ggplot()+ gg(spTransform(BG_DF, CRS("+proj=longlat")),alpha = 0.2) + gg(spTransform(WW_sightings, CRS("+proj=longlat")),alpha =0.2) + gg(spTransform(COAST_simp, CRS("+proj=longlat"))) + gg(inla.spTransform(mesh, CRS("+proj=longlat")))
ggplot()+ gg(BG_DF,alpha = 0.2) + gg(WW_sightings,alpha =0.2) + gg(mesh)

#save.image('final_2009_2016_V4_2_removed_dup.RData')
#load("~/ownCloud/Whale Project/Final Project Scripts/final_2009_2016_V4_2_removed_dup.RData")

# create the barrier model from Haakon Bakka

#First we divide up the mesh accoring to our study area polygon.

tl = length(mesh_barrier$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh_barrier$loc[mesh_barrier$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri, proj4string = CRS(proj4string(COAST_simp)))

# - compute the triangle positions
#normal = over(COAST_simp, posTri, returnList=T)
polygon.triangles = inla.over_sp_mesh(COAST_simp, y = mesh_barrier, type = "centroid", ignore.CRS=T)
polygon.triangles = setdiff(1:tl, polygon.triangles)
# - checking which mesh triangles are inside the normal area
#normal = unlist(normal)
#barrier.triangles = setdiff(1:tl, normal)
plot(posTri[polygon.triangles,])
# poly.triangles needs to make up the barrier area

#mesh$crs = NULL
#poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)
poly.barrier = inla.barrier.polygon(mesh_barrier, polygon.triangles)
proj4string( poly.barrier ) = CRS(proj4string(COAST_simp))
ggplot() + gg(poly.barrier, alpha = 0.5)
# remember this polygon is for the BARRIER - do not use for samplers

#INLAbru does not work yet with LGCP barrier

# Early exploratory models - how much information is in the data?
# Best to aggregate over years for each month.

######### Model files  ############

# create the projector matrix
A_proj = inla.spde.make.A(mesh, coordinates(total_sightings))
A_proj_barrier = inla.spde.make.A(mesh_barrier, coordinates(total_sightings))
# spacetime below
#A_proj = inla.spde.make.A(mesh, coordinates(total_sightings), group = , ngroup = )

dmesh = inla.mesh.dual(mesh)
proj4string(dmesh) = CRS(proj4string( COAST_simp) )
plot(dmesh)

dmesh_barrier = inla.mesh.dual(mesh_barrier)
proj4string(dmesh_barrier) = CRS(proj4string( COAST_simp) )
plot(dmesh_barrier)

# Obtain the intersection between each polygon from the dual mesh with the COAST_ polygon.
library(rgeos)
sum(w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i,], COAST_simp))
    return(gArea(gIntersection(dmesh[i,], COAST_simp)))
  else return(0)
}))
sum(w_barrier <- sapply(1:length(dmesh_barrier), function(i) {
  if (gIntersects(dmesh_barrier[i,], COAST_simp))
    return(gArea(gIntersection(dmesh_barrier[i,], COAST_simp)))
  else return(0)
}))

# Check that the area of the intersection is close to the area of GB_polygons_simplified
sum(w) / gArea(COAST_simp)
sum(w_barrier) / gArea(COAST_simp)

## ----wsummary------------------------------------------------------------
# Check there are integration points with zero weight (sanity check)
table(w>0) # 357 lie on land
table(w_barrier>0) # 690 lie on land

### Bring in covariates
# We saw earlier that some of the satellite rasters had unreasonable values at some pixels
# Need to identify and remove outliers, whilst smoothing the raster.
# chlorophyll dataset has NAs but no outliers, temperature has both, 
 hist(as.numeric(as.matrix(chloro@data)))
 #library(MODIS)
 sum(is.na(as.matrix(SST_1@data))) # c('2006-04-16', '2007-04-16')
 sum(is.na(as.matrix(SST_2@data))) # c('2007-04-16', '2016-11-16') overlap 1 month
 sum(is.na(as.matrix(chloro@data))) # 2006-04-16 - 2016-11-16 70th month is 2012-01-16
 sum(is.na(as.matrix(chloro2@data))) # 2012-01-15 - 2016-11-15
 sum(is.na(as.matrix(topo@data))) # fine
 sum(is.na(as.matrix(vessel_density2@data))) # fine
# 
# # map the dates to the column names
 #tpos <- c('2006-04-01', '2016-11-30') 
 #tpos1 <- c('2006-04-01', '2007-04-15')
 #tpos2 <- c('2007-04-16', '2016-11-30')
 #tpos3 = c('2012-01-15','2016-11-30') 
 
 #library(xtractomatic)
 
#colnames(SST_1@data) = paste('date_',substr(gsub('-','',xtracto_3D("agsstamday", xpos, ypos, tpos1)$time), start = 1, stop  = 6), sep='')
#colnames(SST_2@data) = paste('date_',substr(gsub('-','',xtracto_3D("atsstamday", xpos, ypos, tpos2)$time), start = 1, stop  = 6), sep='')
#colnames(chloro@data) = paste('date_',substr(gsub('-','',xtracto_3D("mhchlamday", xpos, ypos, tpos)$time), start = 1, stop  = 6), sep='')
#colnames(chloro2@data) = paste('date_',substr(gsub('-','',xtracto_3D("erdVH3chlamday", xpos, ypos, tpos3)$time), start = 1, stop  = 6), sep='')
# 
# # save time with pre-loaded names
names <- readRDS("~/ownCloud/Whale Project/names.rds")
names(SST_1) = names$SST1
names(SST_2) = names$SST2
names(chloro) = names$chloro
names(chloro2) = names$chloro2
# 
# library(raster)
 setwd("~/ownCloud/Whale Project/test_shapefile/Smoothed")
 rlist2=list.files(getwd(), pattern="tif$", full.names=FALSE)
 vessel_density = readGDAL(rlist2[1])
 plot(vessel_density)
 crs(vessel_density) = CRS("+init=epsg:3857")
 sum(is.na(as.matrix(vessel_density@data))) # fine
# vessel_density_medres = readGDAL(rlist2[1])
# plot(vessel_density_medres)
# hist(vessel_density_medres@data$band1)
# 
 #### Smooth and clean the remote data - outliers exist in SST2
dev.off()
par(mfrow = c(1,1))
 hist(unlist(chloro@data))
 hist(unlist(chloro2@data))
 hist(unlist(SST_1@data))
 hist(unlist(SST_2@data))
 # remove SST2 values less than 4 and greater than 25
 hist(SST_2@data[,2])
 
 SST_2_3 = SST_2
 # how many outliers?
 sum(SST_2@data < 4| SST_2@data > 25, na.rm = T)/sum(!is.na(SST_2@data)) # 0.36% implausible values
 SST_2@data[SST_2@data < 4 | SST_2@data > 25] = NA
 
 # Now we subset by the months we need - 5,6,7,8,9,10
 colnames(SST_1@data)
 SST_1_TOI = SST_1
 SST_1_TOI@data = SST_1_TOI@data[,substr(colnames(SST_1_TOI@data),10,11) %in% c('05','06','07','08','09','10')]
 colnames(SST_1_TOI@data)
 
 colnames(SST_2@data)
 SST_2_TOI = SST_2
 SST_2_TOI@data = SST_2_TOI@data[,substr(colnames(SST_2_TOI@data),10,11) %in% c('05','06','07','08','09','10')]
 colnames(SST_2_TOI@data) # missing 200806, 201210 - impute as average of previous and next month
 SST_2_TOI@data[,'date_200806'] = apply(cbind(SST_2_TOI@data[,'date_200805'], SST_2_TOI@data[,'date_200807']), 1, mean, na.rm=T)
 SST_2_TOI@data[,'date_201210'] = apply(cbind(SST_2@data[,'date_201209'], SST_2@data[,'date_201211']), 1, mean, na.rm=T)
 SST_2_TOI@data = SST_2_TOI@data[,sort(colnames(SST_2_TOI@data))] # reorder in time order
 
 colnames(chloro@data)
 chloro_TOI = chloro
 chloro_TOI@data = chloro_TOI@data[,substr(colnames(chloro_TOI@data),10,11) %in% c('05','06','07','08','09','10')]
 colnames(chloro_TOI@data) # complete data
 chloro_TOI@data = chloro_TOI@data[,1:36]
 
 colnames(chloro2@data)
 chloro2_TOI = chloro2
 chloro2_TOI@data = chloro2_TOI@data[,substr(colnames(chloro2_TOI@data),10,11) %in% c('05','06','07','08','09','10')]
 colnames(chloro2_TOI@data) # complete data
 
 # how many missings are there each month - look for a complete missing month?
 colSums(is.na(chloro_TOI@data))/dim(chloro_TOI@data)[1] # stable
 colSums(is.na(chloro2_TOI@data))/dim(chloro2_TOI@data)[1] # stable
 colSums(is.na(SST_1_TOI@data))/dim(SST_1_TOI@data)[1] # stable
 colSums(is.na(SST_2_TOI@data))/dim(SST_2_TOI@data)[1] # stable
 
# map the variables onto the observation locations
 
 # I got the projections wrong - change CRS into Long Lat
 chloro_TOI@proj4string = CRS("+init=epsg:4326")
 chloro2_TOI@proj4string = CRS("+init=epsg:4326")
 SST_1_TOI@proj4string = CRS("+init=epsg:4326")
 SST_2_TOI@proj4string = CRS("+init=epsg:4326")
 topo@proj4string = CRS("+init=epsg:4326")

############# Summary - we do not have enough information in the data to estimate the field
# every month, every year. Average over the years to combine information.
xtabs(~ YEAR + MONTH, data = total_sightings@data )
xtabs(~ MONTH, data = total_sightings@data )

 # remove Brian's 1 observation in Nov and 1 in April to put on same time as WW
total_sightings = rbind(WW_sightings, BG_subset)
total_sightings = total_sightings[total_sightings$MONTH %in% c(5,6,7,8,9,10),]

# We can reuse the dmesh projections and simply average over the years for each pixel.
# Need to project monthly averages onto observation locations however
# create monthly average grids

# chloro - same resolution from two satellites so merge. identical(chloro_TOI@coords, chloro2_TOI@coords) == TRUE
chloro_TOI_month = cbind(chloro_TOI@data, chloro2_TOI@data)
chloro_TOI_month = chloro_TOI_month[,-c(1:18)]
chloro_month_avg = list()
chloro_month_avg$May = rowMeans(chloro_TOI_month[,seq(from = 1, to = 48, by = 6)], na.rm = T)
chloro_month_avg$June = rowMeans(chloro_TOI_month[,seq(from = 2, to = 48, by = 6)], na.rm = T)
chloro_month_avg$July = rowMeans(chloro_TOI_month[,seq(from = 3, to = 48, by = 6)], na.rm = T)
chloro_month_avg$Aug = rowMeans(chloro_TOI_month[,seq(from = 4, to = 48, by = 6)], na.rm = T)
chloro_month_avg$Sep = rowMeans(chloro_TOI_month[,seq(from = 5, to = 48, by = 6)], na.rm = T)
chloro_month_avg$Oct = rowMeans(chloro_TOI_month[,seq(from = 6, to = 48, by = 6)], na.rm = T)
chloro_monthly = chloro_TOI
chloro_monthly@data = as.data.frame(chloro_month_avg)
rm(chloro_TOI_month)
rm(chloro_month_avg)

ggplot() + gg(chloro_monthly['May'])

# Make some required files for mapping covariates
# In particular, we obtain the distances between each dual mesh pixel
dmesh_dist = gDistance(dmesh,byid = T)
dmesh_barrier_dist = gDistance(dmesh_barrier,byid = T)

chloro_obs_month = f.interp(sp_points = total_sightings, 
                            sp_grid = chloro_monthly)
chloro_dmesh_month = f.interp.dmesh(dmesh, chloro_monthly, dmesh_dist)
chloro_dmesh_barrier_month = f.interp.dmesh(dmesh_barrier, chloro_monthly, dmesh_barrier_dist)


# A side by side plot of weighted mean vs median is a good way to check for outliers
par(mfrow=c(2,1))
hist(as.matrix(chloro_dmesh_month$weighted_mean@data[,]), main = 'Weighted Mean')
hist(as.matrix(chloro_dmesh_month$median@data[,]),main = 'Median')
par(mfrow=c(1,1))

par(mfrow=c(2,1))
hist(as.matrix(chloro_dmesh_barrier_month$weighted_mean@data[,]), main = 'Weighted Mean')
hist(as.matrix(chloro_dmesh_barrier_month$median@data[,]),main = 'Median')
par(mfrow=c(1,1))

# SST - resolutions different. SST_2 is from 2009 onwards
SST_TOI_month = SST_2_TOI@data
SST_TOI_month = SST_TOI_month[,-c(1:12)]
SST_month_avg = list()
SST_month_avg$May = rowMeans(SST_TOI_month[,seq(from = 1, to = 48, by = 6)], na.rm = T)
SST_month_avg$June = rowMeans(SST_TOI_month[,seq(from = 2, to = 48, by = 6)], na.rm = T)
SST_month_avg$July = rowMeans(SST_TOI_month[,seq(from = 3, to = 48, by = 6)], na.rm = T)
SST_month_avg$Aug = rowMeans(SST_TOI_month[,seq(from = 4, to = 48, by = 6)], na.rm = T)
SST_month_avg$Sep = rowMeans(SST_TOI_month[,seq(from = 5, to = 48, by = 6)], na.rm = T)
SST_month_avg$Oct = rowMeans(SST_TOI_month[,seq(from = 6, to = 48, by = 6)], na.rm = T)
SST_monthly = SST_2_TOI
SST_monthly@data = as.data.frame(SST_month_avg)
rm(SST_TOI_month)
rm(SST_month_avg)

SST_obs_month = f.interp(sp_points = total_sightings, 
                            sp_grid = SST_monthly)
SST_dmesh_month = f.interp.dmesh(dmesh, SST_monthly, dmesh_dist)
SST_dmesh_barrier_month = f.interp.dmesh(dmesh_barrier, SST_monthly, dmesh_barrier_dist)


# A side by side plot of weighted mean vs median is a good way to check for outliers
par(mfrow=c(2,1))
hist(as.matrix(SST_dmesh_month$weighted_mean@data[,]), main = 'Weighted Mean')
hist(as.matrix(SST_dmesh_month$median@data[,]),main = 'Median')
par(mfrow=c(1,1))
par(mfrow=c(2,1))
hist(as.matrix(SST_dmesh_barrier_month$weighted_mean@data[,]), main = 'Weighted Mean')
hist(as.matrix(SST_dmesh_barrier_month$median@data[,]),main = 'Median')
par(mfrow=c(1,1))

# topo
topo_obs = f.interp(sp_points = total_sightings, 
                    sp_grid = topo)
topo_dmesh = f.interp.dmesh(dmesh, topo, dmesh_dist)
topo_dmesh_barrier = f.interp.dmesh(dmesh_barrier, topo, dmesh_barrier_dist)


# Obtain vessel density values at observation locations
vessel_obs = f.interp(sp_points = total_sightings, 
                      sp_grid = vessel_density)
vessel_dmesh = f.interp.dmesh(dmesh, vessel_density, dmesh_dist)
vessel_dmesh_barrier = f.interp.dmesh(dmesh_barrier, vessel_density, dmesh_barrier_dist)


######################################################
# MAP THE POINTS OONTO THE CORRECT MONTH'S COVARS
######################################################
# For INLA we need to convert Month to 1-6
total_sightings@data$MONTH_INLA = total_sightings@data$MONTH - 4

# USE MONTH_INLA TO INDEX TO COVARIATES
colnames(chloro_obs_month) = c('1','2','3','4','5','6')

# create mappings to the covariates
covar_map = sapply(total_sightings@data$MONTH_INLA,
                   FUN = function(x){which(colnames(chloro_obs_month) == x)},
                   simplify = T)
covar_map = unlist(covar_map)

total_sightings@data$chloro = chloro_obs_month[cbind(seq(1,2433),as.numeric(covar_map))]
total_sightings@data$SST = SST_obs_month[cbind(seq(1,2433),as.numeric(covar_map))]


# map to dual mesh
covariates_dmesh = data.frame(chloro = as.numeric(as.matrix(chloro_dmesh_month$weighted_mean@data)),
  SST = as.numeric(as.matrix(SST_dmesh_month$weighted_mean@data)),
  topo = rep(topo_dmesh$weighted_mean@data$V1,times = dim(SST_dmesh_month$weighted_mean)[2]),
  vessel = rep(vessel_dmesh$weighted_mean@data$band1,times = dim(SST_dmesh_month$weighted_mean)[2]),
  MONTH = rep(5:10, each = mesh$n),
  stringsAsFactors = F
)
# map to barrier dual mesh
covariates_dmesh_barrier = data.frame(chloro = as.numeric(as.matrix(chloro_dmesh_barrier_month$weighted_mean@data)),
                              SST = as.numeric(as.matrix(SST_dmesh_barrier_month$weighted_mean@data)),
                              topo = rep(topo_dmesh_barrier$weighted_mean@data$V1,times = dim(SST_dmesh_barrier_month$weighted_mean)[2]),
                              vessel = rep(vessel_dmesh_barrier$weighted_mean@data$band1,times = dim(SST_dmesh_barrier_month$weighted_mean)[2]),
                              MONTH = rep(5:10, each = mesh_barrier$n),
                              stringsAsFactors = F
)

#save.image('Workspace for naive models_remove_dup.RData')
# All of the above has been saved
#load('Workspace for naive models_remove_dup.RData')

no_T = 6

# covariates for the dual mesh for the correct years and months
covariates_pp = covariates_dmesh#[covariates_dmesh$YEAR >= 2009 & covariates_dmesh$MONTH %in% c(5,6,7,8,9,10),]
covariates_pp$MONTH_INLA = covariates_pp$MONTH - 4
covariates_pp_barrier = covariates_dmesh_barrier#[covariates_dmesh$YEAR >= 2009 & covariates_dmesh$MONTH %in% c(5,6,7,8,9,10),]
covariates_pp_barrier$MONTH_INLA = covariates_pp_barrier$MONTH - 4
# spacetime
#covariates_pp = data.frame(R_lag = rep(0, times = (mesh$n*no_T)))
total_sightings@data$topo = as.numeric(topo_obs$V1)
total_sightings@data$vessel = vessel_obs$band1
## Make covariates at site locations
covariates_site = total_sightings@data

# remove unused columns to enable binding to work
#covariates_site = covariates_site[,-c(1,3,4,5)]
#covariates_pp = covariates_pp[,-c(5)]
# Keep MONTH, chloro, SST, vessel, MONTH_INLA and topo

# Remember chloro, vessel and depth all should be log transformed
covariates_pp$chloro = log(covariates_pp$chloro)
covariates_pp_barrier$chloro = log(covariates_pp_barrier$chloro)
covariates_site$chloro = log(covariates_site$chloro)
covariates_pp$vessel = log(covariates_pp$vessel+1)
covariates_pp_barrier$vessel = log(covariates_pp_barrier$vessel+1)
covariates_site$vessel = log(covariates_site$vessel+1)

# All variables are on similar scales so no need to rescale. This also keeps the nice interpretation of SST
  
# Following the advice of YUAN et al, we center each spatiotemporal covariate in 3 ways.
# We center within-month (i.e. average SST across space for each month, and center the spatial SST within each month)
# Finally we compute the interaction between the temporal pattern and the spatial average 
SST_overallmean = mean(covariates_pp$SST)
chloro_overallmean = mean(covariates_pp$chloro)

SST_monthlymeans = as.numeric(by(covariates_pp$SST, as.factor(covariates_pp$MONTH), FUN = mean, na.rm=T ))
chloro_monthlymeans = as.numeric(by(covariates_pp$chloro, as.factor(covariates_pp$MONTH), FUN = mean, na.rm=T ))

# find the dmesh average across time
covariates_pp$meshind = rep(1:mesh$n, times = dim(covariates_pp)[1]/mesh$n)
covariates_pp_barrier$meshind = rep(1:mesh_barrier$n, times = dim(covariates_pp_barrier)[1]/mesh_barrier$n)

SST_spatialmeans = as.numeric(by(covariates_pp$SST, as.factor(covariates_pp$meshind), FUN = mean, na.rm=T ))
chloro_spatialmeans = as.numeric(by(covariates_pp$chloro, as.factor(covariates_pp$meshind), FUN = mean, na.rm=T ))

SST_spatialmeans_barrier = as.numeric(by(covariates_pp_barrier$SST, as.factor(covariates_pp_barrier$meshind), FUN = mean, na.rm=T ))
chloro_spatialmeans_barrier = as.numeric(by(covariates_pp_barrier$chloro, as.factor(covariates_pp_barrier$meshind), FUN = mean, na.rm=T ))

# find the pixel corresponding to each observation
which_dmesh_pixel = gWithin(total_sightings,dmesh, byid = T, returnDense = F)
which_dmesh_pixel = unlist(which_dmesh_pixel)
which_dmesh_pixel = as.numeric(which_dmesh_pixel)

# create spatiotemporal covariates
covariates_pp$SSTmonthavg = rep(SST_monthlymeans, each = mesh$n) - SST_overallmean # repeat each month
covariates_pp$chloromonthavg = rep(chloro_monthlymeans, each = mesh$n) - chloro_overallmean
covariates_pp_barrier$SSTmonthavg = rep(SST_monthlymeans, each = mesh_barrier$n) - SST_overallmean # repeat each month
covariates_pp_barrier$chloromonthavg = rep(chloro_monthlymeans, each = mesh_barrier$n) - chloro_overallmean

covariates_pp$SSTminusmonth = covariates_pp$SST - rep(SST_monthlymeans, each = mesh$n) # how warm is the location compared to monthly average
covariates_pp$chlorominusmonth = covariates_pp$chloro - rep(chloro_monthlymeans, each = mesh$n)
covariates_pp_barrier$SSTminusmonth = covariates_pp_barrier$SST - rep(SST_monthlymeans, each = mesh_barrier$n) # how warm is the location compared to monthly average
covariates_pp_barrier$chlorominusmonth = covariates_pp_barrier$chloro - rep(chloro_monthlymeans, each = mesh_barrier$n)

covariates_pp$SSTspatialavg = rep(SST_spatialmeans, times = dim(covariates_pp)[1]/mesh$n) - SST_overallmean # is the region typically? warm 
covariates_pp$chlorospatialavg = rep(chloro_spatialmeans, times = dim(covariates_pp)[1]/mesh$n) - chloro_overallmean
covariates_pp_barrier$SSTspatialavg = rep(SST_spatialmeans_barrier, times = dim(covariates_pp_barrier)[1]/mesh_barrier$n) - SST_overallmean # is the region typically? warm 
covariates_pp_barrier$chlorospatialavg = rep(chloro_spatialmeans_barrier, times = dim(covariates_pp_barrier)[1]/mesh_barrier$n) - chloro_overallmean

covariates_pp$SSTminusspacetime = covariates_pp$SST - SST_overallmean - covariates_pp$SSTmonthavg - covariates_pp$SSTspatialavg
covariates_pp$chlorominusspacetime = covariates_pp$chloro - chloro_overallmean - covariates_pp$chloromonthavg - covariates_pp$chlorospatialavg
covariates_pp_barrier$SSTminusspacetime = covariates_pp_barrier$SST - SST_overallmean - covariates_pp_barrier$SSTmonthavg - covariates_pp_barrier$SSTspatialavg
covariates_pp_barrier$chlorominusspacetime = covariates_pp_barrier$chloro - chloro_overallmean - covariates_pp_barrier$chloromonthavg - covariates_pp_barrier$chlorospatialavg

# map the spatiotemporal averages onto the points
covariates_site$SSTmonthavg = 0
covariates_site$chloromonthavg = 0

count = 1
for(i in unique(covariates_pp$MONTH))
{
  covariates_site$SSTmonthavg[covariates_site$MONTH == i] = SST_monthlymeans[count]
  covariates_site$chloromonthavg[covariates_site$MONTH == i] = chloro_monthlymeans[count]
  count = count + 1
}

covariates_site$SSTminusmonth = covariates_site$SST - covariates_site$SSTmonthavg # how warm is the location compared to monthly average
covariates_site$chlorominusmonth = covariates_site$chloro - covariates_site$chloromonthavg
covariates_site$SSTmonthavg = covariates_site$SSTmonthavg - SST_overallmean
covariates_site$chloromonthavg = covariates_site$chloromonthavg - chloro_overallmean

covariates_site$SSTspatialavg = SST_spatialmeans[which_dmesh_pixel] - SST_overallmean # is the region typically? warm 
covariates_site$chlorospatialavg = chloro_spatialmeans[which_dmesh_pixel] - chloro_overallmean

covariates_site$SSTminusspacetime = covariates_site$SST - SST_overallmean - covariates_site$SSTmonthavg - covariates_site$SSTspatialavg
covariates_site$chlorominusspacetime = covariates_site$chloro - chloro_overallmean - covariates_site$chloromonthavg - covariates_site$chlorospatialavg

# create a depth variable and take the log to reduce the skew.
# Note we take the absolute value without worry - the topo's above 51 lie on land and thus are not in AOI
covariates_pp$depth = log(abs(covariates_pp$topo-51)+1)
covariates_pp_barrier$depth = log(abs(covariates_pp_barrier$topo-51)+1)
covariates_site$depth = log(abs(covariates_site$topo-51)+1)

# create a reference mean depth variable for comparison
mean_depth_dmesh = mean(log(abs(covariates_pp$topo-51)+1))

covariates_pp$depth = covariates_pp$depth - mean_depth_dmesh
covariates_pp_barrier$depth = covariates_pp_barrier$depth - mean_depth_dmesh
covariates_site$depth = covariates_site$depth - mean_depth_dmesh

# Make some nice plots to demonstrate the barrier vs regular model's effect on the correlation near land
# Also becomes useful to obtain the eigenvectors of the precision matrix for evaluating confounding

# To make it comparable - need to make stationary model on same mesh
mesh_barrier_transformed = inla.spTransform(mesh_barrier, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
mesh_transformed = inla.spTransform(mesh, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
COAST_transformed = spTransform(COAST_plotting,CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))

barrier.model = inla.barrier.pcmatern(mesh_barrier_transformed, barrier.triangles = polygon.triangles, 
                                       prior.range = c(15, 0.01), 
                                       prior.sigma = c(3, 0.1), 
                                       range.fraction = 0.2) # remember we're not in lon/lat, 3e5 is domain size
 # reasonable lower range (median) is half of domain prior.range = c(15e4, .5)
 # reasonable upper scale is sd of 4

simple.model.highres = inla.spde2.pcmatern(mesh_barrier_transformed, 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1)) 

# Haakon Bakka's code
# define a function defining how to plot any spatial fields
colsc_corr <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = c(0.1,1))
}
plot.field = function(field, poly, mesh, pixels, corr = TRUE,...){
  stopifnot(length(field) == mesh$n)
  
  # - choose plotting region to be the same as the study area polygon
  #pixels_plotting = pixels(mesh, mask = poly, nx = 300, ny = 300)
  proj = inla.spde.make.A(mesh, loc = coordinates(pixels) ) 
  #proj = inla.mesh.projector(mesh, xlim = xlim, 
  #                           ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = as.numeric(proj %*% field)
  
  field_df = SpatialPixelsDataFrame(pixels, 
                                    data = data.frame(field = field.proj))
  # - Do the projection
  #image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
  #           xlim = xlim, ylim = ylim, ...)
  if(corr == TRUE)
  {
    plot = ggplot() + gg(field_df) + gg(poly) + colsc_corr()
  }
  if(corr != TRUE)
  {
    plot = ggplot() + gg(field_df) + gg(poly) + colsc(field.proj)
  }
  return(plot)
}



# Compute the correlation between the field at any point and this reference location
local.find.correlation = function(Q, location, mesh) {
  sd = sqrt(diag(inla.qinv(Q)))
  # - the marginal standard deviations
  
  A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
  # - create a fake A matrix, to extract the closest mesh node index
  id.node = which.max(A.tmp[1, ])
  # - index of the closest node
  
  print(paste('The location used was c(', 
              round(mesh$loc[id.node, 1], 4), ', ', 
              round(mesh$loc[id.node, 2], 4), ')' ))
  # - location of the closest node
  # - should be close to the location input
  # - sometimes used to plot a black dot
  
  ## Solve a matrix system to find the column of the covariance matrix
  Inode = rep(0, dim(Q)[1]); Inode[id.node] = 1
  covar.column = solve(Q, Inode)
  corr = drop(matrix(covar.column)) / (sd*sd[id.node])
  return(corr)
}

# plot barrier correlation
theta_max = c(0, 4.13) # range and SD
Q_barrier = inla.rgeneric.q(barrier.model, "Q", theta = c(0, 4.13)) # this is the range found in spatial only model
corr_barrier = local.find.correlation(Q_barrier, loc = c(1231.797,401.001), mesh_barrier_transformed)
pixels_plotting =  pixels(mesh_barrier_transformed, mask = COAST_transformed, nx = 300, ny = 300)
corr_plot_barrier = plot.field(corr_barrier,COAST_transformed, mesh_barrier_transformed, pixels_plotting)
corr_plot_barrier = corr_plot_barrier + gg(SpatialPoints(coords = cbind(1231.797,401.001), 
                                     proj4string = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))) +
  ggtitle('barrier model correlation field')
corr_plot_barrier

# plot stationary correlation
Q = inla.spde2.precision(simple.model.highres, theta = c(4.13,0)) # this is the range found in spatial only model
corr = local.find.correlation(Q, loc = c(1231.797,401.001), mesh_barrier_transformed)
pixels_plotting =  pixels(mesh_barrier_transformed, mask = COAST_transformed, nx = 300, ny = 300)
corr_plot_stationary = plot.field(corr,COAST_transformed, mesh_barrier_transformed, pixels_plotting)
corr_plot_stationary = corr_plot_stationary + gg(SpatialPoints(coords = cbind(1231.797,401.001), 
                                        proj4string = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))) +
  ggtitle('stationary model correlation field')
#plot(poly.barrier, add=T, col='grey', main = 'Barrier Model')
multiplot(corr_plot_barrier, corr_plot_stationary)

# Repeat but for an area not blocked by islands
corr_barrier2 = local.find.correlation(Q_barrier, loc = c(1010,400), mesh_barrier_transformed)
corr_plot_barrier2 = plot.field(corr_barrier2,COAST_transformed, mesh_barrier_transformed, pixels_plotting)
corr_plot_barrier2 = corr_plot_barrier2 + gg(SpatialPoints(coords = cbind(1010,400), 
                                                         proj4string = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))) +
  ggtitle('barrier model correlation field')

corr2 = local.find.correlation(Q, loc = c(1010,400), mesh_barrier_transformed)
corr_plot_stationary2 = plot.field(corr2,COAST_transformed, mesh_barrier_transformed, pixels_plotting)
corr_plot_stationary2 = corr_plot_stationary2 + gg(SpatialPoints(coords = cbind(1010,400), 
                                                               proj4string = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))) +
  ggtitle('stationary model correlation field')
#plot(poly.barrier, add=T, col='grey', main = 'Barrier Model')
multiplot(corr_plot_barrier2, corr_plot_stationary2)


# Now add the sampling intensity layers to the models above...

# Note that the number of boats were based on 2011 numbers. 
# According to Soundwatch reports (2016 in particular), 9am - 6pm is common search time.
# Canadian numbers remained stable 2009-2017 
# US numbers increased dramatically

# 2009 - Canadian = 52, US = 21,
# 2010 - Canadian = 52, US = 24 ,
# 2011 - Canadian = 54, US = 22,
# 2012 - Canadian = 54, US = 25,
# 2013 - Canadian = 54, US = 26,
# 2014 - Canadian = 57, US = 28,
# 2015 - Canadian = 57, US = 39
# 2016 - Canadian = 53, US = 49
vessel_numbers_per_year = data.frame(Can = c(52,52,54,54,54,57,57,53),
                                     US = c(21,24,22,25,26,28,39,49),
                                     Year = 2009:2016)
vessel_numbers_per_year_prop = vessel_numbers_per_year
vessel_numbers_per_year_prop[,c(1,2)] = t(apply(as.matrix(vessel_numbers_per_year[,c(1,2)]), 1, 
                                              FUN = function(x){x/as.matrix(vessel_numbers_per_year[vessel_numbers_per_year$Year == 2011,c(1,2)])}))
vessel_numbers_per_year_prop

vessel_numbers_per_year_diff = vessel_numbers_per_year
vessel_numbers_per_year_diff[,c(1,2)] = t(apply(as.matrix(vessel_numbers_per_year[,c(1,2)]), 1, 
                                                FUN = function(x){x-as.matrix(vessel_numbers_per_year[vessel_numbers_per_year$Year == 2011,c(1,2)])}))
vessel_numbers_per_year_diff
# Additional 130,000 visitors to Lime Kiln State Park from 2009 - 2016,
# 220,000 visitors 2009 - 347,000 visitors 2016
# This adds additional search effort to the region - justifying sensitivity analysis
# We also ignore contributions from kayaks and land-based observers ('keeners')

# We need the total number of days per year/month of acceptible weather conditions.

# Average number of whale watch boats with whale by month (Soundwatch 2011 - 2017)
# 
# Will use this to estimate the total proportion of vessels on water.
monthly_effort = t(cbind(c(NA,6,7,8,5),
                            c(4.75,7,7.25,4,4.75),
                            c(5,7,9,8,8),
                            c(5,8,9,8,6),
                            c(6,7,9,8,5),
                            c(4,7,9,6,7),
                            c(5,6,6,6,6) ))
colnames(monthly_effort) = c('May','Jun','July','Aug','Sep')
rownames(monthly_effort) = as.character(2017:2011)
# These are the reported average number of WW boats with whales per month
# Average these across the years (median as outliers present) to estimate effort
monthly_effort = apply(monthly_effort,2, median, na.rm = T)
monthly_effort = monthly_effort / max(monthly_effort)
dev.off()
par(mfrow=c(1,1))
plot(x = 1:5, y = monthly_effort)
lines(x = 1:5, y = monthly_effort)

# Extrapolate October's effort to equal May's 
monthly_effort[6] = monthly_effort[1]
names(monthly_effort) = c('May','Jun','July','Aug','Sep','Oct')
plot(x = 1:6, y = monthly_effort)
lines(x = 1:6, y = monthly_effort)

# regression modeling approach
monthly_effort2 = data.frame(boats = c(c(NA,6,7,8,5)/7,
                         c(4.75,7,7.25,4,4.75)/7.25,
                         c(5,7,9,8,8)/9,
                         c(5,8,9,8,6)/9,
                         c(6,7,9,8,5)/9,
                         c(4,7,9,6,7)/9,
                         c(5,6,6,6,6)/6 ) )
monthly_effort2$Month = rep(c(5:9), times = 7)
monthly_effort2$Year = as.factor(rep(as.character(c(2011:2017)), each = 5))

library(mgcv)
monthly_effort2_mod = gam(boats ~ s(Month, k = 4, pc = 7) - 1,
                          data = monthly_effort2,
                          offset = rep(1, times = 35))
summary(monthly_effort2_mod)
plot.gam(monthly_effort2_mod)
month_pred_df = data.frame(Month = seq(from = 5, to = 10, length.out = 100))
month_pred_df$boats = predict.gam(monthly_effort2_mod, month_pred_df, se.fit = T)$fit + 1
month_pred_df$se = predict.gam(monthly_effort2_mod, month_pred_df, se.fit = T)$se.fit
month_pred_df$LCL = month_pred_df$boats - 2*month_pred_df$se
month_pred_df$UCL = month_pred_df$boats + 2*month_pred_df$se
month_pred_plot = ggplot(month_pred_df, aes(x = Month, y = boats, ymax = UCL, ymin = LCL)) +
  geom_line() + geom_ribbon(aes(alpha = 0.4)) + 
  ggtitle('A plot showing the estimate of the relative whale watch search effort by month') +
  ylab('Relative search effort')
month_pred_plot
month_pred_df2 = data.frame(Month = 5:10)
month_pred_df2$boats = predict.gam(monthly_effort2_mod, month_pred_df2, se.fit = T)$fit + 1
month_pred_df2$se = predict.gam(monthly_effort2_mod, month_pred_df2, se.fit = T)$se.fit
month_pred_df2$LCL = month_pred_df2$boats - 2*month_pred_df2$se
month_pred_df2$UCL = month_pred_df2$boats + 2*month_pred_df2$se
month_pred_plot2 = ggplot(month_pred_df2, aes(x = Month, y = boats, ymax = UCL, ymin = LCL)) +
  geom_point() +
  geom_errorbar() + 
  ggtitle('A plot showing the estimate of the relative whale watch search effort by month') +
  ylab('Relative search effort')
month_pred_plot2
# How long are the tours?

# Victoria companies (+- 0.5 hour)
# Eagle Wing - 3.5 hours
# Prince of whales 3.5 hours
# Orca Spirit Adventures - 3 hours

# Vancouver companies - largest variability (+- 1 hour)
# Prince of Whales 4 hours (# Vancouver)
# Vancouver Whale Watch 4 hours (# Steveston)
# Wild Whales 4.5 hours (# Ganville island i.e. Vacouver)

# San Juan island companies (+- 0.5 hour)
# San Juan Safaris 3 and 3.5 hours
# Maya's legacy 2.5 - 3 hours

# Sooke companies (+- 0.5 hour)
# Sooke Whale watching - 3 hours

# Salt Spring island
# Salt spring adventures - 3.5 hours

# Port Townsend (+- 0.5 hours)
# Puget Sound Express - 4 hours and 6 hour tour

# Bellingham (+- 1 hour)
# Outer island excursions 4 hours

# Cowichan (+-0.5 hours)
# Ocean Ecoventures 3.5 hours

# Sidney 
# sidneywhalewatching 3-3.5 hours

# 17% error and 25% error respectively
port_coords = SpatialPoints(coords = rbind(
  c(-122.565408, 48.518144),#Anacortes US
  c(-122.522933, 48.748934),#Bellingham US
  c(-122.908956, 48.716899),#Brandt's landing US
  c(-123.614707, 48.746193),#Cowichan CA
  c(-123.004256, 48.611673),#Deer Harbor US
  c(-122.987170, 48.544536),#Friday Harbor US
  c(-122.943187, 48.595250),#Orcas Landing US
  c(-122.817190, 48.162493),#Port Townsend US
  c(-123.166415, 48.611249),#Roche Harbor US
  c(-123.438155, 48.831633),#Salt Spring Island CA
  c(-123.378458, 48.636657),#Sidney CA
  c(-123.174164, 48.574249),#Snug Harbor US
  c(-123.742394, 48.350035),#Sooke CA
  c(-123.200448, 49.120005),#Steveston CA
  c(-123.239652, 49.306778),#Vancouver CA
  c(-123.399219, 48.406696)#Victoria CA
),
proj4string=CRS("+init=epsg:4326"))

boundary <- list(
  as.inla.mesh.segment(COAST_simp),
  NULL)
# create a very high resolution mesh for forming the SI field
mesh_for_SI <- inla.mesh.2d(boundary=boundary,
                            max.edge=c(2000, 2000),
                            min.angle=c(21, 21),
                            max.n=c(48000, 16000), ## Safeguard against large meshes.
                            max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                            cutoff=1500, ## Filter away adjacent points.
                            offset=c(4000, 6000)) ## Offset for extra boundaries, if needed.
mesh_for_SI$crs = CRS(proj4string(COAST_simp))
mesh_for_SI = inla.spTransform(mesh_for_SI, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
polygon.triangles_SI = inla.over_sp_mesh(spTransform(COAST_simp, CRSobj = mesh_for_SI$crs ), y = mesh_for_SI, type = "centroid", ignore.CRS=T)
tl_SI = length(mesh_for_SI$graph$tv[,1])
polygon.triangles_SI = setdiff(1:tl_SI, polygon.triangles_SI)
port_coords_SI = spTransform(port_coords, CRSobj = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
posTri_SI = matrix(0, tl_SI, 2)
for (t in 1:tl_SI){
  temp = mesh_for_SI$loc[mesh_for_SI$graph$tv[t, ], ]
  posTri_SI[t,] = colMeans(temp)[c(1,2)] 
}
posTri_SI = SpatialPoints(posTri_SI, proj4string = CRS(proj4string(COAST_simp)))
plot(posTri_SI[polygon.triangles_SI,])

setwd("~/ownCloud/Whale Project/Data for Joe/Modified Shapefiles")
Sooke_Points = spTransform(readOGR('Sooke_polygon_points.shp'), CRSobj = port_coords_SI@proj4string)
Steve_Points = spTransform(readOGR('Steveston_polygon_points.shp'),CRSobj = port_coords_SI@proj4string)
Van_Points = spTransform (readOGR('Vancouver_polygon_points.shp'),CRSobj = port_coords_SI@proj4string)
Vic_Points = spTransform(readOGR('Victoria_points_after_email.shp'),CRSobj = port_coords_SI@proj4string)
port_lines = list( Points = list(Sooke_Points, Steve_Points, Van_Points, Vic_Points),
                   Range = c(8,10,10,10),
                   ind = c(13:16))
ggplot() + gg(COAST_transformed) + gg(Vic_Points)

trips_per_day = c(3,1,2,2,4,8,3,4,1,1,5,6,1,6,6,93)
trips_per_day_proportion = trips_per_day/sum(trips_per_day)
port_trip_length = c(3,4,3,3.5,3,3,3,5,3,3,3,3,3,4,4,3.5)
port_trip_length_sd = c(0.5,1,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,1,1,.5)
port_trip_range = c(66,66,60,66,60,60,60,66,60,60,66,60,66,100,100,66)

total_boat_hours = sum(trips_per_day*port_trip_length) # 506.5 hours per day 9am - 5pm
total_boat_hours_sd = sqrt(sum(port_trip_length_sd^2 * trips_per_day^2)) # 48 hours
# daily SD is 48/507 ~ 10%

port_df = SpatialPointsDataFrame(port_coords@coords,
                                 data = data.frame(trips_per_day=trips_per_day, 
                                                   trips_per_day_proportion = trips_per_day_proportion,
                                                   country = c('US','US','US','CA','US','US','US','US','US','CA','CA','US','CA','CA','CA','CA')),
                                 proj4string = CRS(proj4string(port_coords)))

port_df_plot = spTransform(port_df, CRSobj = CRS(proj4string(COAST_plotting)))

ggplot() + gg(COAST_plotting) + gg(port_df_plot[], aes(size = sqrt(trips_per_day_proportion), colour = country)) 
ggplot() + gg(mesh_barrier) + gg(port_df_plot[], aes(size = sqrt(trips_per_day_proportion), colour = country))
ggplot() + gg(port_df_plot[], aes(size = sqrt(trips_per_day_proportion), colour = country)) + gg(mesh_barrier)

sum(!gContains(COAST_plotting, port_df_plot, byid = T))
sum(!gContains(COAST_plotting, port_df_plot, byid = T))
# No ports lie on land according to the simplified computational mesh. 
# This is vital for the proper use of INLA's barrier model.

temp = by(port_df@data ,port_df@data$country, 
          FUN = function(x){x$trips_per_day/sum(x$trips_per_day)})

port_df@data$trips_per_day_proportion_by_country[port_df$country == 'US'] = temp$US
port_df@data$trips_per_day_proportion_by_country[port_df$country == 'CA'] = temp$CA 

port_df$trips_per_day[port_df$country == 'US']
port_df$trips_per_day[port_df$country == 'CA']

# 2 of the additional boats in the CA fleet are known to be Steveston. We have updated this in SI script
# need to remove the 2 boats from the growth
vessel_numbers_per_year_diff2 = vessel_numbers_per_year_diff
vessel_numbers_per_year_diff2$Can[6:8] = vessel_numbers_per_year_diff2$Can[6:8] - 2

##### IMPORTANT: Note, due to an error with xtractomatic, the above script corrupts the workspace.
# See the standalone R script wind_data_script.R to generate the data, or read in the final
# workspace below
setwd("~/ownCloud/Whale Project/Final Project Scripts")
wind = readRDS('wind.rds')

sum(duplicated(WW@data[,c('Year','Month','Day')], fromLast = FALSE)) # take the first such sighting
# 4299 repeated sightings in a day
ind_dup = duplicated(WW@data[,c('Year','Month','Day','Pod')], fromLast = FALSE) # take the first such sighting

# BG_subset = BG_subset[,c('YEAR','MONTH','PODS')]
WW_sightings2 = WW[!ind_dup,c('YEAR','MONTH','Day','Pod','TIME2')] # first sighting

# load the dates that have a never-again sighted unidentified SRKW
# See Unidentified SRKW exploratory analysis.R script for details
never_identified_dates = readRDS('never_identified_dates.rds')

# remove the never-identified sightings
WW_identified = WW[which(!(WW@data$Date %in% never_identified_dates)),]
SRKW_days = SRKW_days[!(SRKW_days %in% never_identified_dates)] 

# Assign the remaining 'SRKW' sightings to their nearest identified + verified sighting
for(i in SRKW_days)
{
  # subset sightings by date of unidentified SRKW sighting
  temp = WW_identified@data[WW_identified@data$Date == i,]
  
  # which indices are unidentified
  SRKW_ind = which(temp$Pod == 'SRKW')
  
  # which indices are identified
  non_SRKW_ind = which(temp$Pod != 'SRKW')
  
  for(j in SRKW_ind)
  {
    # find the index of the identified sighting with the closest time to the unidentified sighting
    closest_pod_sighting = which.min( abs( difftime(temp[j,'DateTimes'], temp[non_SRKW_ind,'DateTimes'])) )
    closest_pod_sighting_ind = non_SRKW_ind[closest_pod_sighting]
    
    # assign the 'Pod' status to the closest identified pod
    temp[j,'Pod'] = temp[closest_pod_sighting_ind,'Pod']
  }
  
  # assign the 'Pod' statuses to the original dataset
  WW_identified@data[WW_identified@data$Date == i,] = temp
}

# check for success
sum(WW_identified$Pod == 'SRKW') #success!

# find and remove duplicates
ind_dup_identified = duplicated(WW_identified@data[,c('Year','Month','Day','Pod')], fromLast = FALSE) # take the first such sighting
sum(ind_dup_identified) # 3437 duplicates - take the first
WW_sightings_identified2 = WW_identified[which(!ind_dup_identified),]

# Loop through and find the initial sightings of each pod and the time
WW_sightings_identified2@data$J = ifelse(
  substr(WW_sightings_identified2@data$Pod,1,1) == 'J' |
  substr(WW_sightings_identified2@data$Pod,2,2) == 'J' |
  substr(WW_sightings_identified2@data$Pod,3,3) == 'J',
  1,0
)
# check
head(WW_sightings_identified2@data)
head(WW_sightings_identified2@data[WW_sightings_identified2@data$Pod == 'K',]) # success!

WW_sightings_identified2@data$K = ifelse(
  substr(WW_sightings_identified2@data$Pod,1,1) == 'K' |
    substr(WW_sightings_identified2@data$Pod,2,2) == 'K' |
    substr(WW_sightings_identified2@data$Pod,3,3) == 'K',
  1,0
)
WW_sightings_identified2@data$L = ifelse(
  substr(WW_sightings_identified2@data$Pod,1,1) == 'L' |
    substr(WW_sightings_identified2@data$Pod,2,2) == 'L' |
    substr(WW_sightings_identified2@data$Pod,3,3) == 'L',
  1,0
)

WW_sightings_identified2@data$ALLpods = 0

# loop through the unique dates and find the first sighting time of each pod
for(i in unique(WW_sightings_identified2@data$Date))
{
  temp = WW_sightings_identified2@data[WW_sightings_identified2@data$Date == i,]
  temp$J = cumsum(temp$J)
  temp$K = cumsum(temp$K)
  temp$L = cumsum(temp$L)
  
  # create a variable counting the number of pods discovered
  temp$ALLpods = 
    min(temp$J,1) +
    min(temp$K,1) +
    min(temp$L,1)
  
  WW_sightings_identified2@data[WW_sightings_identified2@data$Date == i,] = temp
}
head(WW_sightings_identified2@data[WW_sightings_identified2@data$Pod == 'K',]) # success!
head(WW_sightings_identified2@data[WW_sightings_identified2@data$Pod == 'JK',]) # success!
head(WW_sightings_identified2@data) #success

# ALLpods gives the daily count of the number of pods discovered, assuming monotonic increasing sightings.
# J == 1 now implies the first sighting!
hist(hour(WW_sightings_identified2@data$DateTimes))

# Now we can use the J, K and L first sighting times to reduce the 
# pod-specific search effort approporiately - since we cannot find
# them again after this time, under the assumption of constant follow.

# Soundwatch 2017 reports that the bulk of whale watch boats are between 9am - 6pm 
# with most between 11am - 4pm and the peak of 1pm - 3pm.
# Let's look at histogram of FIRST sightings
hist(hour(WW[which(!(WW@data$Date %in% never_identified_dates)),]@data$DateTimes),
     main = 'A histogram of times of all and initial WW sightings', xlab = 'Hour of day')
hist(hour(WW_sightings_identified2@data$DateTimes[WW_sightings_identified2@data$J ==1 |
                                                  WW_sightings_identified2@data$K ==1 |
                                                  WW_sightings_identified2@data$L ==1]),
      add=T, col = 'red')
abline(v = median(hour(WW[which(!(WW@data$Date %in% never_identified_dates)),]@data$DateTimes)))
abline(v = median(hour(WW_sightings_identified2@data$DateTimes[WW_sightings_identified2@data$J ==1 |
                                                                 WW_sightings_identified2@data$K ==1 |
                                                                 WW_sightings_identified2@data$L ==1])),
       col = 'red')
summary(hour(WW_sightings_identified2@data$DateTimes[WW_sightings_identified2@data$J ==1 |
                                                     WW_sightings_identified2@data$K ==1 |
                                                     WW_sightings_identified2@data$L ==1]))

summary(hour(WW[which(!(WW@data$Date %in% never_identified_dates)),]@data$DateTimes))
# The Soundwatch report that almost zero sightings take place before 9am is being challenged here
# They report that sunset tours are becoming more popular. POW, 5star whales and eagle wing for example runs them.
# It is largely a new phenonmena.
# Price of whales reports 8:30 am start, Orca Spirit starts 9am, Eagle wing 9am start. 5star whales 10am start.

# Remember, we have data on the alerts. These are often made by land-based and non-WW marine observers.
# Reasonable to consider '9 - 6' as our 'working day' for WW companies?
# Could we use the truncated histogram from 9am as our density to weight the 'lost hours' by?
hist(hour(WW[which(!(WW@data$Date %in% never_identified_dates)),]@data$DateTimes)[hour(WW[which(!(WW@data$Date %in% never_identified_dates)),]@data$DateTimes) >= 9],
     main = 'A histogram of times of all verified WW sightings after 9am', xlab = 'Hour of day',
     xlim = c(8,22))
# This does not reflect the reports of the number of vessels with whale by time of day reports by Soundwatch
# Given that they report a 20 year average, we can be confident in their reports.
# Perhaps companies do not report 'follow-up' sightings with same rate as 'initial sightings'.
# This gives yet more justification for looking only at initital sightings.

# What percentage of initial sightings are made > 6pm
(sum(hour(WW_sightings_identified2@data$DateTimes[WW_sightings_identified2@data$J ==1 |
                                               WW_sightings_identified2@data$K ==1 |
                                               WW_sightings_identified2@data$L ==1])>18) /
  length(hour(WW_sightings_identified2@data$DateTimes[WW_sightings_identified2@data$J ==1 |
                                                     WW_sightings_identified2@data$K ==1 |
                                                     WW_sightings_identified2@data$L ==1])))*100   
# 0.32% of initial sightings are made after 6pm. 
# Assuming the discovery of these pods led to zero loss of search effort that day will have
# a negligible effect on the analysis. Phew...

# What percentage of initial sightings are made < 9am by pod
# IMPORTANT return to this after J_verified has been defined

#table(ifelse(hour(J_verified$DateTimes)>=9,1,0))/length(J_verified$DateTimes)
#table(ifelse(hour(K_verified$DateTimes)>=9,1,0))/length(K_verified$DateTimes)
#table(ifelse(hour(L_verified$DateTimes)>=9,1,0))/length(L_verified$DateTimes)

# between 17% and 22% of initial sightings are made before WW boats are out
# These sightings likely come from land-based observers who are difficult to estimate effort for.
# Ignoring this will underestimate search effort on the right hand side 
# The estimated intensity on the left may be biased upwards. Must report this.

# Use their 20 Year mean
SW_avg_number_boats_by_time = data.frame(avg_boats = c(0,7,12.6,14.5,13,12,16,16,12,6,0),
                                         hour = seq(from = 8.5, to = 18.5, by = 1))
SW_avg_number_boats_by_time$F = cumsum(SW_avg_number_boats_by_time$avg_boats)
SW_avg_number_boats_by_time$F = SW_avg_number_boats_by_time$F/SW_avg_number_boats_by_time$F[11]
SW_avg_number_boats_by_time$F

# create an interpolator function to estimate fraction of lost days
pred_lost_time = approxfun(x = SW_avg_number_boats_by_time$hour+0.5,
                           y = SW_avg_number_boats_by_time$F,
                           yleft = 0, yright = 1)

# note that hour is starting hour

# Loop through the sightings and remove the pod-specific search effort (as a fraction of days) 
# according to CDF.

# check the points lie in the water - move if not.
#######################################

# Check all data lies within COAST_simp
check = gContains(COAST_simp,
                  SpatialPoints(coords = cbind(WW_sightings_identified2@coords[,1], WW_sightings_identified2@coords[,2]),
                                proj4string = CRS(proj4string(COAST))),
                  byid = T)
sum(check)
sum(!check)

# 16 lie on land in WW_sightings_identified2 data
plot(COAST_simp)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col=rgb(.8,0,0,.01),cex=1.2,pch=16)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col='orangered',cex=.2,pch=16) #lies on land - move south
points(coordinates(WW_sightings_identified2)[!check,],  
       pch=16, cex=1.2, col='red')
# lies on land - move west

# Map points to nearest mesh node in hi-res mesh that lies in water
SI_mesh_sp = spTransform(SpatialPoints(coords = mesh_for_SI$loc,
                           proj4string = COAST_transformed@proj4string),
                         CRSobj = COAST_simp@proj4string)
SI_mesh_sp = SI_mesh_sp[gContains(COAST_simp, SI_mesh_sp, byid = T),]
new_coords = gDistance(WW_sightings_identified2[which(!check),],SI_mesh_sp, byid=T)
new_coords = apply(new_coords, 2 , which.min)
new_coords = SI_mesh_sp@coords[new_coords,c(1:2)]
new_coords

WW_sightings_identified2@coords[!check,] = new_coords

plot(COAST_simp)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col=rgb(.8,0,0,.01),cex=1.2,pch=16)
#points(	BG[!check2,c("Alb_East","Alb_North")]-c(0,100),col='orangered',cex=.2,pch=16) #lies on land - move south
points(coordinates(WW_sightings_identified2)[!check,],  
       pch=16, cex=1.2, col='red')

check = gContains(COAST_simp,
                  SpatialPoints(coords = cbind(WW_sightings_identified2@coords[,1], WW_sightings_identified2@coords[,2]),
                                proj4string = CRS(proj4string(COAST))),
                  byid = T)
sum(check)
sum(!check) # done

##############

# First create seperate dataframes for each Pod.
J_verified = WW_sightings_identified2[which(WW_sightings_identified2@data$J == 1),]
ind_dup_J = which(!duplicated(J_verified@data[,c('Date','J')], fromLast = FALSE)) # get the first sighting
J_verified = J_verified[ind_dup_J,]

K_verified = WW_sightings_identified2[which(WW_sightings_identified2@data$K == 1),]
ind_dup_K = which(!duplicated(K_verified@data[,c('Date','K')], fromLast = FALSE))
K_verified = K_verified[ind_dup_K,]

L_verified = WW_sightings_identified2[which(WW_sightings_identified2@data$L == 1),]
ind_dup_L = which(!duplicated(L_verified@data[,c('Date','L')], fromLast = FALSE))
L_verified = L_verified[ind_dup_L,]

# 753, 513 and 529 sightings respectively.

# What percentage of initial sightings are made < 9am by pod
table(ifelse(hour(J_verified$DateTimes)>=9,1,0))/length(J_verified$DateTimes)
table(ifelse(hour(K_verified$DateTimes)>=9,1,0))/length(K_verified$DateTimes)
table(ifelse(hour(L_verified$DateTimes)>=9,1,0))/length(L_verified$DateTimes)
# between 17% and 22% of initial sightings are made before WW boats are out
# These sightings likely come from land-based observers who are difficult to estimate effort for.
# Ignoring this will underestimate search effort on the right hand side 
# The estimated intensity on the left may be biased upwards. Must report this.

# repeat by month due to different sunrises
# compute the proportions of initial sightings made before 9am by month.
prop.table(xtabs(~ ifelse(hour(J_verified$DateTimes)>=9,1,0) + J_verified$MONTH),margin = 2)
prop.table(xtabs(~ ifelse(hour(K_verified$DateTimes)>=9,1,0) + K_verified$MONTH),margin = 2)
prop.table(xtabs(~ ifelse(hour(L_verified$DateTimes)>=9,1,0) + L_verified$MONTH),margin = 2)
# We do indeed see the proportion of initial sightings <9am increase in summer months 
# Probably due to earlier sunrise leading to an increased effectiveness of land-based spotters

sunrise_vec = c(5.5, 5.1, 5.4, 6.0, 6.78, 7.55) # hour of sunrise
# How do the initial sightings come in 
ggplot(data = data.frame(x = c(hour(J_verified$DateTimes), hour(K_verified$DateTimes),hour(L_verified$DateTimes)),
                         month = factor(c(J_verified$MONTH, K_verified$MONTH, L_verified$MONTH)),
                         pod = factor(rep(c('J','K','L'),times = c(dim(J_verified)[1], dim(K_verified)[1], dim(L_verified)[1]))),
                         sunrise = sunrise_vec[c(J_verified$MONTH-4, K_verified$MONTH-4, L_verified$MONTH-4)]),
       aes(x=x,colour = pod, group = pod, fill = pod)) + geom_histogram(breaks=seq(0,23,by=1),alpha = 0.2)+ 
  geom_vline(aes(xintercept = sunrise),data = data.frame(sunrise = sunrise_vec[c(J_verified$MONTH-4, K_verified$MONTH-4, L_verified$MONTH-4)],
                 month = factor(c(J_verified$MONTH, K_verified$MONTH, L_verified$MONTH)))) + facet_wrap(~ month) +
  ggtitle('Number of initial sightings by hour of day and by month',subtitle ='The lines show the average monthly sunrise and 9am.') + xlab('Hour of the day') +
  geom_vline(xintercept = 9) 
# The large number of sightings made in the summer months (especially for pod J), combined with their
# earlier sighting times (they are found earlier in the day in the summer) lead to the search effort being lower.

# The arrivals before 9am are clearly not exponential with constant rate
# Perhaps they could be modeled as coming from an exponential with linearly increasing rate
# Piecewise linear rate function, with maxima at 9am (i.e. constant after 9am)
plot(y = dexp(seq( from = 0, to = 3, length.out = 100), rate = 0.5*seq( from = 0, to = 3, length.out = 100)),x = seq( from = 0, to = 3, length.out = 100))
# We get interesting shape

# Assuming poisson process arrivals with constant rate over the domain. 
# This results in an extra amount of unmeasured search effort that can be estimated.
# We can form this into a sensitivity analysis?

# Based on method of moments derivations, estimate the fraction of total effort attributable to land-observers.
# Ignore truncation for the moment. In reality, we truncate at sunset.
p_1 = prop.table(xtabs(~ ifelse(hour(c(J_verified$DateTimes,K_verified$DateTimes,L_verified$DateTimes))>=9,1,0) + c(J_verified$MONTH,K_verified$MONTH,L_verified$MONTH)),margin = 2)[2,]
p_2 = prop.table(xtabs(~ ifelse(hour(c(J_verified$DateTimes,K_verified$DateTimes,L_verified$DateTimes))>=9 &
                                hour(c(J_verified$DateTimes,K_verified$DateTimes,L_verified$DateTimes))<15,1,0) + c(J_verified$MONTH,K_verified$MONTH,L_verified$MONTH)),margin = 2)[2,]
p_1
p_2
# estimate p_hat using method of moments
sunrise_vec = c(5.5, 5.1, 5.4, 6.0, 6.78, 7.55) # hour of sunrise
daylight_time = 9-sunrise_vec # number of hours of daylight for land-based observers
plot(y=p_1,x=daylight_time) # No obvious trend with number of hours of daylight and observer time.
# Perhaps people wake up at the same time month by month. Try 7 am wake up time
awake_time = 9 - pmax(rep(7,6),sunrise_vec) # looks at the number of hours of sunlight between avg wake up and 9am
plot(y=p_1,x=awake_time)
alpha_hat = -(1/6)*log(1 - p_2/p_1)
p_hat = -log(p_1)/(alpha_hat*(awake_time))
p_hat

# compute the total number of hours after waking up at 7am that land-observers are active

# What percentage of initial sightings are made > 6pm by pod
table(ifelse(hour(J_verified$DateTimes)<18,1,0))/length(J_verified$DateTimes)
table(ifelse(hour(K_verified$DateTimes)<18,1,0))/length(K_verified$DateTimes)
table(ifelse(hour(L_verified$DateTimes)<18,1,0))/length(L_verified$DateTimes)
# between 0.2% and 1.4% of sightings. Small effect

# Are the early sightings (<9am) being made near population centers?
ggplot() + gg(COAST_plotting) +
  gg(J_verified[hour(J_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(K_verified[hour(K_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(L_verified[hour(L_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(J_verified[hour(J_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  gg(K_verified[hour(K_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  gg(L_verified[hour(L_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  ggtitle('A plot of first sightings made before 9am') + xlab('') + ylab('')

# The spatial coverage of early sightings is similar to the other sightings, 
# however they are made closer to land/shore.

# A crude correction for this could be to estimate the land-based search effort 
# as a constant poisson process with identical spatial coverage/intensity
# Then add this to the total WW SI to see how results change.

# Are the early sightings (<9am) being made near population centers?
ggplot() + gg(COAST_plotting) +
  gg(J_verified[hour(J_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(K_verified[hour(K_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(L_verified[hour(L_verified$DateTimes)>=9,],aes(colour = '> 9am')) +
  gg(J_verified[hour(J_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  gg(K_verified[hour(K_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  gg(L_verified[hour(L_verified$DateTimes)<9,],aes(colour = '< 9am')) +
  ggtitle('A plot of first sightings made before 9am') + xlab('') + ylab('')

# Histogram of distance-to-shore
temp.triangles = inla.over_sp_mesh(COAST_simp, y = mesh_barrier, type = "centroid", ignore.CRS=T)
temp.triangles = setdiff(1:tl, polygon.triangles)
# - checking which mesh triangles are inside the normal area
plot(posTri[temp.triangles,])
temp_dist_J_pre9 = apply( gDistance(J_verified[hour(J_verified$DateTimes)<9,], posTri[temp.triangles,], byid = T), 2, min )
temp_dist_J_post9 = apply( gDistance(J_verified[hour(J_verified$DateTimes)>=9,], posTri[temp.triangles,], byid = T), 2, min )
temp_dist_K_pre9 = apply( gDistance(K_verified[hour(K_verified$DateTimes)<9,], posTri[temp.triangles,], byid = T), 2, min )
temp_dist_K_post9 = apply( gDistance(K_verified[hour(K_verified$DateTimes)>=9,], posTri[temp.triangles,], byid = T), 2, min )
temp_dist_L_pre9 = apply( gDistance(L_verified[hour(L_verified$DateTimes)<9,], posTri[temp.triangles,], byid = T), 2, min )
temp_dist_L_post9 = apply( gDistance(L_verified[hour(L_verified$DateTimes)>=9,], posTri[temp.triangles,], byid = T), 2, min )

plot(density(temp_dist_J_pre9), main =''); lines(density(temp_dist_J_post9), add=T,col = 'red'); title('Distance to shore of sightings of J before and after 9am \n Red is after 9am')
plot(density(temp_dist_K_pre9), main =''); lines(density(temp_dist_K_post9), add=T,col = 'red'); title('Distance to shore of sightings of J before and after 9am \n Red is after 9am')
plot(density(temp_dist_L_pre9), main =''); lines(density(temp_dist_L_post9), add=T,col = 'red'); title('Distance to shore of sightings of J before and after 9am \n Red is after 9am')

# The spatial coverage of early sightings is similar to the other sightings, 
# They are made no closer to land/shore.

# A crude correction for this could be to estimate the land-based search effort 
# as a constant poisson process with identical spatial coverage/intensity
# Then add this to the total WW SI to see how results change.

##### Remove duplicated sightings from BG_subset

#save.image('Temp_WS_dup_removed.RData')

################# SAME TREATMENT FOR BG DATA
# Loop through and find the initial sightings of each pod and the time

# remove NA pods first
BG_follows = BG_DF2
BG_follows = BG_follows[which(!is.na(BG_follows@data$PODS)),]
BG_subset_identified = BG_subset[which(!is.na(BG_subset@data$Pod)),]

BG_subset_identified@data$J = ifelse(
  substr(BG_subset_identified@data$Pod,1,1) == 'J' |
    substr(BG_subset_identified@data$Pod,2,2) == 'J' |
    substr(BG_subset_identified@data$Pod,3,3) == 'J',
  1,0
)
# check
head(BG_subset_identified@data)
head(BG_subset_identified@data[BG_subset_identified@data$Pod == 'K',]) # success!

BG_subset_identified@data$K = ifelse(
  substr(BG_subset_identified@data$Pod,1,1) == 'K' |
    substr(BG_subset_identified@data$Pod,2,2) == 'K' |
    substr(BG_subset_identified@data$Pod,3,3) == 'K',
  1,0
)
BG_subset_identified@data$L = ifelse(
  substr(BG_subset_identified@data$Pod,1,1) == 'L' |
    substr(BG_subset_identified@data$Pod,2,2) == 'L' |
    substr(BG_subset_identified@data$Pod,3,3) == 'L',
  1,0
)

BG_subset_identified@data$ALLpods = 0

BG_subset_identified$Date = as.Date(BG_subset_identified$TIME2)
# loop through the unique dates and find the first sighting time of each pod
for(i in unique(BG_subset_identified@data$Date))
{
  temp = BG_subset_identified@data[BG_subset_identified@data$Date == i,]
  temp$J = cumsum(temp$J)
  temp$K = cumsum(temp$K)
  temp$L = cumsum(temp$L)
  
  # create a variable counting the number of pods discovered
  temp$ALLpods = 
    min(temp$J,1) +
    min(temp$K,1) +
    min(temp$L,1)
  
  BG_subset_identified@data[BG_subset_identified@data$Date == i,] = temp
}
head(BG_subset_identified@data[BG_subset_identified@data$Pod == 'K',]) # success!
head(BG_subset_identified@data[BG_subset_identified@data$Pod == 'JK',]) # success!
head(BG_subset_identified@data) #success

# First create seperate dataframes for each Pod. Take the FIRST unique POD sighting
J_verified_Brian = BG_subset_identified[which(BG_subset_identified@data$J == 1),]
ind_dup_J_Brian = which(!duplicated(J_verified_Brian@data[,c('Date','J')], fromLast = FALSE)) # get the first sighting
J_verified_Brian = J_verified_Brian[ind_dup_J_Brian,]

K_verified_Brian = BG_subset_identified[which(BG_subset_identified@data$K == 1),]
ind_dup_K_Brian = which(!duplicated(K_verified_Brian@data[,c('Date','K')], fromLast = FALSE))
K_verified_Brian = K_verified_Brian[ind_dup_K_Brian,]

L_verified_Brian = BG_subset_identified[which(BG_subset_identified@data$L == 1),]
ind_dup_L_Brian = which(!duplicated(L_verified_Brian@data[,c('Date','L')], fromLast = FALSE))
L_verified_Brian = L_verified_Brian[ind_dup_L_Brian,]

##################

# merge the verified sightings from the two sources for each pod to remove WW search effort
temp_J_verified_combined = rbind(J_verified[,c('Date','YEAR','MONTH','DAY','TIME2')], J_verified_Brian[,c('Date','YEAR','MONTH','DAY','TIME2')])
temp_K_verified_combined = rbind(K_verified[,c('Date','YEAR','MONTH','DAY','TIME2')], K_verified_Brian[,c('Date','YEAR','MONTH','DAY','TIME2')])
temp_L_verified_combined = rbind(L_verified[,c('Date','YEAR','MONTH','DAY','TIME2')], L_verified_Brian[,c('Date','YEAR','MONTH','DAY','TIME2')])

# Keep only the first sighting
temp_J_verified_combined_dates = unique(temp_J_verified_combined$Date)
temp_K_verified_combined_dates = unique(temp_K_verified_combined$Date)
temp_L_verified_combined_dates = unique(temp_L_verified_combined$Date)

count_J = 1
for(i in temp_J_verified_combined_dates)
{
  ind_dates_J = which(temp_J_verified_combined$Date == i)
  ind_dates_J_first = which.min( temp_J_verified_combined$TIME2[ind_dates_J] )
  if(count_J == 1)
  {
    ind_dates_J_total = ind_dates_J[ind_dates_J_first]
  }
  if(count_J != 1)
  {
    ind_dates_J_total = c(ind_dates_J_total, ind_dates_J[ind_dates_J_first])  
  }
  count_J = count_J+1
}
dim(temp_J_verified_combined)[1] - length(ind_dates_J_total)
# removed 1 duplicate - correct based on exploratory analysis
temp_J_verified_combined = temp_J_verified_combined[ind_dates_J_total,]

count_K = 1
for(i in temp_K_verified_combined_dates)
{
  ind_dates_K = which(temp_K_verified_combined$Date == i)
  ind_dates_K_first = which.min( temp_K_verified_combined$TIME2[ind_dates_K] )
  if(count_K == 1)
  {
    ind_dates_K_total = ind_dates_K[ind_dates_K_first]
  }
  if(count_K != 1)
  {
    ind_dates_K_total = c(ind_dates_K_total, ind_dates_K[ind_dates_K_first])  
  }
  count_K = count_K + 1
}
dim(temp_K_verified_combined)[1] - length(ind_dates_K_total)
# removed 3 correct
temp_K_verified_combined = temp_K_verified_combined[ind_dates_K_total,]

count_L = 1
for(i in temp_L_verified_combined_dates)
{
  ind_dates_L = which(temp_L_verified_combined$Date == i)
  ind_dates_L_first = which.min( temp_L_verified_combined$TIME2[ind_dates_L] )
  if(count_L == 1)
  {
    ind_dates_L_total = ind_dates_L[ind_dates_L_first]
  }
  if(count_L != 1)
  {
    ind_dates_L_total = c(ind_dates_L_total, ind_dates_L[ind_dates_L_first])  
  }
  count_L = count_L + 1
}
dim(temp_L_verified_combined)[1] - length(ind_dates_L_total)
# removed 51 duplicates - correct due to duplicated L sightings per day from BG
temp_L_verified_combined = temp_L_verified_combined[ind_dates_L_total,]

# how many days (as a fraction) of search effort do we need to remove per month, per pod, per year
J_days_removed_month_year = matrix(0, nrow = 6, ncol = 8)
K_days_removed_month_year = matrix(0, nrow = 6, ncol = 8)
L_days_removed_month_year = matrix(0, nrow = 6, ncol = 8)
count = 1
for(i in 5:10)
{
  count2 = 1
  for(j in 2009:2016)
  {
  ind_J_month = which(temp_J_verified_combined@data$MONTH == i & temp_J_verified_combined@data$YEAR == j)
  J_month = temp_J_verified_combined[ind_J_month,]
  
  ind_K_month = which(temp_K_verified_combined@data$MONTH == i & temp_K_verified_combined@data$YEAR == j)
  K_month = temp_K_verified_combined[ind_K_month,]
  
  ind_L_month = which(temp_L_verified_combined@data$MONTH == i & temp_L_verified_combined@data$YEAR == j)
  L_month = temp_L_verified_combined[ind_L_month,]
  
  initial_times_J = strftime(J_month@data$TIME2, 
                           format = '%H:%M:%S', 
                           tz = 'UTC')
  initial_times_K = strftime(K_month@data$TIME2, 
                             format = '%H:%M:%S', 
                             tz = 'UTC')
  initial_times_L = strftime(L_month@data$TIME2, 
                             format = '%H:%M:%S', 
                             tz = 'UTC')
  
  start_time = strptime('09:00:00','%H:%M:%S',tz = 'UTC')
  
  hours_into_the_day_J = pmax( as.numeric(difftime(strptime(initial_times_J, 
                                                    format = '%H:%M:%S', 
                                                    tz = 'UTC'),start_time, 
                                           units = 'hours')),
                              0 )
  
  hours_into_the_day_K = pmax( as.numeric(difftime(strptime(initial_times_K, 
                                                           format = '%H:%M:%S', 
                                                           tz = 'UTC'),start_time, 
                                                  units = 'hours')),
                              0 )
  hours_into_the_day_L = pmax( as.numeric(difftime(strptime(initial_times_L, 
                                                           format = '%H:%M:%S', 
                                                           tz = 'UTC'),start_time, 
                                                  units = 'hours')),
                              0 )
  # Use the estimated CDF to predict the fraction of daily search effort lost
  J_days_removed_month_year[count,count2] = sum(1 - pred_lost_time(9 + hours_into_the_day_J))
  K_days_removed_month_year[count,count2] = sum(1 - pred_lost_time(9 + hours_into_the_day_K))
  L_days_removed_month_year[count,count2] = sum(1 - pred_lost_time(9 + hours_into_the_day_L))
  
  count2 = count2 + 1
  }

  count = count + 1
}

J_days_removed_month_year
K_days_removed_month_year
L_days_removed_month_year

rowSums(J_days_removed_month_year)
rowSums(K_days_removed_month_year)
rowSums(L_days_removed_month_year)
# Up to 26 days of effort removed for a single pod in a single month!
# Around an 8th of the search effort.
rowSums(J_days_removed_month_year)/8
rowSums(K_days_removed_month_year)/8
rowSums(L_days_removed_month_year)/8
# an average of up to 18 days per month lost!

# Define our uncertainties in various effort layers as probability distributions.
# 1 - Uncertainty regarding relative monthly search effort:
Uncertainty_month = call('rnorm', 6, 
                         mean = month_pred_df2$boats, 
                         sd = month_pred_df2$se)
# Note that this stores the function call as an object. To evaluate (i.e. call) 
# the function, simply run eval(Uncertainty_month)

# 2 - Uncertainty regarding cancelled days due to weather
cancellations_simulator = function(min, max){
  # input matrices of min and max for uniform distribution. 
  # rows are months and columns are years
  tmp = runif(n = prod(dim(min)),
              min = as.numeric(min),
              max = as.numeric(max))
  tmp = matrix(round(tmp), nrow = nrow(min), ncol = ncol(min))
  return(tmp)
}

Uncertainty_cancellations = call('cancellations_simulator',
                                 min = xtabs(~month + year + BSseven, data = wind)[,,2],
                                 max = xtabs(~month + year + BSsix, data = wind)[,,2])

# 3 - Uncertainty regarding trip length from each port
Uncertainty_triptimes = call('runif',16,
                             min = port_trip_length - port_trip_length_sd,
                             max = port_trip_length + port_trip_length_sd)

# 4 - Uncertainty regarding (max) number of trips per day, per port, per year 
# vessel_numbers_per_year_diff is the change in US and Canadian trips per day across the years relative to 2011 baseline
# port_df contains the number of trips per day from each port in 2011, the proportion of trips per day per port, split by country
#sum(port_df$country == 'US') = 9 and sum(port_df$country == 'CA') = 7
number_trips_simulator = function(sizes_year, prob, n){
  # input number of incre
  number_trips = matrix(0, nrow = length(prob), ncol = length(sizes_year))
  for(i in 1:length(sizes_year)){
    tmp_var = as.numeric(rmultinom(n=n, size = abs(sizes_year[i]), prob = prob))
    if(sizes_year[i]<0)
    {
      tmp_var = -1*tmp_var
    }
    number_trips[,i] = tmp_var
  }
  return(number_trips)
}

Uncertainty_number_trips_US = call('number_trips_simulator', 
                                n = 1,
                                sizes_year = vessel_numbers_per_year_diff$US ,
                                prob = port_df@data[port_df@data$country == 'US',]$trips_per_day_proportion)
#Uncertainty_number_trips_CA = call('number_trips_simulator', 
#                                   n = 1,
#                                   sizes_year = vessel_numbers_per_year_diff$Can ,
#                                   prob = port_df@data[port_df@data$country == 'CA',]$trips_per_day_proportion)

Uncertainty_number_trips_CA = call('number_trips_simulator', 
                                   n = 1,
                                   sizes_year = vessel_numbers_per_year_diff2$Can ,
                                   prob = port_df@data[port_df@data$country == 'CA',]$trips_per_day_proportion)


# 5 - Uncertainty regarding number of hours spent with SRKW + wildlife per trip
# Estimate the TOTAL number of hours lost spent tracking whale per trip, per discovery.
#time_stationary_per_trip = call('runif',
#                                n = 1, min = 0.5, max = 0.75)
# NOT NEEDED

# create sampling intensity surface for WW boats
# load neccessary functions
source('sampling_intensity_script.R')

# THIS IS SLOW - RUN AT OWN RISK - Alternatively load pre-compiled file below
SI_WW2 = sampling_intensity_script(barrier_mesh = mesh_for_SI, barrier_triangles = polygon.triangles_SI, 
                                   port_coords = port_coords_SI, port_number_trips = trips_per_day, 
                                   port_countries = port_df$country, port_trip_length = Uncertainty_triptimes, 
                                   port_trip_range = port_trip_range, barrier_range = 4,
                                   number_time_periods = 6, relative_effort_period = Uncertainty_month,
                                   annual_effort_relative = list(US = Uncertainty_number_trips_US,
                                                                 CA = Uncertainty_number_trips_CA), 
                                   number_years = 8, number_days_per_period = c(31,30,31,31,30,31),
                                   cancelled_days_per_period = Uncertainty_cancellations, 
                                   lost_days_from_follow = list(J_days_removed_month_year, 
                                                                K_days_removed_month_year,
                                                                L_days_removed_month_year),
                                   polygons = COAST_transformed, no_MC_samples = 1000, port_sim = 'yes',
                                   port_lines = port_lines)
# scale the WW SI by the number of boat hours on the water

#saveRDS(SI_WW2,'SI_WW2_dup_removed.rds')
setwd("~/ownCloud/Whale Project/Final Project Scripts")
SI_WW2 = readRDS('SI_WW2_dup_removed.rds')

SI_WW = SI_WW2
#Anacortes US #Bellingham US #Brandt's landing US
#Cowichan CA #Deer Harbor US #Friday Harbor US 
#Orcas Landing US #Port Townsend US #Roche Harbor US
#Salt Spring Island CA #Sidney CA #Snug Harbor US #Sooke CA #Steveston CA #Vancouver CA #Victoria CA
multiplot(SI_WW$plotting_fields[[1]]+ ggtitle('Anacortes US sampling field'),
          SI_WW$plotting_fields[[2]]+ ggtitle('Bellingham US sampling field'),
          SI_WW$plotting_fields[[3]]+ ggtitle('Brandts landing US sampling field'),
          SI_WW$plotting_fields[[4]]+ ggtitle('Cowichan CA sampling field'),
          SI_WW$plotting_fields[[5]]+ ggtitle('Deer Harbor US sampling field'),
          SI_WW$plotting_fields[[6]]+ ggtitle('Friday Harbor US sampling field'),
          cols = 3)
multiplot(SI_WW$plotting_fields[[7]]+ ggtitle('Orcas Landing US sampling field'),
          SI_WW$plotting_fields[[8]]+ ggtitle('Port Townsend US sampling field'),
          SI_WW$plotting_fields[[9]]+ ggtitle('Roche Harbor US sampling field'),
          SI_WW$plotting_fields[[10]]+ ggtitle('Salt Spring Island CA sampling field'),
          SI_WW$plotting_fields[[11]]+ ggtitle('Sidney CA sampling field'),
          SI_WW$plotting_fields[[12]]+ ggtitle('Snug Harbor US sampling field'),
          cols = 3)
multiplot(SI_WW$plotting_fields[[13]]+ ggtitle('Sooke CA sampling field'),
          SI_WW$plotting_fields[[14]]+ ggtitle('Steveston CA sampling field'),
          SI_WW$plotting_fields[[15]]+ ggtitle('Vancouver CA sampling field'),
          SI_WW$plotting_fields[[16]]+ ggtitle('Victoria CA sampling field'),
          cols = 2)

total_boat_hours_per_period_per_port_J = array(0, dim = c(dim(SI_WW2$total_boat_hours_per_period_per_port_J[[1]])[1],
                                                        dim(SI_WW2$total_boat_hours_per_period_per_port_J[[1]])[2],
                                                        1000))
total_boat_hours_per_period_per_port_K = array(0, dim = c(dim(SI_WW2$total_boat_hours_per_period_per_port_K[[1]])[1],
                                                          dim(SI_WW2$total_boat_hours_per_period_per_port_K[[1]])[2],
                                                          1000))
total_boat_hours_per_period_per_port_L = array(0, dim = c(dim(SI_WW2$total_boat_hours_per_period_per_port_L[[1]])[1],
                                                          dim(SI_WW2$total_boat_hours_per_period_per_port_L[[1]])[2],
                                                          1000))
for(i in 1:1000)
{
  total_boat_hours_per_period_per_port_J[,,i] = SI_WW2$total_boat_hours_per_period_per_port_J[[i]]
  total_boat_hours_per_period_per_port_K[,,i] = SI_WW2$total_boat_hours_per_period_per_port_K[[i]]
  total_boat_hours_per_period_per_port_L[,,i] = SI_WW2$total_boat_hours_per_period_per_port_L[[i]]
}
sd_total_boat_hours_per_period_per_port_J = apply(total_boat_hours_per_period_per_port_J,
                                                c(1,2), sd, na.rm=T)
sd_total_boat_hours_per_period_per_port_K = apply(total_boat_hours_per_period_per_port_K,
                                                  c(1,2), sd, na.rm=T)
sd_total_boat_hours_per_period_per_port_L = apply(total_boat_hours_per_period_per_port_L,
                                                  c(1,2), sd, na.rm=T)
sd_total_boat_hours_per_period_per_port_J
mean_total_boat_hours_per_period_per_port_J = apply(total_boat_hours_per_period_per_port_J,
                                                c(1,2), mean, na.rm=T)
mean_total_boat_hours_per_period_per_port_K = apply(total_boat_hours_per_period_per_port_K,
                                                    c(1,2), mean, na.rm=T)
mean_total_boat_hours_per_period_per_port_L = apply(total_boat_hours_per_period_per_port_L,
                                                    c(1,2), mean, na.rm=T)
mean_total_boat_hours_per_period_per_port_J

# view the CV of the ports and months
sd_total_boat_hours_per_period_per_port_J / mean_total_boat_hours_per_period_per_port_J
rowSums(sd_total_boat_hours_per_period_per_port_J / mean_total_boat_hours_per_period_per_port_J)
# cv ranges from 0.39 for Victoria to 1.05 to Roche Harbor.

# View the empirical distribution of the values for a given port and month across the MC samples
ggplot(data = data.frame(boat_hours = total_boat_hours_per_period_per_port_J[1,1,]),
       aes(x = boat_hours)) + geom_density() + ggtitle('Anacortes US, May Monte Carlo boat hours')
ggplot(data = data.frame(boat_hours = total_boat_hours_per_period_per_port_J[9,6,]),
       aes(x = boat_hours)) + geom_density() + ggtitle('Roche Harbor US, October Monte Carlo boat hours')
ggplot(data = data.frame(boat_hours = total_boat_hours_per_period_per_port_J[16,3,]),
       aes(x = boat_hours)) + geom_density() + ggtitle('Victoria CA, July Monte Carlo boat hours')

# Estimate sampling intensity surface for Brian's boat
# Idea is to fit a continuous time random walk (zero measurement error) to Brian's 
# boat tracks and use it to predict location every 30 seconds say.
# Then share the total boat hours over the points and count the number of points
# falling inside each dual mesh grid cell.

#save.image('SI_WS_dup_removed.RData')

# longform.points are Brian's movement points
library(crawl)

# reload the longform points - adding the SEG_INDEX variable
# longform.points <- NULL
# for (year in c("09","10","11","12","13","14","15","16")) {
#   # Mac edition  
#   #set.working.dir <- paste("/Users/joe/Documents/Whale Project/Data for Joe/CRP-BG EFFORT DATA/20",year,"/DAYLIGHT_trackpoints",sep="")
#   # PC verison
#   set.working.dir <- paste("~/ownCloud/Whale Project/Data for Joe/CRP-BG EFFORT DATA/20",year,"/DAYLIGHT_trackpoints",sep="")
#   setwd(set.working.dir)
#   temp <- list.files(pattern="*.dbf", full.names=TRUE)
#   filenames <- unlist(lapply(strsplit(unlist(strsplit(temp,".dbf")),"./"),"[[",2))
#   if (year=="09") ind=(1:length(filenames))[-c(36,101)] # take out the run up the east side of Vanc.Isl (36)
#   if (year!="09") ind=(1:length(filenames))
#   ind=(1:length(filenames))
#   for (i in ind) {
#     temp.files=readOGR(".",layer=filenames[i])
#     #assign(paste(substr(filenames[i],7,28),"_", i, sep=""), temp.files)
#     
#     longform.points <- rbind(longform.points , 
#                              data.frame(x= temp.files@coords[,1],y= temp.files@coords[,2], 
#                                         Date=temp.files@data$Day,
#                                         Month = as.integer(substr(as.character(temp.files@data$Day),6,7)),
#                                         Time=substr(as.character(temp.files@data$Date_PDT),12,19),
#                                         Day = substr(as.character(temp.files@data$Date_PDT),9,10),
#                                         Seg_Index = temp.files@data$Seg_Index) )
#     
#     if ((i%%5)==0)	print(paste("year",year,"trackline",i,sep=" - "))
#     print(i)
#     print('out of')
#     print(length(ind))
#   }
# }


longform.points.sp = SpatialPointsDataFrame(coords = cbind(longform.points$x,
                                                           longform.points$y),
                                            data = data.frame(Time = hms(longform.points$Time),
                                                              Date = ymd(longform.points$Date),
                                                              Datetime = as.POSIXlt(paste(longform.points$Date, longform.points$Time)),
                                                              Seg_Index = longform.points$Seg_Index),
                                            proj4string = COAST_simp@proj4string)
               
start_times = rep(0, length(unique(longform.points.sp$Date)))
end_times = rep(0, length(unique(longform.points.sp$Date)))                                               

count = 1
for(i in unique(longform.points.sp$Date))
{
  # find the indices for that date
  ind_BG = which(longform.points.sp@data$Date == i)
  dat_daily = longform.points.sp[ind_BG,]
  # loop through the GPS segments
  for(j in unique(dat_daily@data$Seg_Index))
  {
    # find the indices for that date and GPS segment
    ind_BG2 = which(longform.points.sp@data$Date == i & longform.points.sp@data$Seg_Index == j)
    if(length(ind_BG2) > 2) # remove singletons
    {
      dat_daily2 = longform.points.sp[ind_BG2,]
      
      start_times[count] = as.character(min(dat_daily2$Datetime))
      end_times[count] = as.character(max(dat_daily2$Datetime))
      count = count+1
    }
  }
}

# How long was his journey
sum(difftime(end_times,start_times,units = 'hours')) 
# 6411.575 hours including all segments (i.e. with or without singletons)
# 6400.626 hours if we discard all segments of length < 3.
# lose around 0.1% of total search effort. Can be sure of quality though.
hist(as.numeric(difftime(end_times,start_times,units = 'hours')))
min(difftime(end_times,start_times,units = 'secs'))

# All segments have more than 30 seconds effort

View(cbind(start_times, end_times))


# loop over the 816 dates, and loop over the GPS segment indices.
# fit the CTRW model, predict on 30 second intervals and create SpatialPointsDataframe
# Repeat for the 'follows' and subtract the appropriate percentage of points from the dual mesh cells.
# WARNING SLOW - RECOMMENDED TO LOAD PRE-COMPILED OBJECT

# create objects first to speed up storage

# how many 30 second intervals do we need at most (i.e. rows of a matrix?)

BG_total_boathours = 0
BG_thirtysecond_locations = vector('list',length(unique(longform.points$Seg_Index)))
BG_predicted_paths = vector('list',length(unique(longform.points$Seg_Index)))
count = 1
count2 = 1
for(i in unique(longform.points.sp$Date))
{
  # find the indices for that date
  ind_BG = which(longform.points.sp@data$Date == i)
  dat_daily = longform.points.sp[ind_BG,]
  
  # loop through the GPS segments
  for(j in unique(dat_daily@data$Seg_Index))
  {
    # find the indices for that date and GPS segment
    ind_BG2 = which(longform.points.sp@data$Date == i & longform.points.sp@data$Seg_Index == j)
    if(length(ind_BG2) > 2) # remove singletons and paths with 2 points
    {
      dat_daily2 = longform.points.sp[ind_BG2,]
  
      # how long was his recorded trip?
      diff_time = dseconds(as.interval(min(dat_daily2$Datetime),dat_daily2$Datetime))
      trip_time = max(as.numeric(diff_time))
      BG_total_boathours = BG_total_boathours + trip_time
  
      start_time = min(dat_daily2$Datetime)
      end_time = max(dat_daily2$Datetime)
      
      # Has the boat moved?
      if(sum(abs(diff(dat_daily2@coords)))> 0)
      {
        
        # the below code can fail - catch in loop and retry if so
        function_may_fail = function()
        {
          # fit the model
          BG_mod = crwMLE(data = dat_daily2,
                          Time.name = 'Datetime',
                          time.scale = 'seconds',
                          initialSANN = list(maxit = 200 + (attempt*30)),
                          attempts = 3,
                          #theta = c(log(100),log(0.8)),
                          tryBrownian=T,
                          need.hess = FALSE)
          
          # predict location on a 30 second time scale
          BG_pred = crwPredict(BG_mod, seq.POSIXt(from = start_time, to = end_time, by = '30 sec'), 
                               return.type = 'flat', '30 seconds')
          
          # plot path
          plot_path = crwPredictPlot(BG_pred, 'map')
          # remove observed locations
          BG_pred = BG_pred[BG_pred$locType == 'p',]
          
          return(list(BG_pred, plot_path))
        }
        
        temp = NULL
        attempt = 0
        while(is.null(temp) && attempt <= 30)
        {
          attempt = attempt + 1
          if(attempt > 1)
          {
            print(paste('attempt number',attempt))
          }
          try(
            temp <- function_may_fail()
          )
        }
        
        BG_pred = temp[[1]]
        BG_predicted_paths[[count]] = temp[[2]]

      }
      if(sum(abs(diff(dat_daily2@coords))) == 0)
      {
        # boat is stationary - simply repeat the location appropriately
        no_intervals = round(trip_time / 30) # how many repetitions
        DateTimes_pred = seq.POSIXt(from = start_time, by = '30 sec', length.out = no_intervals)
        
        BG_pred = matrix(0, nrow = no_intervals, ncol = 17)
        colnames(BG_pred) = c("TimeNum","locType","Time","Date","Datetime","Seg_Index","coords.x1","coords.x2","mu.x",     
                              "nu.x","mu.y","nu.y","se.mu.x","se.nu.x","se.mu.y","se.nu.y","speed")
        BG_pred = as.data.frame(BG_pred)
        
        BG_pred$Datetime = DateTimes_pred 
        BG_pred$mu.x = dat_daily2@coords[1,1]
        BG_pred$mu.y = dat_daily2@coords[1,2]
        BG_pred$speed = 0
        BG_pred$Seg_Index = unique(dat_daily2@data$Seg_Index)
        BG_pred$locType = 'p'
        
        # plot path/point
        BG_predicted_paths[[count]] = plot(BG_pred$mu.x, BG_pred$mu.y, main = paste('stationary for',trip_time,'seconds'))
      }
      
      BG_thirtysecond_locations[[count]] = BG_pred
  
      # if(count == 1)
      # {
      #   BG_thirtysecond_locations = BG_pred
      # }
      # if(count!=1)
      # {
      #   BG_thirtysecond_locations = merge(BG_thirtysecond_locations, BG_pred, all.x=T, all.y = T) 
      # }
      
    count = count+1
    }
  }
  print(paste('Iteration',count2,'out of',length(unique(longform.points.sp$Date))))
  count2 = count2 + 1
}

# choose which variables to keep
BG_thirtysecond_locations2 = BG_thirtysecond_locations
for(i in 1:length(BG_thirtysecond_locations2))
{
  BG_thirtysecond_locations2[[i]] = BG_thirtysecond_locations2[[i]][,c('Datetime', 'mu.x','mu.y')]
}

BG_thirtysecond_locations2 = do.call(rbind, BG_thirtysecond_locations2)
# saveRDS(list(BG_total_boathours = BG_total_boathours,
#              BG_thirtysecond_locations = BG_thirtysecond_locations2,
#              BG_predicted_paths = BG_predicted_paths),
#         'BG_paths.rds')

# read in the precompiled stuff to save time
BG_paths = readRDS('BG_paths.rds')
BG_total_boathours = BG_paths$BG_total_boathours
BG_thirtysecond_locations2 = BG_paths$BG_thirtysecond_locations
BG_predicted_paths = BG_paths$BG_predicted_paths
rm(BG_paths)

BG_total_boathours = BG_total_boathours/(60*60) # convert from seconds to hours
# check the math
(dim(BG_thirtysecond_locations2)[1] / 2)/60
BG_total_boathours
#  slight discrepancy (0.2%) due to rounding to nearest 30 second for prediction.

# Convert to a SpatialPoints object and keep ONLY the points lying in AOI and in water.
BG_thirtysecond_locations_sp = SpatialPointsDataFrame(coords = cbind(BG_thirtysecond_locations2$mu.x,
                                                               BG_thirtysecond_locations2$mu.y),
                                                      data = data.frame(Datetime = BG_thirtysecond_locations2$Datetime,
                                                                        month = month(BG_thirtysecond_locations2$Datetime),
                                                                        Date = as.Date(BG_thirtysecond_locations2$Datetime) ),
                                             proj4string = CRS(proj4string(COAST)))

######## Jitter all points lying on land


#######################  START
# Simply count the number of points lying in each dmesh grid.
# First remove all points not lying in AOI and/or on land
BG_remove = gContains(AOI, BG_thirtysecond_locations_sp, byid = T)
sum(!BG_remove) / length(BG_remove) # 15.6% lie outside AOI
BG_thirtysecond_locations_sp = BG_thirtysecond_locations_sp[which(BG_remove),]
ggplot() + gg(COAST_simp) + gg(BG_thirtysecond_locations_sp)

BG_remove = gContains(COAST_simp, BG_thirtysecond_locations_sp, byid = T)
sum(!BG_remove) / length(BG_remove) # 13.7% lie on land

# Check it isn't a resolution issue
BG_remove2 = gContains(COAST_plotting, BG_thirtysecond_locations_sp, byid = T)
sum(!BG_remove2) / length(BG_remove2) # 13.6% lie on land - better

BG_thirtysecond_locations_sp2 = BG_thirtysecond_locations_sp

BG_thirtysecond_locations_sp2$inwater = as.numeric(BG_remove)
BG_thirtysecond_locations_sp2$inwater2 = as.numeric(BG_remove2)

ggplot() + gg(COAST_simp) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater))
ggplot() + gg(COAST_plotting) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater2))

# It appears some of the predicted paths lie slightly too far north east. Jitter
BG_remove3 = BG_thirtysecond_locations_sp2@coords[,1] > 1050000 & BG_thirtysecond_locations_sp2@coords[,1] < 1150000 & !BG_remove2
BG_remove3 = BG_remove3 & BG_thirtysecond_locations_sp2@coords[,2] > 375000 & BG_thirtysecond_locations_sp2@coords[,2] <425000
sum(BG_remove3) / sum(!BG_remove2) # 94 % captured  - probably due to the closed river way we filled in

BG_thirtysecond_locations_sp2@coords[BG_remove3,] = t(apply(BG_thirtysecond_locations_sp2@coords[BG_remove3,],1,
                                                            function(x){x + cbind(0,-3000)}))
BG_remove4 = gContains(COAST_plotting, BG_thirtysecond_locations_sp2, byid = T)
sum(!BG_remove4) / length(BG_remove4) # 5.3% now lie on land - good reduction

BG_thirtysecond_locations_sp2$inwater3 = as.numeric(BG_remove3)
BG_thirtysecond_locations_sp2$inwater4 = as.numeric(BG_remove4)

ggplot() + gg(COAST_plotting) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater3))
ggplot() + gg(COAST_plotting) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater4))
# The 5% lie inside rivers we have filled in! Great :D

# Capture the final section of imputed paths across land and force them around land
BG_remove5 = BG_thirtysecond_locations_sp2@coords[,1] > 1050000 & BG_thirtysecond_locations_sp2@coords[,1] < 1075000 & !BG_remove4
BG_remove5 = BG_remove5 & BG_thirtysecond_locations_sp2@coords[,2] > 375000 & BG_thirtysecond_locations_sp2@coords[,2] <425000
sum(BG_remove5) / sum(!BG_remove4) # 81.3% captured  - probably due to the closed river way we filled in

BG_thirtysecond_locations_sp2@coords[BG_remove5,] = t(apply(BG_thirtysecond_locations_sp2@coords[BG_remove5,],1,
                                                            function(x){x + cbind(-3000,-4000)}))
BG_remove6 = gContains(COAST_plotting, BG_thirtysecond_locations_sp2, byid = T)
sum(!BG_remove6) / length(BG_remove6) # 0.9% now lie on land - as good as we'll get

BG_thirtysecond_locations_sp2$inwater5 = as.numeric(BG_remove5)
BG_thirtysecond_locations_sp2$inwater6 = as.numeric(BG_remove6)

ggplot() + gg(COAST_plotting) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater5))
ggplot() + gg(COAST_plotting) + gg(BG_thirtysecond_locations_sp2,aes(colour = inwater6))

# the remaining points on 'land', lie in side filled in rivers and estuaries. This is fine.

# remove intermediate objects
rm(list = c('BG_remove','BG_remove2','BG_remove3','BG_remove4','BG_remove5','BG_remove6',
            'BG_thirtysecond_locations_sp'))

# remove all observations outside of the time period of interest (i.e. May - October)
BG_thirtysecond_locations_sp3 = BG_thirtysecond_locations_sp2[
  which(BG_thirtysecond_locations_sp2$month < 11 &
          BG_thirtysecond_locations_sp2$month > 4),]

rm(BG_thirtysecond_locations_sp2)
#######################  END


#########
# Datetime stores the date and time for the BG predictions
# We have the TIME2 variable in the temp_X_verified_combined dataframes for the FIRST sightings 
# For EACH POD simply loop through the dates of the first sightings, removing all predicted BG positions 
# AFTER time of initial sighting. That is our search effort for that pod.

temp_J_verified_combined_dates2 = unique(temp_J_verified_combined$Date)
temp_K_verified_combined_dates2 = unique(temp_K_verified_combined$Date)
temp_L_verified_combined_dates2 = unique(temp_L_verified_combined$Date)

BG_thirtysecond_locations_sp_J = BG_thirtysecond_locations_sp3
BG_thirtysecond_locations_sp_K = BG_thirtysecond_locations_sp3
BG_thirtysecond_locations_sp_L = BG_thirtysecond_locations_sp3

BG_thirtysecond_locations_sp_J$KEEP = 1
for(i in temp_J_verified_combined_dates2)
{
  # extract all the tracklines for the date of sighting
  temp_BG_thirtysecond_J = BG_thirtysecond_locations_sp_J[which(BG_thirtysecond_locations_sp_J$Date == i),]
  # If no tracklines for that day. Simply go to next date as KEEP = 1.
  if(dim(temp_BG_thirtysecond_J)[1] > 0)
  {
    # Find the time of initial sighting
    time_of_first_sighting_J = temp_J_verified_combined$TIME2[temp_J_verified_combined$Date == i]
    # find the indices of the predicted boat locations AFTER the sighting
    ind_after_sighting_J = which(temp_BG_thirtysecond_J$Datetime > time_of_first_sighting_J)
    # set the KEEP variable to 1 for only those predicted locations
    BG_thirtysecond_locations_sp_J@data[which(BG_thirtysecond_locations_sp_J$Date == i),]$KEEP[ind_after_sighting_J] = 0
  }
}

BG_thirtysecond_locations_sp_K$KEEP = 1
for(i in temp_K_verified_combined_dates2)
{
  # extract all the tracklines for the date of sighting
  temp_BG_thirtysecond_K = BG_thirtysecond_locations_sp_K[which(BG_thirtysecond_locations_sp_K$Date == i),]
  # If no tracklines for that day. Simply go to next date as KEEP = 1.
  if(dim(temp_BG_thirtysecond_K)[1] > 0)
  {
    # Find the time of initial sighting
    time_of_first_sighting_K = temp_K_verified_combined$TIME2[temp_K_verified_combined$Date == i]
    # find the indices of the predicted boat locations AFTER the sighting
    ind_after_sighting_K = which(temp_BG_thirtysecond_K$Datetime > time_of_first_sighting_K)
    # set the KEEP variable to 1 for only those predicted locations
    BG_thirtysecond_locations_sp_K@data[which(BG_thirtysecond_locations_sp_K$Date == i),]$KEEP[ind_after_sighting_K] = 0
  }
}

BG_thirtysecond_locations_sp_L$KEEP = 1
for(i in temp_L_verified_combined_dates2)
{
  # extract all the tracklines for the date of sighting
  temp_BG_thirtysecond_L = BG_thirtysecond_locations_sp_L[which(BG_thirtysecond_locations_sp_L$Date == i),]
  # If no tracklines for that day. Simply go to next date as KEEP = 1.
  if(dim(temp_BG_thirtysecond_L)[1] > 0)
  {
    # Find the time of initial sighting
    time_of_first_sighting_L = temp_L_verified_combined$TIME2[temp_L_verified_combined$Date == i]
    # find the indices of the predicted boat locations AFTER the sighting
    ind_after_sighting_L = which(temp_BG_thirtysecond_L$Datetime > time_of_first_sighting_L)
    # set the KEEP variable to 1 for only those predicted locations
    BG_thirtysecond_locations_sp_L@data[which(BG_thirtysecond_locations_sp_L$Date == i),]$KEEP[ind_after_sighting_L] = 0
  }
}

dim(BG_thirtysecond_locations_sp3)[1] /(2*60) # 4980 boat hours
sum(BG_thirtysecond_locations_sp_J$KEEP) /(2*60) # 2341 boat hours
sum(BG_thirtysecond_locations_sp_K$KEEP) /(2*60) # 2997 boat hours
sum(BG_thirtysecond_locations_sp_L$KEEP) /(2*60) # 2701 boat hours

BG_thirtysecond_locations_sp_J = BG_thirtysecond_locations_sp_J[which(BG_thirtysecond_locations_sp_J$KEEP == 1),]
BG_thirtysecond_locations_sp_K = BG_thirtysecond_locations_sp_K[which(BG_thirtysecond_locations_sp_K$KEEP == 1),]
BG_thirtysecond_locations_sp_L = BG_thirtysecond_locations_sp_L[which(BG_thirtysecond_locations_sp_L$KEEP == 1),]

#########

# count how many points fall inside each dual mesh polygon for each pod
which_dmesh_pixel2_J = gWithin(BG_thirtysecond_locations_sp_J,dmesh, byid = T, returnDense = F)
which_dmesh_pixel2_J = unlist(which_dmesh_pixel2_J)
which_dmesh_pixel2_J = as.numeric(which_dmesh_pixel2_J)

which_dmesh_barrier_pixel2_J = gWithin(BG_thirtysecond_locations_sp_J,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_barrier_pixel2_J = unlist(which_dmesh_barrier_pixel2_J)
which_dmesh_barrier_pixel2_J = as.numeric(which_dmesh_barrier_pixel2_J)

SI_Brian_J = rep(0, times = mesh$n)
SI_Brian_J[as.numeric(names(table(which_dmesh_pixel2_J)))] = as.numeric(table(which_dmesh_pixel2_J))
SI_Brian_J[w == 0] = 0 # boats can't be on land
sum(SI_Brian_J>0) / length(SI_Brian_J) # 34% of pixels are visited by Brian

SI_Brian_barrier_J = rep(0, times = mesh_barrier$n)
SI_Brian_barrier_J[as.numeric(names(table(which_dmesh_barrier_pixel2_J)))] = as.numeric(table(which_dmesh_barrier_pixel2_J))
SI_Brian_barrier_J[w_barrier == 0] = 0 # boats can't be on land
sum(SI_Brian_barrier_J>0) / length(SI_Brian_barrier_J) # 41% of pixels are visited by Brian

rm(list = c('which_dmesh_pixel2_J','which_dmesh_barrier_pixel2_J'))

# Pod K
which_dmesh_pixel2_K = gWithin(BG_thirtysecond_locations_sp_K,dmesh, byid = T, returnDense = F)
which_dmesh_pixel2_K = unlist(which_dmesh_pixel2_K)
which_dmesh_pixel2_K = as.numeric(which_dmesh_pixel2_K)

which_dmesh_barrier_pixel2_K = gWithin(BG_thirtysecond_locations_sp_K,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_barrier_pixel2_K = unlist(which_dmesh_barrier_pixel2_K)
which_dmesh_barrier_pixel2_K = as.numeric(which_dmesh_barrier_pixel2_K)

SI_Brian_K = rep(0, times = mesh$n)
SI_Brian_K[as.numeric(names(table(which_dmesh_pixel2_K)))] = as.numeric(table(which_dmesh_pixel2_K))
SI_Brian_K[w == 0] = 0 # boats can't be on land
sum(SI_Brian_K>0) / length(SI_Brian_K) # 36% of pixels are visited by Brian

SI_Brian_barrier_K = rep(0, times = mesh_barrier$n)
SI_Brian_barrier_K[as.numeric(names(table(which_dmesh_barrier_pixel2_K)))] = as.numeric(table(which_dmesh_barrier_pixel2_K))
SI_Brian_barrier_K[w_barrier == 0] = 0 # boats can't be on land
sum(SI_Brian_barrier_K>0) / length(SI_Brian_barrier_K) # 43% of pixels are visited by Brian

rm(list = c('which_dmesh_pixel2_K','which_dmesh_barrier_pixel2_K'))

# Pod L

which_dmesh_pixel2_L = gWithin(BG_thirtysecond_locations_sp_L,dmesh, byid = T, returnDense = F)
which_dmesh_pixel2_L = unlist(which_dmesh_pixel2_L)
which_dmesh_pixel2_L = as.numeric(which_dmesh_pixel2_L)

which_dmesh_barrier_pixel2_L = gWithin(BG_thirtysecond_locations_sp_L,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_barrier_pixel2_L = unlist(which_dmesh_barrier_pixel2_L)
which_dmesh_barrier_pixel2_L = as.numeric(which_dmesh_barrier_pixel2_L)

SI_Brian_L = rep(0, times = mesh$n)
SI_Brian_L[as.numeric(names(table(which_dmesh_pixel2_L)))] = as.numeric(table(which_dmesh_pixel2_L))
SI_Brian_L[w == 0] = 0 # boats can't be on land
sum(SI_Brian_L>0) / length(SI_Brian_L) # 36% of pixels are visited by Brian

SI_Brian_barrier_L = rep(0, times = mesh_barrier$n)
SI_Brian_barrier_L[as.numeric(names(table(which_dmesh_barrier_pixel2_L)))] = as.numeric(table(which_dmesh_barrier_pixel2_L))
SI_Brian_barrier_L[w_barrier == 0] = 0 # boats can't be on land
sum(SI_Brian_barrier_L>0) / length(SI_Brian_barrier_L) # 43% of pixels are visited by Brian

rm(list = c('which_dmesh_pixel2_L','which_dmesh_barrier_pixel2_L'))

# repeat for each month
SI_Brian_monthly_J = matrix(0, nrow = mesh$n, ncol = 6)
SI_Brian_barrier_monthly_J = matrix(0, nrow = mesh_barrier$n, ncol = 6)
count = 1
for( i in 5:10)
{
  ind_month = which(BG_thirtysecond_locations_sp_J$month == i)
  
  # count how many points fall inside each dual mesh polygon each month
  which_dmesh_pixel2_month_J = gWithin(BG_thirtysecond_locations_sp_J[ind_month,],dmesh, byid = T, returnDense = F)
  which_dmesh_pixel2_month_J = unlist(which_dmesh_pixel2_month_J)
  which_dmesh_pixel2_month_J = as.numeric(which_dmesh_pixel2_month_J)
  
  which_dmesh_barrier_pixel2_month_J = gWithin(BG_thirtysecond_locations_sp_J[ind_month,],dmesh_barrier, byid = T, returnDense = F)
  which_dmesh_barrier_pixel2_month_J = unlist(which_dmesh_barrier_pixel2_month_J)
  which_dmesh_barrier_pixel2_month_J = as.numeric(which_dmesh_barrier_pixel2_month_J)
  
  SI_Brian_monthly_J[as.numeric(names(table(which_dmesh_pixel2_month_J))),count] = as.numeric(table(which_dmesh_pixel2_month_J))
  SI_Brian_monthly_J[w == 0,count] = 0 # boats can't be on land
  
  SI_Brian_barrier_monthly_J[as.numeric(names(table(which_dmesh_barrier_pixel2_month_J))),count] = as.numeric(table(which_dmesh_barrier_pixel2_month_J))
  SI_Brian_barrier_monthly_J[w_barrier == 0,count] = 0 # boats can't be on land
  
  print(paste('iteration',count,'out of 6'))
  count = count+1
}

# Pod K
SI_Brian_monthly_K = matrix(0, nrow = mesh$n, ncol = 6)
SI_Brian_barrier_monthly_K = matrix(0, nrow = mesh_barrier$n, ncol = 6)
count = 1
for( i in 5:10)
{
  ind_month = which(BG_thirtysecond_locations_sp_K$month == i)
  
  # count how many points fall inside each dual mesh polygon each month
  which_dmesh_pixel2_month_K = gWithin(BG_thirtysecond_locations_sp_K[ind_month,],dmesh, byid = T, returnDense = F)
  which_dmesh_pixel2_month_K = unlist(which_dmesh_pixel2_month_K)
  which_dmesh_pixel2_month_K = as.numeric(which_dmesh_pixel2_month_K)
  
  which_dmesh_barrier_pixel2_month_K = gWithin(BG_thirtysecond_locations_sp_K[ind_month,],dmesh_barrier, byid = T, returnDense = F)
  which_dmesh_barrier_pixel2_month_K = unlist(which_dmesh_barrier_pixel2_month_K)
  which_dmesh_barrier_pixel2_month_K = as.numeric(which_dmesh_barrier_pixel2_month_K)
  
  SI_Brian_monthly_K[as.numeric(names(table(which_dmesh_pixel2_month_K))),count] = as.numeric(table(which_dmesh_pixel2_month_K))
  SI_Brian_monthly_K[w == 0,count] = 0 # boats can't be on land
  
  SI_Brian_barrier_monthly_K[as.numeric(names(table(which_dmesh_barrier_pixel2_month_K))),count] = as.numeric(table(which_dmesh_barrier_pixel2_month_K))
  SI_Brian_barrier_monthly_K[w_barrier == 0,count] = 0 # boats can't be on land
  
  print(paste('iteration',count,'out of 6'))
  count = count+1
}

# Pod L
SI_Brian_monthly_L = matrix(0, nrow = mesh$n, ncol = 6)
SI_Brian_barrier_monthly_L = matrix(0, nrow = mesh_barrier$n, ncol = 6)
count = 1
for( i in 5:10)
{
  ind_month = which(BG_thirtysecond_locations_sp_L$month == i)
  
  # count how many points fall inside each dual mesh polygon each month
  which_dmesh_pixel2_month_L = gWithin(BG_thirtysecond_locations_sp_L[ind_month,],dmesh, byid = T, returnDense = F)
  which_dmesh_pixel2_month_L = unlist(which_dmesh_pixel2_month_L)
  which_dmesh_pixel2_month_L = as.numeric(which_dmesh_pixel2_month_L)
  
  which_dmesh_barrier_pixel2_month_L = gWithin(BG_thirtysecond_locations_sp_L[ind_month,],dmesh_barrier, byid = T, returnDense = F)
  which_dmesh_barrier_pixel2_month_L = unlist(which_dmesh_barrier_pixel2_month_L)
  which_dmesh_barrier_pixel2_month_L = as.numeric(which_dmesh_barrier_pixel2_month_L)
  
  SI_Brian_monthly_L[as.numeric(names(table(which_dmesh_pixel2_month_L))),count] = as.numeric(table(which_dmesh_pixel2_month_L))
  SI_Brian_monthly_L[w == 0,count] = 0 # boats can't be on land
  
  SI_Brian_barrier_monthly_L[as.numeric(names(table(which_dmesh_barrier_pixel2_month_L))),count] = as.numeric(table(which_dmesh_barrier_pixel2_month_L))
  SI_Brian_barrier_monthly_L[w_barrier == 0,count] = 0 # boats can't be on land
  
  print(paste('iteration',count,'out of 6'))
  count = count+1
}

# plot Brian's search effort
pixels_plotting =  pixels(mesh_barrier_transformed, mask = COAST_transformed, nx = 300, ny = 300)

plot.field(SI_Brian_barrier_J, poly = COAST_transformed, mesh = mesh_barrier_transformed, 
           pixels = pixels_plotting, corr = FALSE)
plot.field(log(SI_Brian_barrier_J+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
           pixels = pixels_plotting, corr = FALSE)

plot.field(SI_Brian_J, poly = COAST_transformed, mesh = mesh_transformed, 
           pixels = pixels_plotting, corr = FALSE)
plot.field(log(SI_Brian_J+1), poly = COAST_transformed, mesh = mesh_transformed, 
           pixels = pixels_plotting, corr = FALSE)

months_vec = c('May','June','July','August','September','October')
# Plot by month
 multiplot(   
  plot.field(log(SI_Brian_barrier_monthly_J[,1]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[1])),
  plot.field(log(SI_Brian_barrier_monthly_J[,2]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[2])),
  plot.field(log(SI_Brian_barrier_monthly_J[,3]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[3])), 
  plot.field(log(SI_Brian_barrier_monthly_J[,4]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[4])), 
  plot.field(log(SI_Brian_barrier_monthly_J[,5]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[5])), 
  plot.field(log(SI_Brian_barrier_monthly_J[,6]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
             pixels = pixels_plotting, corr = FALSE) + ggtitle(paste('Brian Search Effort for Pod J in', months_vec[6])),
layout = matrix(c(1,2,3,4,5,6),nrow = 3, byrow = T))

# each point is worth 30 seconds of search time. Scale into hours
SI_Brian_J = (SI_Brian_J/2)/60
SI_Brian_barrier_J = (SI_Brian_barrier_J/2)/60
SI_Brian_monthly_J = (SI_Brian_monthly_J/2)/60
SI_Brian_barrier_monthly_J = (SI_Brian_barrier_monthly_J/2)/60
colSums(SI_Brian_monthly_J)

SI_Brian_K = (SI_Brian_K/2)/60
SI_Brian_barrier_K = (SI_Brian_barrier_K/2)/60
SI_Brian_monthly_K = (SI_Brian_monthly_K/2)/60
SI_Brian_barrier_monthly_K = (SI_Brian_barrier_monthly_K/2)/60
colSums(SI_Brian_monthly_K)

SI_Brian_L = (SI_Brian_L/2)/60
SI_Brian_barrier_L = (SI_Brian_barrier_L/2)/60
SI_Brian_monthly_L = (SI_Brian_monthly_L/2)/60
SI_Brian_barrier_monthly_L = (SI_Brian_barrier_monthly_L/2)/60
colSums(SI_Brian_monthly_L)
sum(SI_Brian_L)
sum(colSums(SI_Brian_monthly_L)) # woo
#save.image('SI_WS_dup_removed2.RData')
load('SI_WS_dup_removed2.RData')

# Load one sample of the search effort for model selection
# Computing the mean using the sampling intensity MC compiler R script
# leads to a non-smooth search effort if we take the mean
load('WW_SI_compiled_remove_dup.RData')

# Some NAs exist. Fill them in by using the dinterp function that uses local median 
sum(SI_WW_complete_J_mean,na.rm=T)  
sum(mean_total_boat_hours_per_period_per_port_J)
# not close in magnitude - due to reasons below

pid <- sapply(slot(dmesh, "polygons"), function(x) slot(x, "ID")) 
SI_dmesh = SpatialPolygonsDataFrame(dmesh,
                                    data.frame(cbind(SI_WW_complete_J_mean, #SI_WW_complete_J_var,
                                                     SI_WW_complete_K_mean, #SI_WW_complete_K_var,
                                                     SI_WW_complete_L_mean, #SI_WW_complete_L_var,
                                                     SI_WW_complete_monthly_J_mean, #SI_WW_complete_monthly_J_var,
                                                     SI_WW_complete_monthly_K_mean, #SI_WW_complete_monthly_K_var,
                                                     SI_WW_complete_monthly_L_mean),#, SI_WW_complete_monthly_L_var),
                                               row.names = pid))
pid <- sapply(slot(dmesh_barrier, "polygons"), function(x) slot(x, "ID")) 
SI_dmesh_barrier = SpatialPolygonsDataFrame(dmesh_barrier,
                                            data.frame(cbind(SI_WW_complete_J_mean_barrier, #SI_WW_complete_J_var_barrier@data,
                                                             SI_WW_complete_K_mean_barrier, #SI_WW_complete_K_var_barrier@data,
                                                             SI_WW_complete_L_mean_barrier, #SI_WW_complete_L_var_barrier@data,
                                                             SI_WW_complete_monthly_J_mean_barrier, #SI_WW_complete_monthly_J_var_barrier@data,
                                                             SI_WW_complete_monthly_K_mean_barrier, #SI_WW_complete_monthly_K_var_barrier@data,
                                                             SI_WW_complete_monthly_L_mean_barrier),#, SI_WW_complete_monthly_L_var_barrier@data),
                                                       row.names = pid))

# WARNING SLOW
SI_dmesh = f.interp.dmesh(dmesh, SI_dmesh,dmesh_dist)$weighted_mean@data
SI_dmesh_barrier = f.interp.dmesh(dmesh_barrier, SI_dmesh_barrier,dmesh_barrier_dist)$weighted_mean@data

names_SI_dmesh = c('SI_WW_complete_J_mean', #'SI_WW_complete_J_var',
                   'SI_WW_complete_K_mean', #'SI_WW_complete_K_var',
                   'SI_WW_complete_L_mean', #'SI_WW_complete_L_var',
                   rep('SI_WW_complete_monthly_J_mean',6), 
                   #rep('SI_WW_complete_monthly_J_var',6),
                   rep('SI_WW_complete_monthly_K_mean',6), 
                   #rep('SI_WW_complete_monthly_K_var',6),
                   rep('SI_WW_complete_monthly_L_mean',6))#,
#rep('SI_WW_complete_monthly_L_var',6))

names_SI_dmesh_barrier = c('SI_WW_complete_J_mean_barrier', #'SI_WW_complete_J_var_barrier',
                           'SI_WW_complete_K_mean_barrier', #'SI_WW_complete_K_var_barrier',
                           'SI_WW_complete_L_mean_barrier', #'SI_WW_complete_L_var_barrier',
                           rep('SI_WW_complete_monthly_J_mean_barrier',6), 
                           #rep('SI_WW_complete_monthly_J_var_barrier',6),
                           rep('SI_WW_complete_monthly_K_mean_barrier',6), 
                           #rep('SI_WW_complete_monthly_K_var_barrier',6),
                           rep('SI_WW_complete_monthly_L_mean_barrier',6))#,
#rep('SI_WW_complete_monthly_L_var_barrier',6))

sum(SI_WW_complete_J_mean,na.rm=T)  
sum(mean_total_boat_hours_per_period_per_port_J)
sum(SI_dmesh$SI_WW_complete_J_mean,na.rm=T)
sum(apply(SI_dmesh, 2 ,FUN=function(x){x*(w/mean(w))}),na.rm=T)

#save.image('Temp_WS_dup_removed.RData')
load('Temp_WS_dup_removed.RData')

# Sampling_Intesity_MC_compiler does not in fact integrate the dmesh pixels over the SI.
# It computes the weighted average and thus we need to scale by the area of the pixels
SI_dmesh = apply(SI_dmesh, 2 ,FUN=function(x){x*(w/mean(w))})
SI_dmesh_barrier = apply(SI_dmesh_barrier, 2 ,FUN=function(x){x*(w_barrier/mean(w_barrier))})

for(i in unique(names_SI_dmesh))
{
  ind_name = which(names_SI_dmesh == i)
  assign(i, SI_dmesh[,ind_name])
}
for(i in unique(names_SI_dmesh_barrier))
{
  ind_name = which(names_SI_dmesh_barrier == i)
  assign(i, SI_dmesh_barrier[,ind_name])
}

rm(SI_dmesh, SI_dmesh_barrier)
gc()

summary(SI_WW_complete_monthly_J_mean)
summary(SI_WW_complete_J_mean)

# save temp file
#save.image('temp_checking.RData')

# Now set the WW search effort to zero on land
sum(SI_WW_complete_J_mean)
sum(SI_WW_complete_monthly_J_mean)
sum(mean_total_boat_hours_per_period_per_port_J)
# not close in magnitude still
SI_WW_complete_J_mean[w==0] = 0
SI_WW_complete_K_mean[w==0] = 0
SI_WW_complete_L_mean[w==0] = 0

SI_WW_complete_monthly_J_mean[w==0,] = 0
SI_WW_complete_monthly_K_mean[w==0,] = 0
SI_WW_complete_monthly_L_mean[w==0,] = 0

SI_WW_complete_J_mean_barrier[w_barrier==0] = 0
SI_WW_complete_K_mean_barrier[w_barrier==0] = 0
SI_WW_complete_L_mean_barrier[w_barrier==0] = 0

SI_WW_complete_monthly_J_mean_barrier[w_barrier==0,] = 0
SI_WW_complete_monthly_K_mean_barrier[w_barrier==0,] = 0
SI_WW_complete_monthly_L_mean_barrier[w_barrier==0,] = 0

sum(SI_WW_complete_J_mean)#,na.rm=T)  # as expected not much of a decrease

# rescale to equal the mean
rescale_fun = function(input, scale)
{
  if(length(scale)>1){# we should have a matrix input
    if(dim(input)[2] != length(scale))
    {
      stop('The length of the scale argument should match the number of cols of input')
    }
    for(i in 1:length(scale))
    {
      input[,i] = (input[,i] * scale[i]) / sum(input[,i], na.rm=T)
    }
  }
  if(length(scale)==1){# we should have a vector input
    if(!is.vector(input))
    {
      stop('The input should be a vector')
    }
    input = (input * scale) / sum(input, na.rm=T)}
  return(input)
}

SI_WW_complete_J_mean = rescale_fun(SI_WW_complete_J_mean, sum(mean_total_boat_hours_per_period_per_port_J) )
SI_WW_complete_K_mean = rescale_fun(SI_WW_complete_K_mean, sum(mean_total_boat_hours_per_period_per_port_K) )
SI_WW_complete_L_mean = rescale_fun(SI_WW_complete_L_mean, sum(mean_total_boat_hours_per_period_per_port_L) )

SI_WW_complete_monthly_J_mean = rescale_fun(SI_WW_complete_monthly_J_mean, colSums(mean_total_boat_hours_per_period_per_port_J))
SI_WW_complete_monthly_K_mean = rescale_fun(SI_WW_complete_monthly_K_mean, colSums(mean_total_boat_hours_per_period_per_port_K))
SI_WW_complete_monthly_L_mean = rescale_fun(SI_WW_complete_monthly_L_mean, colSums(mean_total_boat_hours_per_period_per_port_L))

SI_WW_complete_J_mean_barrier = rescale_fun(SI_WW_complete_J_mean_barrier, sum(mean_total_boat_hours_per_period_per_port_J) )
SI_WW_complete_K_mean_barrier = rescale_fun(SI_WW_complete_K_mean_barrier, sum(mean_total_boat_hours_per_period_per_port_K) )
SI_WW_complete_L_mean_barrier = rescale_fun(SI_WW_complete_L_mean_barrier, sum(mean_total_boat_hours_per_period_per_port_L) )

SI_WW_complete_monthly_J_mean_barrier = rescale_fun(SI_WW_complete_monthly_J_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_J))
SI_WW_complete_monthly_K_mean_barrier = rescale_fun(SI_WW_complete_monthly_K_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_K))
SI_WW_complete_monthly_L_mean_barrier = rescale_fun(SI_WW_complete_monthly_L_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_L))

# How does this compare with Brian's search effort in the period?
colSums(SI_Brian_monthly_J)
colSums(SI_WW_complete_monthly_J_mean)

colSums(SI_WW_complete_monthly_J_mean) / colSums(SI_Brian_monthly_J)
colSums(SI_WW_complete_monthly_K_mean) / colSums(SI_Brian_monthly_K)
colSums(SI_WW_complete_monthly_L_mean) / colSums(SI_Brian_monthly_L)
# WW companies provide between 31 - 60 times the search effort per month

# check with a plot
May_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,1]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1)) + ggtitle('May WW log Search Effort J pod')
Jun_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,2]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1)) + ggtitle('June WW log Search Effort J pod')
Jul_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,3]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1))+ ggtitle('July WW log Search Effort J pod')
Aug_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,4]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1)) + ggtitle('August WW log Search Effort J pod')
Sep_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,5]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1)) + ggtitle('September WW log Search Effort J pod')
Oct_mean_WW_SI_plot =  plot.field(log(SI_WW_complete_monthly_J_mean[,6]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean)+1)) + ggtitle('October WW log Search Effort J pod')
multiplot(May_mean_WW_SI_plot, Jun_mean_WW_SI_plot, Jul_mean_WW_SI_plot, Aug_mean_WW_SI_plot,
          Sep_mean_WW_SI_plot, Oct_mean_WW_SI_plot, 
          layout = matrix(c(1,2,3,4,5,6),nrow = 3, byrow = T))

May_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,1], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean)) + ggtitle('May WW Search Effort J pod')
Jun_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,2], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean)) + ggtitle('June WW Search Effort J pod')
Jul_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,3], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean))+ ggtitle('July WW Search Effort J pod')
Aug_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,4], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean)) + ggtitle('August WW Search Effort J pod')
Sep_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,5], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean)) + ggtitle('September WW Search Effort J pod')
Oct_mean_WW_SI_plot =  plot.field(SI_WW_complete_monthly_J_mean[,6], poly = COAST_transformed, mesh = mesh_transformed, 
                                  pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_WW_complete_monthly_J_mean)) + ggtitle('October WW Search Effort J pod')
multiplot(May_mean_WW_SI_plot, Jun_mean_WW_SI_plot, Jul_mean_WW_SI_plot, Aug_mean_WW_SI_plot,
          Sep_mean_WW_SI_plot, Oct_mean_WW_SI_plot, 
          layout = matrix(c(1,2,3,4,5,6),nrow = 3, byrow = T))

rm(May_mean_WW_SI_plot, May_sd_WW_SI_plot, Jun_mean_WW_SI_plot, Jun_sd_WW_SI_plot, Jul_mean_WW_SI_plot, Jul_sd_WW_SI_plot, Aug_mean_WW_SI_plot, Aug_sd_WW_SI_plot, 
   Sep_mean_WW_SI_plot, Sep_sd_WW_SI_plot, Oct_mean_WW_SI_plot, Oct_sd_WW_SI_plot )


# Combine the search efforts together. These are on boat hour scales and have been adjusted by area already
# BG has had counts, thus bigger area -> bigger counts. WW has had effort density layer integrated over pixels.

ggplot(data = data.frame(y = c(colSums(SI_WW_complete_monthly_J_mean), colSums(SI_WW_complete_monthly_K_mean),colSums(SI_WW_complete_monthly_L_mean)),
                         x = factor(rep(months_vec, times = 3), levels = months_vec),
                         pod = factor(rep(c('J','K','L'),each = 6))),
       aes(y=y,x=x,colour = pod, group = pod)) + geom_point() + geom_line() +
  ggtitle('Monthly Total Corrected Whale Watch Search Effort Per Pod') + ylab('Boat Hours')

# Why does the total search effort decrease in the Summer months for J pod?
# J pod is found most frequently and
ggplot(data = data.frame(x = c(hour(J_verified$DateTimes), hour(K_verified$DateTimes),hour(L_verified$DateTimes)),
                         month = factor(c(J_verified$MONTH, K_verified$MONTH, L_verified$MONTH)),
                         pod = factor(rep(c('J','K','L'),times = c(dim(J_verified)[1], dim(K_verified)[1], dim(L_verified)[1])))),
       aes(x=x,colour = pod, group = pod, fill = pod)) + geom_histogram(breaks=seq(0,23,by=1),alpha = 0.2) + facet_wrap(~ month) +
  ggtitle('Number of initial sightings by hour of day and by month',subtitle ='Split per Pod') + xlab('Hour of the day') +
  geom_vline(xintercept = 9)
# The large number of sightings made in the summer months (especially for pod J), combined with their
# earlier sighting times (they are found earlier in the day in the summer) lead to the search effort being lower.

#### ERROR IN CODE - NEED TO MAP SIGHTINGS COVARIATES AGAIN
total_sightings_J = temp_J_verified_combined[,c('YEAR','MONTH')] #rbind(J_verified[,c('YEAR','MONTH','J','K','L')],J_verified_Brian[,c('YEAR','MONTH','J','K','L')])
total_sightings_K = temp_K_verified_combined[,c('YEAR','MONTH')] #rbind(K_verified[,c('YEAR','MONTH','J','K','L')],K_verified_Brian[,c('YEAR','MONTH','J','K','L')])
total_sightings_L = temp_L_verified_combined[,c('YEAR','MONTH')] #rbind(L_verified[,c('YEAR','MONTH','J','K','L')],L_verified_Brian[,c('YEAR','MONTH','J','K','L')])

# Combine the search efforts together. These are on boat hour scales and have been adjusted by area already
# BG has had counts, thus bigger area -> bigger counts. WW has had effort density layer integrated over pixels.
SI_complete_monthly_J_mean = SI_WW_complete_monthly_J_mean + SI_Brian_monthly_J
SI_complete_monthly_K_mean = SI_WW_complete_monthly_K_mean + SI_Brian_monthly_K
SI_complete_monthly_L_mean = SI_WW_complete_monthly_L_mean + SI_Brian_monthly_L

SI_complete_monthly_J_mean_barrier = SI_WW_complete_monthly_J_mean_barrier + SI_Brian_barrier_monthly_J
SI_complete_monthly_K_mean_barrier = SI_WW_complete_monthly_K_mean_barrier + SI_Brian_barrier_monthly_K
SI_complete_monthly_L_mean_barrier = SI_WW_complete_monthly_L_mean_barrier + SI_Brian_barrier_monthly_L

SI_complete_J_mean = SI_WW_complete_J_mean + SI_Brian_J
SI_complete_K_mean = SI_WW_complete_K_mean + SI_Brian_K
SI_complete_L_mean = SI_WW_complete_L_mean + SI_Brian_L

SI_complete_J_mean_barrier = SI_WW_complete_J_mean_barrier + SI_Brian_barrier_J
SI_complete_K_mean_barrier = SI_WW_complete_K_mean_barrier + SI_Brian_barrier_K
SI_complete_L_mean_barrier = SI_WW_complete_L_mean_barrier + SI_Brian_barrier_L

# Repeat for the total (combined effort)
May_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,1], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean)) + ggtitle('May  Search Effort J pod')
Jun_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,2], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean)) + ggtitle('June  Search Effort J pod')
Jul_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,3], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean))+ ggtitle('July  Search Effort J pod')
Aug_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,4], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean)) + ggtitle('August  Search Effort J pod')
Sep_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,5], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean)) + ggtitle('September  Search Effort J pod')
Oct_mean_SI_plot =  plot.field(SI_complete_monthly_J_mean[,6], poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(as.matrix(SI_complete_monthly_J_mean)) + ggtitle('October  Search Effort J pod')
multiplot(May_mean_SI_plot, Jun_mean_SI_plot, Jul_mean_SI_plot, Aug_mean_SI_plot,
          Sep_mean_SI_plot, Oct_mean_SI_plot, 
          layout = matrix(c(1,2,3,4,5,6),nrow = 3, byrow = T))

May_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,1]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + ggtitle('May Total log Search Effort J pod')
Jun_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,2]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + ggtitle('June Total log Search Effort J pod')
Jul_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,3]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1))+ ggtitle('July Total log Search Effort J pod')
Aug_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,4]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + ggtitle('August Total log Search Effort J pod')
Sep_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,5]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + ggtitle('September Total log Search Effort J pod')
Oct_mean_SI_plot =  plot.field(log(SI_complete_monthly_J_mean[,6]+1), poly = COAST_transformed, mesh = mesh_transformed, 
                               pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + ggtitle('October Total log Search Effort J pod')
multiplot(May_mean_SI_plot, Jun_mean_SI_plot, Jul_mean_SI_plot, Aug_mean_SI_plot,
          Sep_mean_SI_plot, Oct_mean_SI_plot, 
          layout = matrix(c(1,2,3,4,5,6),nrow = 3, byrow = T))

rm(May_mean_SI_plot, Jun_mean_SI_plot, Jul_mean_SI_plot,  Aug_mean_SI_plot,
   Sep_mean_SI_plot, Oct_mean_SI_plot )

# Format the data in the manner needed for INLA

no_T = 6

total_sightings_J$POD_INLA = 1
total_sightings_K$POD_INLA = 2
total_sightings_L$POD_INLA = 3

total_sightings_J@data$MONTH_INLA = total_sightings_J@data$MONTH - 4
total_sightings_K@data$MONTH_INLA = total_sightings_K@data$MONTH - 4
total_sightings_L@data$MONTH_INLA = total_sightings_L@data$MONTH - 4

# check all observations lie between May - Oct
summary(total_sightings_J@data$MONTH_INLA)
summary(total_sightings_K@data$MONTH_INLA)
summary(total_sightings_L@data$MONTH_INLA)
total_sightings_L = total_sightings_L[which(total_sightings_L@data$MONTH_INLA<7),]
# Use MONTH_INLA to map covariates 

 # Obtain vessel density values at observation locations
 vessel_obs_J = f.interp(sp_points = total_sightings_J, 
                        sp_grid = vessel_density)
 vessel_obs_K = f.interp(sp_points = total_sightings_K, 
                         sp_grid = vessel_density)
 vessel_obs_L = f.interp(sp_points = total_sightings_L, 
                         sp_grid = vessel_density)

 par(mfrow=c(1,1))
 # plot the vessel density of the area vs those experienced by the whale.
 # remember to weight by the area of the pixels.
 ggplot(data = data.frame(x = vessel_dmesh$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),1],
                          weight = w[gContains(COAST_simp, dmesh, byid=T)]/sum(w[gContains(COAST_simp, dmesh, byid=T)])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = vessel_obs_J$band1, weight = 1/length(vessel_obs_J$band1)),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   xlim(c(0, 250)) + ggtitle('Histogram of vessel density within AOI and at SRKW J pod locations' ) + xlab('Vessel Density') + guides(colour = F)
 
 # 2 Chloro_2
 chloro_obs_J = f.interp(sp_points = total_sightings_J, 
                       sp_grid = chloro_monthly)
 chloro_obs_K = f.interp(sp_points = total_sightings_K, 
                        sp_grid = chloro_monthly)
 chloro_obs_L = f.interp(sp_points = total_sightings_L, 
                          sp_grid = chloro_monthly)
 # plot the chloro2 value of the area vs those experienced by the whale.
 # Remember to subset the pixels ONLY CONTAINED WITHIN COAST_SIMP to avoid land-based values.
 # weight by the area of pixels
 
 ggplot(data = data.frame(x = as.numeric(as.matrix(chloro_dmesh_month$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),])),
                          weight = rep(w[gContains(COAST_simp, dmesh, byid=T)]/(dim(chloro_dmesh_month$weighted_mean@data)[2]*sum(w[gContains(COAST_simp, dmesh, byid=T)])), times = dim(chloro_dmesh_month$weighted_mean@data)[2])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = as.numeric(as.matrix(chloro_obs_J)), weight = 1/length(as.numeric(as.matrix(chloro_obs_J)))),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('Histogram of chloro 2 density within AOI and at SRKW J pod locations' ) + xlab('Chloro 2 Density') + guides(colour = F)
 ggplot(data = data.frame(x = log(as.numeric(as.matrix(chloro_dmesh_month$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),]))),
                          weight = rep(w[gContains(COAST_simp, dmesh, byid=T)]/(dim(chloro_dmesh_month$weighted_mean@data)[2]*sum(w[gContains(COAST_simp, dmesh, byid=T)])), times = dim(chloro_dmesh_month$weighted_mean@data)[2])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = log(as.numeric(as.matrix(chloro_obs_J))), weight = 1/length(as.numeric(as.matrix(chloro_obs_J)))),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('Histogram of log chloro 2 density within AOI and at SRKW J pod locations' ) + xlab('Chloro 2 Density') + guides(colour = F)
 
 # 4 SST_2
 SST_obs_J = f.interp(sp_points = total_sightings_J, 
                      sp_grid = SST_monthly)
 SST_obs_K = f.interp(sp_points = total_sightings_K, 
                        sp_grid = SST_monthly)
 SST_obs_L = f.interp(sp_points = total_sightings_L, 
                        sp_grid = SST_monthly)

# # plot the SST_2 value of the area vs those experienced by the whale.
# # Remember to subset the pixels ONLY CONTAINED WITHIN COAST_SIMP to avoid land-based values.
# # weight by the area of pixels
# 
 ggplot(data = data.frame(x = as.numeric(as.matrix(SST_dmesh_month$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),])),
                          weight = rep(w[gContains(COAST_simp, dmesh, byid=T)]/(dim(SST_dmesh_month$weighted_mean@data)[2]*sum(w[gContains(COAST_simp, dmesh, byid=T)])), times = dim(SST_dmesh_month$weighted_mean@data)[2])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = as.numeric(as.matrix(SST_obs_J)), weight = 1/length(as.numeric(as.matrix(SST_obs_J)))),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('SST_2 density within AOI and at SRKW J pod locations' ) + xlab('SST_2 Density') + guides(colour = F)
 ggplot(data = data.frame(x = log(as.numeric(as.matrix(SST_dmesh_month$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),]))),
                          weight = rep(w[gContains(COAST_simp, dmesh, byid=T)]/(dim(SST_dmesh_month$weighted_mean@data)[2]*sum(w[gContains(COAST_simp, dmesh, byid=T)])), times = dim(SST_dmesh_month$weighted_mean@data)[2])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = log(as.numeric(as.matrix(SST_obs_J))), weight = 1/length(as.numeric(as.matrix(SST_obs_J)))),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('Log SST_2 density within AOI and at SRKW J pod locations' ) + xlab('SST_2 Density') + guides(colour = F)

# # 4 topo
 topo_obs_J = f.interp(sp_points = total_sightings_J, 
                      sp_grid = topo)
 topo_obs_K = f.interp(sp_points = total_sightings_K, 
                     sp_grid = topo)
 topo_obs_L = f.interp(sp_points = total_sightings_L, 
                     sp_grid = topo)

# # plot the topo value of the area vs those experienced by the whale.
# # Remember to subset the pixels ONLY CONTAINED WITHIN COAST_SIMP to avoid land-based values.
# # weight by area of pixels
 ggplot(data = data.frame(x = topo_dmesh$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),1],
                          weight = w[gContains(COAST_simp, dmesh, byid=T)]/sum(w[gContains(COAST_simp, dmesh, byid=T)])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = topo_obs_J$V1, weight = 1/length(topo_obs_J$V1)),
                                                                   aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('Topo density within AOI and at SRKW J pod locations' ) + xlab('Topo Density') + guides(colour = F)
 
 ggplot(data = data.frame(x = log(abs(topo_dmesh$weighted_mean@data[gContains(COAST_simp, dmesh, byid=T),1]-51)),
                          weight = w[gContains(COAST_simp, dmesh, byid=T)]/sum(w[gContains(COAST_simp, dmesh, byid=T)])),
        aes(x=x, weight = weight)) + geom_density() + geom_density(data = data.frame(x = log(abs(topo_obs_J$V1-51)), weight = 1/length(topo_obs_J$V1)),
                                                                  aes(x=x, weight = weight, colour = 'red')) +
   ggtitle('Log depth density within AOI and at SRKW J pod locations' ) + xlab('Depth Density') + guides(colour = F)
 

######################################################
# MAP THE POINTS OONTO THE CORRECT MONTH'S COVARS
######################################################

colnames(chloro_obs_J) = c('1','2','3','4','5','6')
colnames(chloro_obs_K) = c('1','2','3','4','5','6')
colnames(chloro_obs_L) = c('1','2','3','4','5','6')

colnames(SST_obs_J) = c('1','2','3','4','5','6')
colnames(SST_obs_K) = c('1','2','3','4','5','6')
colnames(SST_obs_L) = c('1','2','3','4','5','6')

# create mappings to the covariates
covar_map_J = sapply(total_sightings_J@data$MONTH_INLA,
                   FUN = function(x){which(colnames(chloro_obs_J) == x)},
                   simplify = T)
covar_map_J = unlist(covar_map_J)

covar_map_K = sapply(total_sightings_K@data$MONTH_INLA,
                     FUN = function(x){which(colnames(chloro_obs_K) == x)},
                     simplify = T)
covar_map_K = unlist(covar_map_K)

covar_map_L = sapply(total_sightings_L@data$MONTH_INLA,
                     FUN = function(x){which(colnames(chloro_obs_L) == x)},
                     simplify = T)
covar_map_L = unlist(covar_map_L)
# 
 total_sightings_J@data$chloro = chloro_obs_J[cbind(seq(1,dim(chloro_obs_J)[1]),as.numeric(covar_map_J))]
 total_sightings_J@data$SST = SST_obs_J[cbind(seq(1,dim(chloro_obs_J)[1]),as.numeric(covar_map_J))]
 total_sightings_J@data$topo = as.numeric(topo_obs_J$V1)
 total_sightings_J@data$vessel = vessel_obs_J$band1

 total_sightings_K@data$chloro = chloro_obs_K[cbind(seq(1,dim(chloro_obs_K)[1]),as.numeric(covar_map_K))]
 total_sightings_K@data$SST = SST_obs_K[cbind(seq(1,dim(chloro_obs_K)[1]),as.numeric(covar_map_K))]
 total_sightings_K@data$topo = as.numeric(topo_obs_K$V1)
 total_sightings_K@data$vessel = vessel_obs_K$band1
 
 total_sightings_L@data$chloro = chloro_obs_L[cbind(seq(1,dim(chloro_obs_L)[1]),as.numeric(covar_map_L))]
 total_sightings_L@data$SST = SST_obs_L[cbind(seq(1,dim(chloro_obs_L)[1]),as.numeric(covar_map_L))]
 total_sightings_L@data$topo = as.numeric(topo_obs_L$V1)
 total_sightings_L@data$vessel = vessel_obs_L$band1
#
 
 ## Corrected total_sightings
 #save.image('SI_WS2_dup_removed.RData')

# create the projector matrix for each pod
A_proj_J = inla.spde.make.A(mesh, coordinates(total_sightings_J),
                          group = total_sightings_J$MONTH_INLA,
                          n.group = no_T,
                          repl = total_sightings_J$POD_INLA,
                          n.repl = 3)
A_proj_barrier_J = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_J),
                                  group = total_sightings_J$MONTH_INLA,
                                  n.group = no_T,
                                  repl = total_sightings_J$POD_INLA,
                                  n.repl = 3)
imat_J = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                          group = rep(1:no_T,each = mesh$n),
                          n.group = no_T,
                          repl = rep(1, times = mesh$n*no_T),
                          n.repl = 3)
imat_barrier_J = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                          group = rep(1:no_T,each = mesh_barrier$n),
                          n.group = no_T,
                          repl = rep(1, times = mesh_barrier$n*no_T),
                          n.repl = 3)

A_proj_K = inla.spde.make.A(mesh, coordinates(total_sightings_K),
                            group = total_sightings_K$MONTH_INLA,
                            n.group = no_T,
                            repl = total_sightings_K$POD_INLA,
                            n.repl = 3)
A_proj_barrier_K = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_K),
                                    group = total_sightings_K$MONTH_INLA,
                                    n.group = no_T,
                                    repl = total_sightings_K$POD_INLA,
                                    n.repl = 3)
imat_K = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                          group = rep(1:no_T,each = mesh$n),
                          n.group = no_T,
                          repl = rep(2, times = mesh$n*no_T),
                          n.repl = 3)
imat_barrier_K = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                                  group = rep(1:no_T,each = mesh_barrier$n),
                                  n.group = no_T,
                                  repl = rep(2, times = mesh_barrier$n*no_T),
                                  n.repl = 3)

A_proj_L = inla.spde.make.A(mesh, coordinates(total_sightings_L),
                            group = total_sightings_L$MONTH_INLA,
                            n.group = no_T,
                            repl = total_sightings_L$POD_INLA,
                            n.repl = 3)
A_proj_barrier_L = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_L),
                                    group = total_sightings_L$MONTH_INLA,
                                    n.group = no_T,
                                    repl = total_sightings_L$POD_INLA,
                                    n.repl = 3)
imat_L = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                          group = rep(1:no_T,each = mesh$n),
                          n.group = no_T,
                          repl = rep(3, times = mesh$n*no_T),
                          n.repl = 3)
imat_barrier_L = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                                  group = rep(1:no_T,each = mesh_barrier$n),
                                  n.group = no_T,
                                  repl = rep(3, times = mesh_barrier$n*no_T),
                                  n.repl = 3)
# First mesh$n are zeros and the remaining dim(BlackSmokePrefData2)[1] equal R
w2 = rep(w, times=no_T)
w2_barrier = rep(w_barrier, times=no_T)

y.pp_J <- c(rep(0, mesh$n*no_T), rep(1, dim(total_sightings_J)[1]))
y.pp_J[1:(mesh$n*no_T)][w2==0] = NA # set the points on land to NA as they are not integration points

y.pp_barrier_J <- c(rep(0, mesh_barrier$n*no_T), rep(1, dim(total_sightings_J)[1]))
y.pp_barrier_J[1:(mesh_barrier$n*no_T)][w2_barrier==0] = NA

y.pp_K <- c(rep(0, mesh$n*no_T), rep(1, dim(total_sightings_K)[1]))
y.pp_K[1:(mesh$n*no_T)][w2==0] = NA # set the points on land to NA as they are not integration points

y.pp_barrier_K <- c(rep(0, mesh_barrier$n*no_T), rep(1, dim(total_sightings_K)[1]))
y.pp_barrier_K[1:(mesh_barrier$n*no_T)][w2_barrier==0] = NA

y.pp_L <- c(rep(0, mesh$n*no_T), rep(1, dim(total_sightings_L)[1]))
y.pp_L[1:(mesh$n*no_T)][w2==0] = NA # set the points on land to NA as they are not integration points

y.pp_barrier_L <- c(rep(0, mesh_barrier$n*no_T), rep(1, dim(total_sightings_L)[1]))
y.pp_barrier_L[1:(mesh_barrier$n*no_T)][w2_barrier==0] = NA
# spacetime below
#y.pp <- c(rep(0, mesh$n*no_T),  rep(1, dim(total_sightings)[1]))

## ----expected------------------------------------------------------------
e.pp_J <- c(w2, rep(0, dim(total_sightings_J)[1])) 
e.pp_barrier_J <- c(w2_barrier, rep(0, dim(total_sightings_J)[1])) 

e.pp_K <- c(w2, rep(0, dim(total_sightings_K)[1])) 
e.pp_barrier_K <- c(w2_barrier, rep(0, dim(total_sightings_K)[1])) 

e.pp_L <- c(w2, rep(0, dim(total_sightings_L)[1])) 
e.pp_barrier_L <- c(w2_barrier, rep(0, dim(total_sightings_L)[1])) 

# scale w so the intercept is less affected by the prior (closer to 0)
#e.pp <- c(w2/max(w2), rep(0, dim(total_sightings)[1])) 
simple.model = inla.spde2.pcmatern(mesh_transformed, 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1)) 

s_index = inla.spde.make.index(name = "simple.field",
                               n.spde = simple.model$n.spde,
                               n.group = no_T, n.repl = 3)
s_index_barrier = inla.spde.make.index(name = "barrier.field",
                                       n.spde = barrier.model$f$n,
                                       n.group = no_T, n.repl = 3)

## ----App-----------------------------------------------------------------
# join the Identity matrix for the mesh nodes with the projection matrix for the obs
A.pp_constr_J <- rBind(imat_J, A_proj_J)
A.pp_constr_barrier_J <- rBind(imat_barrier_J, A_proj_barrier_J)

A.pp_constr_K <- rBind(imat_K, A_proj_K)
A.pp_constr_barrier_K <- rBind(imat_barrier_K, A_proj_barrier_K)

A.pp_constr_L <- rBind(imat_L, A_proj_L)
A.pp_constr_barrier_L <- rBind(imat_barrier_L, A_proj_barrier_L)

### Modifying the Search Effort Layers
## USE THE SCRIPT 2009_2016_WW_SI_modification.R to change WW search effort!
# load the updated files now and run the below lines of code 
# load()

## Link the covariates to the observations
# loop over the pods
covariates_site_J = total_sightings_J@data

# remove unused columns to enable binding to work
covariates_site_J = covariates_site_J[,c("YEAR","MONTH","J","K","L","POD_INLA","chloro", "SST", "vessel","topo")]
#covariates_pp = covariates_pp[,-c(5)]
# Keep MONTH, chloro, SST, vessel, MONTH_INLA and topo

#############################
# Smooth the WW search effort layer
plot.field(log(SI_complete_monthly_J_mean_barrier[,1]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
           pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean_barrier)+1)) + ggtitle('May WW log Search Effort J pod')
# very pixellated. Use smoother SI_WW from coarser mesh

SI_WW_complete_monthly_J_mean_old =  apply(SI_WW_complete_monthly_J_mean, 2 ,FUN=function(x){x/(w/mean(w))})
SI_WW_complete_monthly_K_mean_old =  apply(SI_WW_complete_monthly_K_mean, 2 ,FUN=function(x){x/(w/mean(w))})
SI_WW_complete_monthly_L_mean_old =  apply(SI_WW_complete_monthly_L_mean, 2 ,FUN=function(x){x/(w/mean(w))})

SI_WW_complete_monthly_J_mean_old[w==0,] = 0
SI_WW_complete_monthly_K_mean_old[w==0,] = 0
SI_WW_complete_monthly_L_mean_old[w==0,] = 0

proj_barrier_nonbarrier = inla.spde.make.A(mesh_transformed,mesh_barrier_transformed$loc[,c(1:2)])
SI_WW_complete_monthly_J_mean_barrier_new = 
  as.matrix(proj_barrier_nonbarrier %*% SI_WW_complete_monthly_J_mean_old)
plot.field(log(SI_WW_complete_monthly_J_mean_barrier_new[,1]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
           pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_WW_complete_monthly_J_mean_barrier_new)+1)) + ggtitle('May WW log Search Effort J pod')
# fixed
SI_WW_complete_monthly_J_mean_barrier = SI_WW_complete_monthly_J_mean_barrier_new
SI_WW_complete_monthly_K_mean_barrier = as.matrix(proj_barrier_nonbarrier %*% SI_WW_complete_monthly_K_mean_old)
SI_WW_complete_monthly_L_mean_barrier = as.matrix(proj_barrier_nonbarrier %*% SI_WW_complete_monthly_L_mean_old)

# scale by area
SI_WW_complete_monthly_J_mean_barrier =  apply(SI_WW_complete_monthly_J_mean_barrier, 2 ,FUN=function(x){x*(w_barrier/mean(w_barrier))})
SI_WW_complete_monthly_K_mean_barrier =  apply(SI_WW_complete_monthly_K_mean_barrier, 2 ,FUN=function(x){x*(w_barrier/mean(w_barrier))})
SI_WW_complete_monthly_L_mean_barrier =  apply(SI_WW_complete_monthly_L_mean_barrier, 2 ,FUN=function(x){x*(w_barrier/mean(w_barrier))})

SI_WW_complete_monthly_J_mean_barrier[w_barrier==0,] = 0
SI_WW_complete_monthly_K_mean_barrier[w_barrier==0,] = 0
SI_WW_complete_monthly_L_mean_barrier[w_barrier==0,] = 0

# rescale
SI_WW_complete_monthly_J_mean_barrier = rescale_fun(SI_WW_complete_monthly_J_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_J))
SI_WW_complete_monthly_K_mean_barrier = rescale_fun(SI_WW_complete_monthly_K_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_K))
SI_WW_complete_monthly_L_mean_barrier = rescale_fun(SI_WW_complete_monthly_L_mean_barrier, colSums(mean_total_boat_hours_per_period_per_port_L))

SI_complete_monthly_J_mean_barrier = SI_WW_complete_monthly_J_mean_barrier + SI_Brian_barrier_monthly_J
SI_complete_monthly_K_mean_barrier = SI_WW_complete_monthly_K_mean_barrier + SI_Brian_barrier_monthly_K
SI_complete_monthly_L_mean_barrier = SI_WW_complete_monthly_L_mean_barrier + SI_Brian_barrier_monthly_L

plot.field(log(SI_complete_monthly_J_mean_barrier[,1]+1), poly = COAST_transformed, mesh = mesh_barrier_transformed, 
           pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean_barrier)+1)) + ggtitle('May WW log Search Effort J pod')
##############################

# Remember chloro, vessel and depth all should be log transformed
# covariates_pp and covariates_pp_barrier have already been formed

covariates_site_J$chloro = log(covariates_site_J$chloro)
covariates_site_J$vessel = log(covariates_site_J$vessel+1)

# find the pixel corresponding to each observation
which_dmesh_pixel_J = gWithin(total_sightings_J,dmesh, byid = T, returnDense = F)
which_dmesh_pixel_J = unlist(which_dmesh_pixel_J)
which_dmesh_pixel_J = as.numeric(which_dmesh_pixel_J)

# map the spatiotemporal averages onto the points
covariates_site_J$SSTmonthavg = 0
covariates_site_J$chloromonthavg = 0

count = 1
for(i in unique(covariates_pp$MONTH))
{
  covariates_site_J$SSTmonthavg[covariates_site_J$MONTH == i] = SST_monthlymeans[count]
  covariates_site_J$chloromonthavg[covariates_site_J$MONTH == i] = chloro_monthlymeans[count]
  count = count + 1
}

covariates_site_J$SSTminusmonth = covariates_site_J$SST - covariates_site_J$SSTmonthavg # how warm is the location compared to monthly average
covariates_site_J$chlorominusmonth = covariates_site_J$chloro - covariates_site_J$chloromonthavg
covariates_site_J$SSTmonthavg = covariates_site_J$SSTmonthavg - SST_overallmean
covariates_site_J$chloromonthavg = covariates_site_J$chloromonthavg - chloro_overallmean

covariates_site_J$SSTspatialavg = SST_spatialmeans[which_dmesh_pixel_J] - SST_overallmean # is the region typically? warm 
covariates_site_J$chlorospatialavg = chloro_spatialmeans[which_dmesh_pixel_J] - chloro_overallmean

covariates_site_J$SSTminusspacetime = covariates_site_J$SST - SST_overallmean - covariates_site_J$SSTmonthavg - covariates_site_J$SSTspatialavg
covariates_site_J$chlorominusspacetime = covariates_site_J$chloro - chloro_overallmean - covariates_site_J$chloromonthavg - covariates_site_J$chlorospatialavg

# create a depth variable and take the log to reduce the skew.
# Note we take the absolute value without worry - the topo's above 51 lie on land and thus are not in AOI
covariates_site_J$depth = log(abs(covariates_site_J$topo-51)+1)
covariates_site_J$depth = covariates_site_J$depth - mean_depth_dmesh
covariates_site_J$topo = covariates_site_J$topo - mean(covariates_pp$topo)

# add the sampling intensity
covariates_pp_J = covariates_pp
covariates_pp_J$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_J_mean))
covariates_pp_J$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_J_var)))

covariates_pp_barrier_J = covariates_pp_barrier
covariates_pp_barrier_J$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_J_mean_barrier))
covariates_pp_barrier_J$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_J_var_barrier)))

covariates_site_J$Search_Effort = 0
covariates_site_J$Search_Effort_SD = 0
covariates_site_J$MONTH_INLA = covariates_site_J$MONTH - 4
# Join the covariates of the locations PP and the nodes PP together(nodes first)
#covariates_pp = covariates_pp[,-c(7)]
temp_names = intersect(names(covariates_pp_J), names(covariates_site_J))
covariates_PP_join_J = rbind.data.frame(covariates_pp_J[,temp_names], 
                                        covariates_site_J[,temp_names])
covariates_PP_join_J$POD_INLA = 1

covariates_PP_join_barrier_J = rbind.data.frame(covariates_pp_barrier_J[,temp_names], 
                                                covariates_site_J[,temp_names])
covariates_PP_join_barrier_J$POD_INLA = 1
#Pod K

covariates_site_K = total_sightings_K@data

# remove unused columns to enable binding to work
covariates_site_K = covariates_site_K[,c("YEAR","MONTH","J","K","L","POD_INLA","chloro", "SST", "vessel","topo")]
#covariates_pp = covariates_pp[,-c(5)]
# Keep MONTH, chloro, SST, vessel, MONTH_INLA and topo

# Remember chloro, vessel and depth all should be log transformed
# covariates_pp and covariates_pp_barrier have already been formed

covariates_site_K$chloro = log(covariates_site_K$chloro)
covariates_site_K$vessel = log(covariates_site_K$vessel+1)

# find the pixel corresponding to each observation
which_dmesh_pixel_K = gWithin(total_sightings_K,dmesh, byid = T, returnDense = F)
which_dmesh_pixel_K = unlist(which_dmesh_pixel_K)
which_dmesh_pixel_K = as.numeric(which_dmesh_pixel_K)

# map the spatiotemporal averages onto the points
covariates_site_K$SSTmonthavg = 0
covariates_site_K$chloromonthavg = 0

count = 1
for(i in unique(covariates_pp$MONTH))
{
  covariates_site_K$SSTmonthavg[covariates_site_K$MONTH == i] = SST_monthlymeans[count]
  covariates_site_K$chloromonthavg[covariates_site_K$MONTH == i] = chloro_monthlymeans[count]
  count = count + 1
}

covariates_site_K$SSTminusmonth = covariates_site_K$SST - covariates_site_K$SSTmonthavg # how warm is the location compared to monthly average
covariates_site_K$chlorominusmonth = covariates_site_K$chloro - covariates_site_K$chloromonthavg
covariates_site_K$SSTmonthavg = covariates_site_K$SSTmonthavg - SST_overallmean
covariates_site_K$chloromonthavg = covariates_site_K$chloromonthavg - chloro_overallmean

covariates_site_K$SSTspatialavg = SST_spatialmeans[which_dmesh_pixel_K] - SST_overallmean # is the region typically? warm 
covariates_site_K$chlorospatialavg = chloro_spatialmeans[which_dmesh_pixel_K] - chloro_overallmean

covariates_site_K$SSTminusspacetime = covariates_site_K$SST - SST_overallmean - covariates_site_K$SSTmonthavg - covariates_site_K$SSTspatialavg
covariates_site_K$chlorominusspacetime = covariates_site_K$chloro - chloro_overallmean - covariates_site_K$chloromonthavg - covariates_site_K$chlorospatialavg

# create a depth variable and take the log to reduce the skew.
# Note we take the absolute value without worry - the topo's above 51 lie on land and thus are not in AOI
covariates_site_K$depth = log(abs(covariates_site_K$topo-51)+1)
covariates_site_K$depth = covariates_site_K$depth - mean_depth_dmesh
covariates_site_K$topo = covariates_site_K$topo - mean(covariates_pp$topo)

# add the sampling intensity
covariates_pp_K = covariates_pp
covariates_pp_K$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_K_mean))
covariates_pp_K$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_K_var)))

covariates_pp_barrier_K = covariates_pp_barrier
covariates_pp_barrier_K$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_K_mean_barrier))
covariates_pp_barrier_K$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_K_var_barrier)))

covariates_site_K$Search_Effort = 0
covariates_site_K$Search_Effort_SD = 0
covariates_site_K$MONTH_INLA = covariates_site_K$MONTH - 4
# Join the covariates of the locations PP and the nodes PP together(nodes first)
#covariates_pp = covariates_pp[,-c(7)]
covariates_PP_join_K = rbind.data.frame(covariates_pp_K[,temp_names], 
                                        covariates_site_K[,temp_names])
covariates_PP_join_K$POD_INLA = 2

covariates_PP_join_barrier_K = rbind.data.frame(covariates_pp_barrier_K[,temp_names], 
                                                covariates_site_K[,temp_names])
covariates_PP_join_barrier_K$POD_INLA = 2
#Pod L

covariates_site_L = total_sightings_L@data

# remove unused columns to enable binding to work
covariates_site_L = covariates_site_L[,c("YEAR","MONTH","J","K","L","POD_INLA","chloro", "SST", "vessel","topo")]
#covariates_pp = covariates_pp[,-c(5)]
# Keep MONTH, chloro, SST, vessel, MONTH_INLA and topo

# Remember chloro, vessel and depth all should be log transformed
# covariates_pp and covariates_pp_barrier have already been formed

covariates_site_L$chloro = log(covariates_site_L$chloro)
covariates_site_L$vessel = log(covariates_site_L$vessel+1)

# find the pixel corresponding to each observation
which_dmesh_pixel_L = gWithin(total_sightings_L,dmesh, byid = T, returnDense = F)
which_dmesh_pixel_L = unlist(which_dmesh_pixel_L)
which_dmesh_pixel_L = as.numeric(which_dmesh_pixel_L)

# map the spatiotemporal averages onto the points
covariates_site_L$SSTmonthavg = 0
covariates_site_L$chloromonthavg = 0

count = 1
for(i in unique(covariates_pp$MONTH))
{
  covariates_site_L$SSTmonthavg[covariates_site_L$MONTH == i] = SST_monthlymeans[count]
  covariates_site_L$chloromonthavg[covariates_site_L$MONTH == i] = chloro_monthlymeans[count]
  count = count + 1
}

covariates_site_L$SSTminusmonth = covariates_site_L$SST - covariates_site_L$SSTmonthavg # how warm is the location compared to monthly average
covariates_site_L$chlorominusmonth = covariates_site_L$chloro - covariates_site_L$chloromonthavg
covariates_site_L$SSTmonthavg = covariates_site_L$SSTmonthavg - SST_overallmean
covariates_site_L$chloromonthavg = covariates_site_L$chloromonthavg - chloro_overallmean

covariates_site_L$SSTspatialavg = SST_spatialmeans[which_dmesh_pixel_L] - SST_overallmean # is the region typically? warm 
covariates_site_L$chlorospatialavg = chloro_spatialmeans[which_dmesh_pixel_L] - chloro_overallmean

covariates_site_L$SSTminusspacetime = covariates_site_L$SST - SST_overallmean - covariates_site_L$SSTmonthavg - covariates_site_L$SSTspatialavg
covariates_site_L$chlorominusspacetime = covariates_site_L$chloro - chloro_overallmean - covariates_site_L$chloromonthavg - covariates_site_L$chlorospatialavg

# create a depth variable and take the log to reduce the skew.
# Note we take the absolute value without worry - the topo's above 51 lie on land and thus are not in AOI
covariates_site_L$depth = log(abs(covariates_site_L$topo-51)+1)
covariates_site_L$depth = covariates_site_L$depth - mean_depth_dmesh
covariates_site_L$topo = covariates_site_L$topo - mean(covariates_pp$topo)

# add the sampling intensity
covariates_pp_L = covariates_pp
covariates_pp_L$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_L_mean))
covariates_pp_L$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_L_var)))

covariates_pp_barrier_L = covariates_pp_barrier
covariates_pp_barrier_L$Search_Effort = as.numeric(as.matrix(SI_complete_monthly_L_mean_barrier))
covariates_pp_barrier_L$Search_Effort_SD = sqrt(as.numeric(as.matrix(SI_complete_monthly_L_var_barrier)))

covariates_site_L$Search_Effort = 0
covariates_site_L$Search_Effort_SD = 0
covariates_site_L$MONTH_INLA = covariates_site_L$MONTH - 4
# Join the covariates of the locations PP and the nodes PP together(nodes first)
#covariates_pp = covariates_pp[,-c(7)]
covariates_PP_join_L = rbind.data.frame(covariates_pp_L[,temp_names], 
                                        covariates_site_L[,temp_names])
covariates_PP_join_L$POD_INLA = 3

covariates_PP_join_barrier_L = rbind.data.frame(covariates_pp_barrier_L[,temp_names], 
                                                covariates_site_L[,temp_names])
covariates_PP_join_barrier_L$POD_INLA = 3

barrier.model = inla.barrier.pcmatern(mesh_barrier_transformed, barrier.triangles = polygon.triangles, 
                                      prior.range = c(15, 0.01), 
                                      prior.sigma = c(3, 0.1), 
                                      range.fraction = 0.2) # remember we're not in lon/lat, 3e5 is domain size
# reasonable lower range (median) is half of domain prior.range = c(15e4, .5)
# reasonable upper scale is sd of 4

simple.model = inla.spde2.pcmatern(mesh_transformed, 
                                   prior.range = c(15, 0.01), 
                                   prior.sigma = c(3, 0.1)) 


# create stack matrix for INLA manually
# create the data stack for the PP
max_effort = max(c(covariates_PP_join_J$Search_Effort,
                   covariates_PP_join_K$Search_Effort,
                   covariates_PP_join_L$Search_Effort),
                 na.rm=T)
max_effort_barrier = max(c(covariates_PP_join_barrier_J$Search_Effort,
                   covariates_PP_join_barrier_K$Search_Effort,
                   covariates_PP_join_barrier_L$Search_Effort),
                 na.rm=T)

A.pp_J = A.pp_constr_J
A.pp_K = A.pp_constr_K
A.pp_L = A.pp_constr_L

A.pp_barrier_J = A.pp_constr_barrier_J
A.pp_barrier_K = A.pp_constr_barrier_K
A.pp_barrier_L = A.pp_constr_barrier_L

covariates_PP_join_J$Intercept = 1
covariates_PP_join_J$Search_Effort_SD_scaled = covariates_PP_join_J$Search_Effort_SD / max_effort
PP_stack_J = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                        A=list(1,A.pp_J), tag='pp',
                        effects=list(covariates_PP_join_J, 
                                     s_index)) # i are the indices for the SPDE 
print('stack pp J complete')

covariates_PP_join_barrier_J$Intercept = 1
covariates_PP_join_barrier_J$Search_Effort_SD_scaled = covariates_PP_join_barrier_J$Search_Effort_SD / max_effort_barrier
PP_stack_barrier_J = inla.stack(data=list(y=y.pp_barrier_J, e=covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier),#(e.pp_barrier_J/max(e.pp_barrier_J))*(covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier) ),
                                A=list(1,A.pp_barrier_J), tag='pp',
                                effects=list(covariates_PP_join_barrier_J, 
                                             s_index_barrier)) # i are the indices for the SPDE 
print('stack pp barrier J complete')

covariates_PP_join_K$Intercept = 1
covariates_PP_join_K$Search_Effort_SD_scaled = covariates_PP_join_K$Search_Effort_SD / max_effort
PP_stack_K = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                        A=list(1,A.pp_K), tag='pp',
                        effects=list(covariates_PP_join_K, 
                                     s_index)) # i are the indices for the SPDE 
print('stack pp K complete')

covariates_PP_join_barrier_K$Intercept = 1
covariates_PP_join_barrier_K$Search_Effort_SD_scaled = covariates_PP_join_barrier_K$Search_Effort_SD / max_effort_barrier
PP_stack_barrier_K = inla.stack(data=list(y=y.pp_barrier_K, e=covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                A=list(1,A.pp_barrier_K), tag='pp',
                                effects=list(covariates_PP_join_barrier_K, 
                                             s_index_barrier)) # i are the indices for the SPDE 
print('stack pp barrier K complete')

covariates_PP_join_L$Intercept = 1
covariates_PP_join_L$Search_Effort_SD_scaled = covariates_PP_join_L$Search_Effort_SD / max_effort
PP_stack_L = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
                        A=list(1,A.pp_L), tag='pp',
                        effects=list(covariates_PP_join_L, 
                                     s_index)) # i are the indices for the SPDE 
print('stack pp L complete')

covariates_PP_join_barrier_L$Intercept = 1
covariates_PP_join_barrier_L$Search_Effort_SD_scaled = covariates_PP_join_barrier_L$Search_Effort_SD / max_effort_barrier
PP_stack_barrier_L = inla.stack(data=list(y=y.pp_barrier_L, e=covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier),#(e.pp_barrier_L/max(e.pp_barrier_L))*(covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier) ),
                                A=list(1,A.pp_barrier_L), tag='pp',
                                effects=list(covariates_PP_join_barrier_L, 
                                             s_index_barrier)) # i are the indices for the SPDE 
print('stack pp barrier L complete')

# revisede.RData version has removed the scaling by e.pp which was erroneous
#save.image('SI_WS4_revisedSI_revisede.RData')
#save.image('SI_WS4_revisedSI_dup_removed.RData')
load('SI_WS4_revisedSI_dup_removed.RData')

model_fitting_files = list(PP_stack_J = PP_stack_J,
                           PP_stack_K = PP_stack_K,
                           PP_stack_L = PP_stack_L,
                           PP_stack_barrier_J = PP_stack_barrier_J,
                           PP_stack_barrier_K = PP_stack_barrier_K,
                           PP_stack_barrier_L = PP_stack_barrier_L,
                           barrier.model = barrier.model,
                           simple.model = simple.model)


#saveRDS(model_fitting_files, 'Model_fitting_files_scaled.rds')
#saveRDS(model_fitting_files, 'Model_fitting_files_scaled_new.rds')
#saveRDS(model_fitting_files, 'Model_fitting_files_scaled_dup_removed.rds')


# Save the files needed for simulating and plotting
#saveRDS(covariates_pp_J,'covariates_pp.rds') # same for inla mesh
#saveRDS(covariates_pp_J,'covariates_pp_dup_removed.rds')
#saveRDS(covariates_pp_barrier_J,'covariates_pp_barrier_dup_removed.rds')


#### Compare the search effort from Brian vs WW companies!
colSums(SI_WW_complete_monthly_J_mean)
colSums(SI_Brian_monthly_J)

# average monthly WW search effort
mean(colSums(SI_WW_complete_monthly_L_mean)) # ~ 20,000 boat hours

#### Find the mesh vertices that lie in the water - needed for restricted spatial regression constraints
polygon.triangles_mesh = inla.over_sp_mesh(spTransform(COAST_simp, CRSobj = mesh$crs ), y = mesh, type = "vertex", ignore.CRS=T)
plot(mesh$loc[polygon.triangles_mesh,])

#### Add an spde for depth
## Need to recreate depth variable due to error
# covariates_PP_join_J$depth = covariates_PP_join_J$topo
# #covariates_PP_join_J$depth[covariates_PP_join_J$depth>=51] = 51
# covariates_PP_join_J$depth = log(abs(covariates_PP_join_J$depth - 52)+1)
# hist(covariates_PP_join_J$depth[covariates_PP_join_J$depth >0] )
# covariates_PP_join_J$depth = scale(covariates_PP_join_J$depth)
# hist(covariates_PP_join_J$depth)

# covariates_PP_join_K$depth = covariates_PP_join_K$topo
# #covariates_PP_join_K$depth[covariates_PP_join_K$depth>=51] = 51
# covariates_PP_join_K$depth = log(abs(covariates_PP_join_K$depth - 52))
# hist(covariates_PP_join_K$depth[covariates_PP_join_K$depth >0] )
# covariates_PP_join_K$depth = scale(covariates_PP_join_K$depth)
# 
# covariates_PP_join_L$depth = covariates_PP_join_L$topo
# #covariates_PP_join_L$depth[covariates_PP_join_L$depth>=51] = 51
# covariates_PP_join_L$depth = log(abs(covariates_PP_join_L$depth - 52))
# hist(covariates_PP_join_L$depth[covariates_PP_join_L$depth >0] )
# covariates_PP_join_L$depth = scale(covariates_PP_join_L$depth)

depth_mesh = inla.mesh.1d(loc = seq(from = min(covariates_PP_join_J$depth,na.rm=T)-1,
                                    to = max(covariates_PP_join_J$depth,na.rm=T)+1, length.out = 20),
                          degree = 2)
depth_spde = inla.spde2.pcmatern(depth_mesh, constr = T, prior.range = c(5, 0.95), prior.sigma = c(1,0.1))
depth_index = inla.spde.make.index('depth_smooth', n.spde = depth_spde$n.spde)

A.pp_J_depth = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_J$depth)
PP_stack_J_depth = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                              A=list(1, A.pp_J, A.pp_J_depth), tag='pp',
                              effects=list(covariates_PP_join_J, 
                                           s_index, depth_index)) # i are the indices for the SPDE 
print('stack pp J complete')

A.pp_K_depth = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_K$depth)
PP_stack_K_depth = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                              A=list(1,A.pp_K, A.pp_K_depth), tag='pp',
                              effects=list(covariates_PP_join_K, 
                                           s_index, depth_index)) # i are the indices for the SPDE 
print('stack pp K complete')

A.pp_L_depth = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_L$depth)
PP_stack_L_depth = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
                              A=list(1,A.pp_L, A.pp_L_depth), tag='pp',
                              effects=list(covariates_PP_join_L, 
                                           s_index, depth_index)) # i are the indices for the SPDE 
print('stack pp L complete')

model_fitting_files_depth = list(PP_stack_J_depth = PP_stack_J_depth,
                                 PP_stack_K_depth = PP_stack_K_depth,
                                 PP_stack_L_depth = PP_stack_L_depth,
                                 depth_spde = depth_spde)


#saveRDS(model_fitting_files_depth, 'Model_fitting_files_scaled_depthspde.rds')
#saveRDS(model_fitting_files_depth, 'Model_fitting_files_scaled_depthspde_dup_removed.rds')

## Add a second spde for pods K and L to allow for pod-specific spatial effects

# Add a second spde index and projector matrix
s_index_J = inla.spde.make.index(name = "simple.field.J",
                                 n.spde = simple.model$n.spde,
                                 n.group = no_T, n.repl = 3)
s_index_J$simple.field.J[s_index_J$simple.field.J %in% which(w == 0)] = NA
s_index_barrier_J = inla.spde.make.index(name = "barrier.field.J",
                                         n.spde = barrier.model$f$n,
                                         n.group = no_T, n.repl = 3)
s_index_barrier_J$barrier.field.J[s_index_barrier_J$barrier.field.J %in% which(w_barrier == 0)] = NA

s_index_K = inla.spde.make.index(name = "simple.field.K",
                                 n.spde = simple.model$n.spde,
                                 n.group = no_T, n.repl = 3)
s_index_K$simple.field.K[s_index_K$simple.field.K %in% which(w == 0)] = NA
s_index_barrier_K = inla.spde.make.index(name = "barrier.field.K",
                                         n.spde = barrier.model$f$n,
                                         n.group = no_T, n.repl = 3)
s_index_barrier_K$barrier.field.K[s_index_barrier_K$barrier.field.K %in% which(w_barrier == 0)] = NA

s_index_L = inla.spde.make.index(name = "simple.field.L",
                                 n.spde = simple.model$n.spde,
                                 n.group = no_T, n.repl = 3)
s_index_barrier_L = inla.spde.make.index(name = "barrier.field.L",
                                         n.spde = barrier.model$f$n,
                                         n.group = no_T, n.repl = 3)
s_index_L2 = inla.spde.make.index(name = "simple.field.L2",
                                 n.spde = simple.model$n.spde,
                                 n.group = no_T, n.repl = 3)
s_index_L2$simple.field.L2[s_index_L2$simple.field.L2 %in% which(w == 0)] = NA
s_index_barrier_L2 = inla.spde.make.index(name = "barrier.field.L2",
                                         n.spde = barrier.model$f$n,
                                         n.group = no_T, n.repl = 3)
s_index_barrier_L2$barrier.field.L2[s_index_barrier_L2$barrier.field.L2 %in% which(w_barrier == 0)] = NA



PP_stack_J_podSPDE = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                                A=list(1,A.pp_J, A.pp_J), tag='pp',
                                effects=list(covariates_PP_join_J, 
                                             s_index, s_index_J)) # i are the indices for the SPDE 
print('stack pp J complete')

PP_stack_barrier_J_podSPDE = inla.stack(data=list(y=y.pp_barrier_J, e=covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier),#(e.pp_barrier_J/max(e.pp_barrier_J))*(covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier) ),
                                        A=list(1,A.pp_barrier_J, A.pp_barrier_J), tag='pp',
                                        effects=list(covariates_PP_join_barrier_J, 
                                                     s_index_barrier, s_index_barrier_J)) # i are the indices for the SPDE 
print('stack pp barrier J complete')

PP_stack_K_podSPDE = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                                A=list(1,A.pp_K, A.pp_K), tag='pp',
                                effects=list(covariates_PP_join_K, 
                                             s_index, s_index_K)) # i are the indices for the SPDE 
print('stack pp K complete')

PP_stack_barrier_K_podSPDE = inla.stack(data=list(y=y.pp_barrier_K, e=covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                        A=list(1,A.pp_barrier_K, A.pp_barrier_K), tag='pp',
                                        effects=list(covariates_PP_join_barrier_K, 
                                                     s_index_barrier, s_index_barrier_K)) # i are the indices for the SPDE 
print('stack pp barrier K complete')

PP_stack_L_podSPDE = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                                A=list(1,A.pp_L, A.pp_L, A.pp_L), tag='pp',
                                effects=list(covariates_PP_join_L, 
                                             s_index, s_index_L, s_index_L2)) # i are the indices for the SPDE 
print('stack pp L complete')

PP_stack_barrier_L_podSPDE = inla.stack(data=list(y=y.pp_barrier_L, e=covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                        A=list(1,A.pp_barrier_L, A.pp_barrier_L, A.pp_barrier_L), tag='pp',
                                        effects=list(covariates_PP_join_barrier_L, 
                                                     s_index_barrier, s_index_barrier_L, s_index_barrier_L2)) # i are the indices for the SPDE 
print('stack pp barrier L complete')

model_fitting_files_podSPDE = list(PP_stack_J_podSPDE = PP_stack_J_podSPDE,
                                   PP_stack_K_podSPDE = PP_stack_K_podSPDE,
                                   PP_stack_L_podSPDE = PP_stack_L_podSPDE,
                                   PP_stack_barrier_J_podSPDE = PP_stack_barrier_J_podSPDE,
                                   PP_stack_barrier_K_podSPDE = PP_stack_barrier_K_podSPDE,
                                   PP_stack_barrier_L_podSPDE = PP_stack_barrier_L_podSPDE)

#saveRDS(model_fitting_files_podSPDE, 'model_fitting_files_podSPDE.rds')
saveRDS(model_fitting_files_podSPDE, 'model_fitting_files_podSPDE_dup_removed2.rds')

# Add a new pod index and a new depth spde
PP_stack_J_depth_podSPDE = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                              A=list(1, A.pp_J, A.pp_J_depth), tag='pp',
                              effects=list(covariates_PP_join_J, 
                                           s_index, depth_index)) # i are the indices for the SPDE 
print('stack pp J complete')

A.pp_J_depth_barrier = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_barrier_J$depth)
PP_stack_barrier_J_depth_podSPDE = inla.stack(data=list(y=y.pp_barrier_J, e=covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                        A=list(1,A.pp_barrier_J, A.pp_J_depth_barrier), tag='pp',
                                        effects=list(covariates_PP_join_barrier_J, 
                                                     s_index_barrier, depth_index)) # i are the indices for the SPDE 
print('stack pp barrier J complete')

A.pp_K_depth = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_K$depth)
PP_stack_K_depth_podSPDE = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                              A=list(1,A.pp_K, A.pp_K, A.pp_K_depth), tag='pp',
                              effects=list(covariates_PP_join_K, 
                                           s_index, s_index_K, depth_index)) # i are the indices for the SPDE 
print('stack pp K complete')
A.pp_K_depth_barrier = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_barrier_K$depth)
PP_stack_barrier_K_depth_podSPDE = inla.stack(data=list(y=y.pp_barrier_K, e=covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                              A=list(1,A.pp_barrier_K, A.pp_barrier_K, A.pp_K_depth_barrier), tag='pp',
                                              effects=list(covariates_PP_join_barrier_K, 
                                                           s_index_barrier, s_index_barrier_K, depth_index)) # i are the indices for the SPDE 
print('stack pp barrier K complete')

A.pp_L_depth = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_L$depth)
PP_stack_L_depth_podSPDE = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
                              A=list(1,A.pp_L, A.pp_L, A.pp_L_depth), tag='pp',
                              effects=list(covariates_PP_join_L, 
                                           s_index, s_index_L, depth_index)) # i are the indices for the SPDE 
print('stack pp L complete')
A.pp_L_depth_barrier = inla.spde.make.A(mesh = depth_mesh, loc = covariates_PP_join_barrier_L$depth)
PP_stack_barrier_L_depth_podSPDE = inla.stack(data=list(y=y.pp_barrier_L, e=covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier),#(e.pp_barrier_L/max(e.pp_barrier_L))*(covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier) ),
                                              A=list(1,A.pp_barrier_L,A.pp_barrier_L, A.pp_L_depth_barrier), tag='pp',
                                              effects=list(covariates_PP_join_barrier_L, 
                                                           s_index_barrier, s_index_barrier_L, depth_index)) # i are the indices for the SPDE 
print('stack pp barrier L complete')

model_fitting_files_depth_podSPDE = list(PP_stack_J_depth_podSPDE = PP_stack_J_depth_podSPDE,
                                 PP_stack_K_depth_podSPDE = PP_stack_K_depth_podSPDE,
                                 PP_stack_L_depth_podSPDE = PP_stack_L_depth_podSPDE,
                                 PP_stack_barrier_J_depth_podSPDE = PP_stack_barrier_J_depth_podSPDE,
                                 PP_stack_barrier_K_depth_podSPDE = PP_stack_barrier_K_depth_podSPDE,
                                 PP_stack_barrier_L_depth_podSPDE = PP_stack_barrier_L_depth_podSPDE,
                                 depth_spde = depth_spde)


#saveRDS(model_fitting_files_depth, 'Model_fitting_files_scaled_depthspde.rds')
saveRDS(model_fitting_files_depth_podSPDE, 'Model_fitting_files_scaled_depthspde_podspde_dup_removed.rds')


# Save the files needed for simulating and plotting
#saveRDS(covariates_pp_J,'covariates_pp.rds') # same for inla mesh
saveRDS(mesh,'mesh.rds') 
saveRDS(mesh_barrier,'mesh_barrier.rds')
saveRDS(mesh_transformed,'mesh_transformed.rds') 
saveRDS(mesh_barrier_transformed,'mesh_barrier_transformed.rds') 

proj = inla.spde.make.A(mesh_transformed, loc = coordinates(pixels_plotting))
proj_barrier = inla.spde.make.A(mesh_barrier_transformed, loc = coordinates(pixels_plotting))

# New files for posterior plots
pixels_mesh = inlabru::pixels(mesh_transformed ,mask = spTransform(COAST_simp, CRSobj = COAST_transformed@proj4string), nx = 300, ny = 300 )
proj_meshpoly = inla.spde.make.A(mesh_transformed, loc = coordinates(pixels_mesh))
proj_barrier_meshpoly = inla.spde.make.A(mesh_barrier_transformed, loc = coordinates(pixels_mesh))

saveRDS(proj,'proj_dup_removed.rds') 
saveRDS(proj_barrier,'proj_barrier_dup_removed.rds') 
saveRDS(pixels_plotting, 'pixels.rds')
saveRDS(proj_meshpoly,'proj_meshpoly.rds') 
saveRDS(proj_barrier_meshpoly,'proj_barrier_meshpoly.rds') 
saveRDS(pixels_mesh, 'pixels_meshpoly.rds')
saveRDS(COAST_transformed, 'COAST_transformed.rds')

# Save matrices of search effort for posterior predictive checks
saveRDS(matrix(covariates_pp_J$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F),'E_J_dup_removed.rds')
saveRDS(matrix(covariates_pp_K$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F),'E_K_dup_removed.rds')
saveRDS(matrix(covariates_pp_L$Search_Effort/max_effort, nrow = mesh$n, ncol = no_T, byrow = F),'E_L_dup_removed.rds')

# Save matrices of search effort for posterior predictive checks
saveRDS(matrix(covariates_pp_barrier_J$Search_Effort/max_effort_barrier, nrow = mesh_barrier$n, ncol = no_T, byrow = F),'E_J_dup_removed_barrier.rds')
saveRDS(matrix(covariates_pp_barrier_K$Search_Effort/max_effort_barrier, nrow = mesh_barrier$n, ncol = no_T, byrow = F),'E_K_dup_removed_barrier.rds')
saveRDS(matrix(covariates_pp_barrier_L$Search_Effort/max_effort_barrier, nrow = mesh_barrier$n, ncol = no_T, byrow = F),'E_L_dup_removed_barrier.rds')

#### Add the constrained files
# create the projector matrix for each pod
A_proj_constr_spacetime_J = inla.spde.make.A(mesh, coordinates(total_sightings_J),
                                             group = total_sightings_J$MONTH_INLA,
                                             n.group = no_T)
A_proj_constr_spacetime_barrier_J = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_J),
                                                     group = total_sightings_J$MONTH_INLA,
                                                     n.group = no_T)
imat_constr_spacetime_J = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                                           group = rep(1:no_T,each = mesh$n),
                                           n.group = no_T)
imat_constr_spacetime_barrier_J = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                                                   group = rep(1:no_T,each = mesh_barrier$n),
                                                   n.group = no_T)

A_proj_constr_spacetime_K = inla.spde.make.A(mesh, coordinates(total_sightings_K),
                                             group = total_sightings_K$MONTH_INLA,
                                             n.group = no_T)
A_proj_constr_spacetime_barrier_K = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_K),
                                                     group = total_sightings_K$MONTH_INLA,
                                                     n.group = no_T)
imat_constr_spacetime_K = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                                           group = rep(1:no_T,each = mesh$n),
                                           n.group = no_T)
imat_constr_spacetime_barrier_K = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                                                   group = rep(1:no_T,each = mesh_barrier$n),
                                                   n.group = no_T)

A_proj_constr_spacetime_L = inla.spde.make.A(mesh, coordinates(total_sightings_L),
                                             group = total_sightings_L$MONTH_INLA,
                                             n.group = no_T)
A_proj_constr_spacetime_barrier_L = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_L),
                                                     group = total_sightings_L$MONTH_INLA,
                                                     n.group = no_T)
imat_constr_spacetime_L = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)),
                                           group = rep(1:no_T,each = mesh$n),
                                           n.group = no_T)
imat_constr_spacetime_barrier_L = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)),
                                                   group = rep(1:no_T,each = mesh_barrier$n),
                                                   n.group = no_T)

simple.model.constr.spacetime = inla.spde2.pcmatern(mesh_transformed,
                                                    prior.range = c(15,0.01),
                                                    prior.sigma = c(3,0.1),
                                                    n.iid.group = 6)
simple.model.constr.spacetime.barrier = inla.spde2.pcmatern(mesh_barrier_transformed,
                                                            prior.range = c(15,0.01),
                                                            prior.sigma = c(3,0.1),
                                                            n.iid.group = 6)
s_index_constr_spacetime = inla.spde.make.index('simple.field.constr.spacetime', n.spde = simple.model.constr.spacetime$n.spde)
s_index_constr_spacetime_barrier = inla.spde.make.index('simple.field.constr.spacetime.barrier', n.spde = simple.model.constr.spacetime.barrier$n.spde)

A_proj_constr_J = inla.spde.make.A(mesh, coordinates(total_sightings_J))

A_proj_constr_barrier_J = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_J))

imat_constr_J = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)))

imat_constr_barrier_J = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)))

A_proj_constr_K = inla.spde.make.A(mesh, coordinates(total_sightings_K))

A_proj_constr_barrier_K = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_K))

imat_constr_K = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)))

imat_constr_barrier_K = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)))

A_proj_constr_L = inla.spde.make.A(mesh, coordinates(total_sightings_L))

A_proj_constr_barrier_L = inla.spde.make.A(mesh_barrier, coordinates(total_sightings_L))

imat_constr_L = inla.spde.make.A(mesh, cbind(rep(mesh$loc[,1],times = no_T),rep(mesh$loc[,2],times = no_T)))

imat_constr_barrier_L = inla.spde.make.A(mesh_barrier, cbind(rep(mesh_barrier$loc[,1],times = no_T),rep(mesh_barrier$loc[,2],times = no_T)))

simple.model.constr = inla.spde2.pcmatern(mesh_transformed,
                                          prior.range = c(15,0.01),
                                          prior.sigma = c(3,0.1),
                                          n.iid.group = 1)
simple.model.constr.barrier = inla.spde2.pcmatern(mesh_barrier_transformed,
                                                  prior.range = c(15,0.01),
                                                  prior.sigma = c(3,0.1),
                                                  n.iid.group = 1)
s_index_constr = inla.spde.make.index('simple.field.constr', n.spde = simple.model.constr$n.spde)
s_index_constr_barrier = inla.spde.make.index('simple.field.constr.barrier', n.spde = simple.model.constr.barrier$n.spde)

# join the Identity matrix for the mesh nodes with the projection matrix for the obs
A.pp_constr_J <- rBind(imat_constr_J, A_proj_constr_J)
A.pp_constr_barrier_J <- rBind(imat_constr_barrier_J, A_proj_constr_barrier_J)

A.pp_constr_K <- rBind(imat_constr_K, A_proj_constr_K)
A.pp_constr_barrier_K <- rBind(imat_constr_barrier_K, A_proj_constr_barrier_K)

A.pp_constr_L <- rBind(imat_constr_L, A_proj_constr_L)
A.pp_constr_barrier_L <- rBind(imat_constr_barrier_L, A_proj_constr_barrier_L)

PP_stack_constr_J = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                               A=list(1,A.pp_constr_J), tag='pp',
                               effects=list(covariates_PP_join_J, 
                                            s_index_constr)) # i are the indices for the SPDE 
print('stack pp J complete')

PP_stack_constr_barrier_J = inla.stack(data=list(y=y.pp_barrier_J, e=covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier),#(e.pp_barrier_J/max(e.pp_barrier_J))*(covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier) ),
                                       A=list(1,A.pp_constr_barrier_J), tag='pp',
                                       effects=list(covariates_PP_join_barrier_J, 
                                                    s_index_constr_barrier)) # i are the indices for the SPDE 
print('stack pp barrier J complete')

PP_stack_constr_K = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                               A=list(1,A.pp_constr_K), tag='pp',
                               effects=list(covariates_PP_join_K, 
                                            s_index_constr)) # i are the indices for the SPDE 
print('stack pp K complete')

PP_stack_constr_barrier_K = inla.stack(data=list(y=y.pp_barrier_K, e=covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                       A=list(1,A.pp_constr_barrier_K), tag='pp',
                                       effects=list(covariates_PP_join_barrier_K, 
                                                    s_index_constr_barrier)) # i are the indices for the SPDE 
print('stack pp barrier K complete')

PP_stack_constr_L = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
                               A=list(1,A.pp_constr_L), tag='pp',
                               effects=list(covariates_PP_join_L, 
                                            s_index_constr)) # i are the indices for the SPDE 
print('stack pp L complete')

PP_stack_constr_barrier_L = inla.stack(data=list(y=y.pp_barrier_L, e=covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier),#(e.pp_barrier_L/max(e.pp_barrier_L))*(covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier) ),
                                       A=list(1,A.pp_constr_barrier_L), tag='pp',
                                       effects=list(covariates_PP_join_barrier_L, 
                                                    s_index_constr_barrier)) # i are the indices for the SPDE 
print('stack pp barrier L complete')

# ## Need to recreate depth variable due to error
# PP_stack_constr_J$effects$data$depth = PP_stack_constr_J$effects$data$topo
# #PP_stack_constr_J$effects$data$depth[PP_stack_constr_J$effects$data$depth>=51] = 51
# PP_stack_constr_J$effects$data$depth = log(abs(PP_stack_constr_J$effects$data$depth - 52))
# hist(PP_stack_constr_J$effects$data$depth[PP_stack_constr_J$effects$data$depth >0] )
# PP_stack_constr_J$effects$data$depth = scale(PP_stack_constr_J$effects$data$depth)
# hist(PP_stack_constr_J$effects$data$depth)
# 
# PP_stack_constr_K$effects$data$depth = PP_stack_constr_K$effects$data$topo
# #PP_stack_constr_K$effects$data$depth[PP_stack_constr_K$effects$data$depth>=51] = 51
# PP_stack_constr_K$effects$data$depth = log(abs(PP_stack_constr_K$effects$data$depth - 52))
# hist(PP_stack_constr_K$effects$data$depth[PP_stack_constr_K$effects$data$depth >0] )
# PP_stack_constr_K$effects$data$depth = scale(PP_stack_constr_K$effects$data$depth)
# 
# PP_stack_constr_L$effects$data$depth = PP_stack_constr_L$effects$data$topo
# #PP_stack_constr_L$effects$data$depth[PP_stack_constr_L$effects$data$depth>=51] = 51
# PP_stack_constr_L$effects$data$depth = log(abs(PP_stack_constr_L$effects$data$depth - 52))
# hist(PP_stack_constr_L$effects$data$depth[PP_stack_constr_L$effects$data$depth >0] )
# PP_stack_constr_L$effects$data$depth = scale(PP_stack_constr_L$effects$data$depth)
# 
# PP_stack_constr_barrier_J$effects$data$depth = PP_stack_constr_barrier_J$effects$data$topo
# #PP_stack_constr_barrier_J$effects$data$depth[PP_stack_constr_barrier_J$effects$data$depth>=51] = 51
# PP_stack_constr_barrier_J$effects$data$depth = log(abs(PP_stack_constr_barrier_J$effects$data$depth - 52))
# hist(PP_stack_constr_barrier_J$effects$data$depth[PP_stack_constr_barrier_J$effects$data$depth >0] )
# PP_stack_constr_barrier_J$effects$data$depth = scale(PP_stack_constr_barrier_J$effects$data$depth)
# 
# PP_stack_constr_barrier_K$effects$data$depth = PP_stack_constr_barrier_K$effects$data$topo
# #PP_stack_constr_barrier_K$effects$data$depth[PP_stack_constr_barrier_K$effects$data$depth>=51] = 51
# PP_stack_constr_barrier_K$effects$data$depth = log(abs(PP_stack_constr_barrier_K$effects$data$depth - 52))
# hist(PP_stack_constr_barrier_K$effects$data$depth[PP_stack_constr_barrier_K$effects$data$depth >0] )
# PP_stack_constr_barrier_K$effects$data$depth = scale(PP_stack_constr_barrier_K$effects$data$depth)
# 
# PP_stack_constr_barrier_L$effects$data$depth = PP_stack_constr_barrier_L$effects$data$topo
# #PP_stack_constr_barrier_L$effects$data$depth[PP_stack_constr_barrier_L$effects$data$depth>=51] = 51
# PP_stack_constr_barrier_L$effects$data$depth = log(abs(PP_stack_constr_barrier_L$effects$data$depth - 52))
# hist(PP_stack_constr_barrier_L$effects$data$depth[PP_stack_constr_barrier_L$effects$data$depth >0] )
# PP_stack_constr_barrier_L$effects$data$depth = scale(PP_stack_constr_barrier_L$effects$data$depth)

A.pp_constr_spacetime_J <- rBind(imat_constr_spacetime_J, A_proj_constr_spacetime_J)
A.pp_constr_spacetime_barrier_J <- rBind(imat_constr_spacetime_barrier_J, A_proj_constr_spacetime_barrier_J)

A.pp_constr_spacetime_K <- rBind(imat_constr_spacetime_K, A_proj_constr_spacetime_K)
A.pp_constr_spacetime_barrier_K <- rBind(imat_constr_spacetime_barrier_K, A_proj_constr_spacetime_barrier_K)

A.pp_constr_spacetime_L <- rBind(imat_constr_spacetime_L, A_proj_constr_spacetime_L)
A.pp_constr_spacetime_barrier_L <- rBind(imat_constr_spacetime_barrier_L, A_proj_constr_spacetime_barrier_L)

PP_stack_constr_spacetime_J = inla.stack(data=list(y=y.pp_J, e=covariates_PP_join_J$Search_Effort/max_effort),#(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
                                         A=list(1,A.pp_constr_spacetime_J), tag='pp',
                                         effects=list(covariates_PP_join_J, 
                                                      s_index_constr_spacetime)) # i are the indices for the SPDE 
print('stack pp J complete')

PP_stack_constr_spacetime_barrier_J = inla.stack(data=list(y=y.pp_barrier_J, e=covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier),#(e.pp_barrier_J/max(e.pp_barrier_J))*(covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier) ),
                                                 A=list(1,A.pp_constr_spacetime_barrier_J), tag='pp',
                                                 effects=list(covariates_PP_join_barrier_J, 
                                                              s_index_constr_spacetime_barrier)) # i are the indices for the SPDE 
print('stack pp barrier J complete')

PP_stack_constr_spacetime_K = inla.stack(data=list(y=y.pp_K, e=covariates_PP_join_K$Search_Effort/max_effort),#(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
                                         A=list(1,A.pp_constr_spacetime_K), tag='pp',
                                         effects=list(covariates_PP_join_K, 
                                                      s_index_constr_spacetime)) # i are the indices for the SPDE 
print('stack pp K complete')

PP_stack_constr_spacetime_barrier_K = inla.stack(data=list(y=y.pp_barrier_K, e=covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier),#(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
                                                 A=list(1,A.pp_constr_spacetime_barrier_K), tag='pp',
                                                 effects=list(covariates_PP_join_barrier_K, 
                                                              s_index_constr_spacetime_barrier)) # i are the indices for the SPDE 
print('stack pp barrier K complete')

PP_stack_constr_spacetime_L = inla.stack(data=list(y=y.pp_L, e=covariates_PP_join_L$Search_Effort/max_effort),#(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
                                         A=list(1,A.pp_constr_spacetime_L), tag='pp',
                                         effects=list(covariates_PP_join_L, 
                                                      s_index_constr_spacetime)) # i are the indices for the SPDE 
print('stack pp L complete')

PP_stack_constr_spacetime_barrier_L = inla.stack(data=list(y=y.pp_barrier_L, e=covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier),#(e.pp_barrier_L/max(e.pp_barrier_L))*(covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier) ),
                                                 A=list(1,A.pp_constr_spacetime_barrier_L), tag='pp',
                                                 effects=list(covariates_PP_join_barrier_L, 
                                                              s_index_constr_spacetime_barrier)) # i are the indices for the SPDE 
print('stack pp barrier L complete')

# ## Need to recreate depth variable due to error
# PP_stack_constr_spacetime_J$effects$data$depth = PP_stack_constr_spacetime_J$effects$data$topo
# #PP_stack_constr_spacetime_J$effects$data$depth[PP_stack_constr_spacetime_J$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_J$effects$data$depth = log(abs(PP_stack_constr_spacetime_J$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_J$effects$data$depth[PP_stack_constr_spacetime_J$effects$data$depth >0] )
# PP_stack_constr_spacetime_J$effects$data$depth = scale(PP_stack_constr_spacetime_J$effects$data$depth)
# hist(PP_stack_constr_spacetime_J$effects$data$depth)
# 
# PP_stack_constr_spacetime_K$effects$data$depth = PP_stack_constr_spacetime_K$effects$data$topo
# #PP_stack_constr_spacetime_K$effects$data$depth[PP_stack_constr_spacetime_K$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_K$effects$data$depth = log(abs(PP_stack_constr_spacetime_K$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_K$effects$data$depth[PP_stack_constr_spacetime_K$effects$data$depth >0] )
# PP_stack_constr_spacetime_K$effects$data$depth = scale(PP_stack_constr_spacetime_K$effects$data$depth)
# 
# PP_stack_constr_spacetime_L$effects$data$depth = PP_stack_constr_spacetime_L$effects$data$topo
# #PP_stack_constr_spacetime_L$effects$data$depth[PP_stack_constr_spacetime_L$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_L$effects$data$depth = log(abs(PP_stack_constr_spacetime_L$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_L$effects$data$depth[PP_stack_constr_spacetime_L$effects$data$depth >0] )
# PP_stack_constr_spacetime_L$effects$data$depth = scale(PP_stack_constr_spacetime_L$effects$data$depth)
# 
# PP_stack_constr_spacetime_barrier_J$effects$data$depth = PP_stack_constr_spacetime_barrier_J$effects$data$topo
# #PP_stack_constr_spacetime_barrier_J$effects$data$depth[PP_stack_constr_spacetime_barrier_J$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_barrier_J$effects$data$depth = log(abs(PP_stack_constr_spacetime_barrier_J$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_barrier_J$effects$data$depth[PP_stack_constr_spacetime_barrier_J$effects$data$depth >0] )
# PP_stack_constr_spacetime_barrier_J$effects$data$depth = scale(PP_stack_constr_spacetime_barrier_J$effects$data$depth)
# 
# PP_stack_constr_spacetime_barrier_K$effects$data$depth = PP_stack_constr_spacetime_barrier_K$effects$data$topo
# #PP_stack_constr_spacetime_barrier_K$effects$data$depth[PP_stack_constr_spacetime_barrier_K$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_barrier_K$effects$data$depth = log(abs(PP_stack_constr_spacetime_barrier_K$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_barrier_K$effects$data$depth[PP_stack_constr_spacetime_barrier_K$effects$data$depth >0] )
# PP_stack_constr_spacetime_barrier_K$effects$data$depth = scale(PP_stack_constr_spacetime_barrier_K$effects$data$depth)
# 
# PP_stack_constr_spacetime_barrier_L$effects$data$depth = PP_stack_constr_spacetime_barrier_L$effects$data$topo
# #PP_stack_constr_spacetime_barrier_L$effects$data$depth[PP_stack_constr_spacetime_barrier_L$effects$data$depth>=51] = 51
# PP_stack_constr_spacetime_barrier_L$effects$data$depth = log(abs(PP_stack_constr_spacetime_barrier_L$effects$data$depth - 52))
# hist(PP_stack_constr_spacetime_barrier_L$effects$data$depth[PP_stack_constr_spacetime_barrier_L$effects$data$depth >0] )
# PP_stack_constr_spacetime_barrier_L$effects$data$depth = scale(PP_stack_constr_spacetime_barrier_L$effects$data$depth)

# Join sighting covariates together
X_sightings = rbind(covariates_site_J[,temp_names][,c(3,c(7:15))],
                    covariates_site_K[,temp_names][,c(3,c(7:15))],
                    covariates_site_L[,temp_names][,c(3,c(7:15))])
#X_sightings$depth = log(abs(X_sightings$topo - 52))

model_fitting_files_new_constr = list(PP_stack_constr_J = PP_stack_constr_J,
                                      PP_stack_constr_K = PP_stack_constr_K,
                                      PP_stack_constr_L = PP_stack_constr_L,
                                      PP_stack_constr_barrier_J = PP_stack_constr_barrier_J,
                                      PP_stack_constr_barrier_K = PP_stack_constr_barrier_K,
                                      PP_stack_constr_barrier_L = PP_stack_constr_barrier_L,
                                      PP_stack_constr_spacetime_J = PP_stack_constr_spacetime_J,
                                      PP_stack_constr_spacetime_K = PP_stack_constr_spacetime_K,
                                      PP_stack_constr_spacetime_L = PP_stack_constr_spacetime_L,
                                      PP_stack_constr_spacetime_barrier_J = PP_stack_constr_spacetime_barrier_J,
                                      PP_stack_constr_spacetime_barrier_K = PP_stack_constr_spacetime_barrier_K,
                                      PP_stack_constr_spacetime_barrier_L = PP_stack_constr_spacetime_barrier_L,
                                      A_constr_sightings = rbind(A_proj_constr_J, A_proj_constr_K, A_proj_constr_L),
                                      A_constr_spacetime_sightings = rbind(A_proj_constr_spacetime_J, A_proj_constr_spacetime_K, A_proj_constr_spacetime_L),
                                      A_constr_sightings_barrier = rbind(A_proj_constr_barrier_J, A_proj_constr_barrier_K, A_proj_constr_barrier_L),
                                      A_constr_spacetime_sightings_barrier = rbind(A_proj_constr_spacetime_barrier_J, A_proj_constr_spacetime_barrier_K, A_proj_constr_spacetime_barrier_L),
                                      X_sightings = X_sightings,
                                      barrier.model = barrier.model,
                                      simple.model = simple.model)
#saveRDS(model_fitting_files_new_constr, 'Model_fitting_files_scaled_new_constr.rds')
saveRDS(model_fitting_files_new_constr, 'Model_fitting_files_scaled_new_constr_dup_removed.rds')

# Plot the WW sightings per month
May_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 5),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings May') +
  xlab('Eastings') + ylab('Northings')
June_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 6),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings June') +
  xlab('Eastings') + ylab('Northings')
July_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 7),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings July') +
  xlab('Eastings') + ylab('Northings')
August_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 8),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings August') +
  xlab('Eastings') + ylab('Northings')
September_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 9),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings September') +
  xlab('Eastings') + ylab('Northings')
October_WW_plot = ggplot() + gg(WW_sightings2[which(WW_sightings2$MONTH == 10),]) +
  gg(COAST_plotting) + ggtitle('Initial WW Sightings October') +
  xlab('Eastings') + ylab('Northings')
multiplot(May_WW_plot,June_WW_plot,July_WW_plot,August_WW_plot,September_WW_plot,October_WW_plot ,cols = 2)

## Check the chances of observing all three pods for each month
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 5], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 5]),
        FUN = function(x){sum(x == 3)})/31)
# 0.016

hist(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 6], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 6]),
        FUN = function(x){sum(x == 3)}), 
     xlab = 'Number of days', 
     main = 'Number of days in June in which all three pods are sighted per month/year')
# estimate probability of hitting bound
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 6], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 6]),
        FUN = function(x){sum(x == 3)})/30)
# 0.0292

hist(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 7], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 7]),
        FUN = function(x){sum(x == 3)}), 
     xlab = 'Number of days', 
     main = 'Number of days in July in which all three pods are sighted per month/year')
# estimate probability of hitting bound
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 7], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 7]),
        FUN = function(x){sum(x == 3)})/31)
# 0.234

hist(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 8], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 8]),
        FUN = function(x){sum(x == 3)}), 
     xlab = 'Number of days', 
     main = 'Number of days in August in which all three pods are sighted per month/year')
# estimate probability of hitting bound
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 8], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 8]),
        FUN = function(x){sum(x == 3)})/31)
# 0.282

hist(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 9], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 9]),
        FUN = function(x){sum(x == 3)}), 
     xlab = 'Number of days', 
     main = 'Number of days in September in which all three pods are sighted per month/year')
# estimate probability of hitting bound
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 9], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 9]),
        FUN = function(x){sum(x == 3)})/30)
# 0.25

hist(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 10], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 10]),
        FUN = function(x){sum(x == 3)}), 
     xlab = 'Number of days', 
     main = 'Number of days in October in which all three pods are sighted per month/year')
# estimate probability of hitting bound
mean(by(WW_sightings_identified2$ALLpods[WW_sightings_identified2$MONTH == 10], 
        list(Year=WW_sightings_identified2$YEAR[WW_sightings_identified2$MONTH == 10]),
        FUN = function(x){sum(x == 3)})/31)
# 0.09

# Check the bigger area polygons have higher intensity
plot.field(log(SI_complete_monthly_J_mean[,1]+1), poly = COAST_transformed, mesh = mesh_transformed, 
           pixels = pixels_plotting, corr = FALSE) + colsc(log(as.matrix(SI_complete_monthly_J_mean)+1)) + 
  ggtitle('May WW Search Effort J pod') + 
  gg(spTransform(dmesh,CRSobj = pixels_plotting@proj4string))
# good

## NOTE IMPORTANT!! THE pixels_plotting will cause issues when plotting posterior estimates
## It has been defined on a different polygon used to define the mesh
## Need to use COAST_simp as a mask as this is what was used for defining computational mesh
## Otherwise we show plots on land - these are erroneously estimated to have high intensity.

# create mapping functions to map observations to the dual mesh pixels
# Use this for leave-one-out-pixel cross validation

# find the pixel corresponding to each observation
library(rgeos)
which_dmesh_pixel_J_barrier = gWithin(total_sightings_J,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_pixel_J_barrier = unlist(which_dmesh_pixel_J_barrier)
which_dmesh_pixel_J_barrier = as.numeric(which_dmesh_pixel_J_barrier)

which_dmesh_pixel_K_barrier = gWithin(total_sightings_K,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_pixel_K_barrier = unlist(which_dmesh_pixel_K_barrier)
which_dmesh_pixel_K_barrier = as.numeric(which_dmesh_pixel_K_barrier)

which_dmesh_pixel_L_barrier = gWithin(total_sightings_L,dmesh_barrier, byid = T, returnDense = F)
which_dmesh_pixel_L_barrier = unlist(which_dmesh_pixel_L_barrier)
which_dmesh_pixel_L_barrier = as.numeric(which_dmesh_pixel_L_barrier)

# assess search effort for pixels that contain a point
which_dmesh_pixel_J_month = which_dmesh_pixel_J + ( (total_sightings_J$MONTH_INLA - 1) * mesh$n )
which_dmesh_pixel_K_month = which_dmesh_pixel_K + ( (total_sightings_K$MONTH_INLA - 1) * mesh$n )
which_dmesh_pixel_L_month = which_dmesh_pixel_L + ( (total_sightings_L$MONTH_INLA - 1) * mesh$n )

which_dmesh_pixel_J_barrier_month = which_dmesh_pixel_J_barrier + ( (total_sightings_J$MONTH_INLA - 1) * mesh_barrier$n )
which_dmesh_pixel_K_barrier_month = which_dmesh_pixel_K_barrier + ( (total_sightings_K$MONTH_INLA - 1) * mesh_barrier$n )
which_dmesh_pixel_L_barrier_month = which_dmesh_pixel_L_barrier + ( (total_sightings_L$MONTH_INLA - 1) * mesh_barrier$n )

summary(as.numeric(SI_complete_monthly_J_mean)[which_dmesh_pixel_J_month])
summary(as.numeric(SI_complete_monthly_K_mean)[which_dmesh_pixel_K_month])
summary(as.numeric(SI_complete_monthly_L_mean)[which_dmesh_pixel_L_month])

hist(as.numeric(SI_complete_J_mean), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_J_mean)[which_dmesh_pixel_J_month], freq=F, add=T, border='red')
hist(as.numeric(SI_complete_K_mean), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_K_mean)[which_dmesh_pixel_K_month], freq=F, add=T, border='red')
hist(as.numeric(SI_complete_L_mean), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_L_mean)[which_dmesh_pixel_L_month], freq=F, add=T, border='red')

summary(as.numeric(SI_complete_monthly_J_mean_barrier)[which_dmesh_pixel_J_barrier_month])
summary(as.numeric(SI_complete_monthly_K_mean_barrier)[which_dmesh_pixel_K_barrier_month])
summary(as.numeric(SI_complete_monthly_L_mean_barrier)[which_dmesh_pixel_L_barrier_month])

hist(as.numeric(SI_complete_J_mean_barrier), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_J_mean_barrier)[which_dmesh_pixel_J_barrier_month], freq=F, add=T, border='red')
hist(as.numeric(SI_complete_K_mean_barrier), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_K_mean_barrier)[which_dmesh_pixel_K_barrier_month], freq=F, add=T, border='red', breaks=30)
hist(as.numeric(SI_complete_L_mean_barrier), freq = F, breaks=100, xlim = c(0,250))
hist(as.numeric(SI_complete_monthly_L_mean_barrier)[which_dmesh_pixel_L_barrier_month], freq=F, add=T, border='red', breaks=30)

# How many points fall in each polygon
pixel_counts_J_month = as.numeric( 1:(mesh$n*6) %in% which_dmesh_pixel_J_month )
pixel_counts_J_month[pixel_counts_J_month != 0] = as.numeric(table(pmatch(which_dmesh_pixel_J_month,1:(mesh$n*6), duplicates.ok = T)))
pixel_counts_J_month = matrix(pixel_counts_J_month, nrow=mesh$n, ncol = 6, byrow = F)

pixel_counts_K_month = as.numeric( 1:(mesh$n*6) %in% which_dmesh_pixel_K_month )
pixel_counts_K_month[pixel_counts_K_month != 0] = as.numeric(table(pmatch(which_dmesh_pixel_K_month,1:(mesh$n*6), duplicates.ok = T)))
pixel_counts_K_month = matrix(pixel_counts_K_month, nrow=mesh$n, ncol = 6, byrow = F)

pixel_counts_L_month = as.numeric( 1:(mesh$n*6) %in% which_dmesh_pixel_L_month )
pixel_counts_L_month[pixel_counts_L_month != 0] = as.numeric(table(pmatch(which_dmesh_pixel_L_month,1:(mesh$n*6), duplicates.ok = T)))
pixel_counts_L_month = matrix(pixel_counts_L_month, nrow=mesh$n, ncol = 6, byrow = F)

pixel_counts = list(pixel_counts_J_month = pixel_counts_J_month,
                    pixel_counts_K_month = pixel_counts_K_month,
                    pixel_counts_L_month = pixel_counts_L_month)
saveRDS(pixel_counts, 'pixel_counts_monthly_dup_removed.rds')

table(pmatch(which_dmesh_pixel_J_month,1:(mesh$n*6), duplicates.ok = T))
table(pmatch(1:(mesh$n*6), which_dmesh_pixel_J_month, duplicates.ok = T))
# How many pixels have k points in them
table(table(pmatch(which_dmesh_pixel_J_month,1:(mesh$n*6), duplicates.ok = T)))
table(table(pmatch(which_dmesh_pixel_K_month,1:(mesh$n*6), duplicates.ok = T)))
table(table(pmatch(which_dmesh_pixel_L_month,1:(mesh$n*6), duplicates.ok = T)))

table(table(pmatch(which_dmesh_pixel_J_barrier_month,1:(mesh_barrier$n*6), duplicates.ok = T)))
table(table(pmatch(which_dmesh_pixel_K_barrier_month,1:(mesh_barrier$n*6), duplicates.ok = T)))
table(table(pmatch(which_dmesh_pixel_L_barrier_month,1:(mesh_barrier$n*6), duplicates.ok = T)))

# Use spatstat package for exploratory analysis
library(spatstat)
library(maptools)

ppp_J = ppp(x = coordinates(total_sightings_J)[,1],
            y = coordinates(total_sightings_J)[,2],
            window = as.owin(W = COAST_simp))#,
            #marks = covariates_site_J)
plot(ppp_J, add=F, use.marks=F, pch=15)
ppp_J = rescale(ppp_J, s = 1000, unitname = 'km')
ppp_J
sum(duplicated(total_sightings_J@coords))

ppp_K = ppp(x = coordinates(total_sightings_K)[,1],
            y = coordinates(total_sightings_K)[,2],
            window = as.owin(W = COAST_simp))#,
            #marks = covariates_site_K)
plot(ppp_K, add=F, use.marks=F, pch=15)
ppp_K = rescale(ppp_K, s = 1000, unitname = 'km')
sum(duplicated(total_sightings_K@coords))

ppp_L = ppp(x = coordinates(total_sightings_L)[,1],
            y = coordinates(total_sightings_L)[,2],
            window = as.owin(W = COAST_simp))#,
            #marks = covariates_site_L)
plot(ppp_L, add=F, use.marks=F, pch=15)
ppp_L = rescale(ppp_L, s = 1000, unitname = 'km')
sum(duplicated(total_sightings_L@coords))

ppp_total = superimpose(ppp_J,ppp_K,ppp_L)

ppp_J = ppp(x = coordinates(total_sightings_J)[,1],
            y = coordinates(total_sightings_J)[,2],
            window = as.owin(W = COAST_simp),
            marks = covariates_site_J)
ppp_J = rescale(ppp_J, s = 1000, unitname = 'km')

ppp_K = ppp(x = coordinates(total_sightings_K)[,1],
            y = coordinates(total_sightings_K)[,2],
            window = as.owin(W = COAST_simp),
            marks = covariates_site_K)
ppp_K = rescale(ppp_K, s = 1000, unitname = 'km')

ppp_L = ppp(x = coordinates(total_sightings_L)[,1],
            y = coordinates(total_sightings_L)[,2],
            window = as.owin(W = COAST_simp),
            marks = covariates_site_L)
ppp_L = rescale(ppp_L, s = 1000, unitname = 'km')

# create a function that creates an 'im' file for weighting search effort by
im.maker = function(field, poly, mesh, pixels, corr = TRUE,...){
  stopifnot(length(field) == mesh$n)
  
  # - choose plotting region to be the same as the study area polygon
  #pixels_plotting = pixels(mesh, mask = poly, nx = 300, ny = 300)
  proj = inla.spde.make.A(mesh, loc = coordinates(pixels) ) 
  #proj = inla.mesh.projector(mesh, xlim = xlim, 
  #                           ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = as.numeric(proj %*% field)
  
  field_df = SpatialPixelsDataFrame(pixels, 
                                    data = data.frame(field = field.proj))
  field_df = as( as( field_df ,"SpatialGridDataFrame"),'im')
  return(field_df)
}

May.im = im.maker(SI_complete_monthly_J_mean[,1] +
                      SI_complete_monthly_K_mean[,1] +
                      SI_complete_monthly_L_mean[,1], poly = COAST_transformed, mesh = mesh_transformed, 
           pixels = pixels_plotting, corr = FALSE)
May.im$v[is.na(May.im$v)] = 0.5
May.im$v[May.im$v < 0.5] = 0.5

SST.im.May = im.maker(SST_dmesh_month$weighted_mean$May,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.May[is.na(SST.im.May$v)] = mean(SST.im.May$v)
chloro.im.May = im.maker(chloro_dmesh_month$weighted_mean$May,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.May[is.na(chloro.im.May$v)] = mean(chloro.im.May$v)

June.im = im.maker(SI_complete_monthly_J_mean[,2] +
                    SI_complete_monthly_K_mean[,2] +
                    SI_complete_monthly_L_mean[,2], poly = COAST_transformed, mesh = mesh_transformed, 
                  pixels = pixels_plotting, corr = FALSE)
June.im$v[is.na(June.im$v)] = 0.5
June.im$v[June.im$v < 0.5] = 0.5

SST.im.June = im.maker(SST_dmesh_month$weighted_mean$June,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.June[is.na(SST.im.June$v)] = mean(SST.im.June$v)
chloro.im.June = im.maker(chloro_dmesh_month$weighted_mean$June,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.June[is.na(chloro.im.June$v)] = mean(chloro.im.June$v)

July.im = im.maker(SI_complete_monthly_J_mean[,3] +
                    SI_complete_monthly_K_mean[,3] +
                    SI_complete_monthly_L_mean[,3], poly = COAST_transformed, mesh = mesh_transformed, 
                  pixels = pixels_plotting, corr = FALSE)
July.im$v[is.na(July.im$v)] = 0.5
July.im$v[July.im$v < 0.5] = 0.5

SST.im.July = im.maker(SST_dmesh_month$weighted_mean$July,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.July[is.na(SST.im.July$v)] = mean(SST.im.July$v)
chloro.im.July = im.maker(chloro_dmesh_month$weighted_mean$July,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.July[is.na(chloro.im.July$v)] = mean(chloro.im.July$v)

August.im = im.maker(SI_complete_monthly_J_mean[,4] +
                    SI_complete_monthly_K_mean[,4] +
                    SI_complete_monthly_L_mean[,4], poly = COAST_transformed, mesh = mesh_transformed, 
                  pixels = pixels_plotting, corr = FALSE)
August.im$v[is.na(August.im$v)] = 0.5
August.im$v[August.im$v < 0.5] = 0.5

SST.im.August = im.maker(SST_dmesh_month$weighted_mean$Aug,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.August[is.na(SST.im.August$v)] = mean(SST.im.August$v)
chloro.im.August = im.maker(chloro_dmesh_month$weighted_mean$Aug,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.August[is.na(chloro.im.August$v)] = mean(chloro.im.August$v)

September.im = im.maker(SI_complete_monthly_J_mean[,5] +
                    SI_complete_monthly_K_mean[,5] +
                    SI_complete_monthly_L_mean[,5], poly = COAST_transformed, mesh = mesh_transformed, 
                  pixels = pixels_plotting, corr = FALSE)
September.im$v[is.na(September.im$v)] = 0.5
September.im$v[September.im$v < 0.5] = 0.5

SST.im.September = im.maker(SST_dmesh_month$weighted_mean$Sep,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.September[is.na(SST.im.September$v)] = mean(SST.im.September$v)
chloro.im.September = im.maker(chloro_dmesh_month$weighted_mean$Sep,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.September[is.na(chloro.im.September$v)] = mean(chloro.im.September$v)

October.im = im.maker(SI_complete_monthly_J_mean[,6] +
                    SI_complete_monthly_K_mean[,6] +
                    SI_complete_monthly_L_mean[,6], poly = COAST_transformed, mesh = mesh_transformed, 
                  pixels = pixels_plotting, corr = FALSE)
October.im$v[is.na(October.im$v)] = 0.5
October.im$v[October.im$v < 0.5] = 0.5

SST.im.October = im.maker(SST_dmesh_month$weighted_mean$Oct,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
SST.im.October[is.na(SST.im.October$v)] = mean(SST.im.October$v)
chloro.im.October = im.maker(chloro_dmesh_month$weighted_mean$Oct,poly = COAST_transformed, mesh = mesh_transformed, pixels = pixels_plotting)
chloro.im.October[is.na(chloro.im.October$v)] = mean(chloro.im.October$v)

# test for departures complete spatial randomness for May J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_May = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==1],
              covariate = May.im,#SST.im, 
              #baseline = May.J.im, #horvitz = T,
              subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
              maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_May, xlim = quantile(May.im, probs=c(0.0,0.9)), xlab = 'May Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in May')
pred_May = predict(test_May)
plot(pred_May) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==1])
# Now estimate a kernel density estimate (spatially).
kden_May = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==1],diggle=T)
plot(kden_May)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_May, 
      kden_May)
# No - more is at play
plot(eval.im(kden_May - pred_May)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_May = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==1],
               SST.im.May,
               baseline = May.im,
               bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_May, xlim = quantile(SST.im.May, probs = c(0.1,0.8)),xlab = 'May SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in May')
test3_May = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==1],
               chloro.im.May,
               baseline = May.im,
               bw =1)
# How about chlorophyll?
plot(test3_May, xlim = quantile(chloro.im.May, probs = c(0.1,0.8)),xlab = 'May Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in May')

# test for departures complete spatial randomness for June J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_June = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==2],
                  covariate = June.im,#SST.im, 
                  #baseline = June.J.im, #horvitz = T,
                  subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
                  maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_June, xlim = quantile(June.im, probs=c(0.0,0.9)), xlab = 'June Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in June')
pred_June = predict(test_June)
plot(pred_June) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==2])
# Now estimate a kernel density estimate (spatially).
kden_June = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==2],diggle=T)
plot(kden_June)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_June, kden_June)
# No - more is at play
plot(eval.im(kden_June - pred_June)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_June = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==2],
                   SST.im.June,
                   baseline = June.im,
                   bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_June, xlim = quantile(SST.im.June, probs = c(0.1,0.8)),xlab = 'June SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in June')
test3_June = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==2],
                   chloro.im.June,
                   baseline = June.im,
                   bw =1)
# How about chlorophyll?
plot(test3_June, xlim = quantile(chloro.im.June, probs = c(0.1,0.8)),xlab = 'June Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in June')

# test for departures complete spatial randomness for July J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_July = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==3],
                   covariate = July.im,#SST.im, 
                   #baseline = July.J.im, #horvitz = T,
                   subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
                   maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_July, xlim = quantile(July.im, probs=c(0.0,0.9)), xlab = 'July Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in July')
pred_July = predict(test_July)
plot(pred_July) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==3])
# Now estimate a kernel density estimate (spatially).
kden_July = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==3],diggle=T)
plot(kden_July)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_July,kden_July) 
# No - more is at play
plot(eval.im(kden_July - pred_July)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_July = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==3],
                    SST.im.July,
                    baseline = July.im,
                    bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_July, xlim = quantile(SST.im.July, probs = c(0.1,0.8)),xlab = 'July SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in July')
test3_July = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==3],
                    chloro.im.July,
                    baseline = July.im,
                    bw =1)
# How about chlorophyll?
plot(test3_July, xlim = quantile(chloro.im.July, probs = c(0.1,0.8)),xlab = 'July Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in July')

# test for departures complete spatial randomness for August J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_August = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4],
                     covariate = August.im,#SST.im, 
                     #baseline = August.J.im, #horvitz = T,
                     subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
                     maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_August, xlim = quantile(August.im, probs=c(0.0,0.9)), xlab = 'August Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in August')
pred_August = predict(test_August)
plot(pred_August) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4])
# Now estimate a kernel density estimate (spatially).
kden_August = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4],diggle=T)
plot(kden_August)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_August, 
      density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4],5,diggle=T))
# No - more is at play
plot(eval.im(kden_August - pred_August)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_August = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4],
                      SST.im.August,
                      baseline = August.im,
                      bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_August, xlim = quantile(SST.im.August, probs = c(0.1,0.8)),xlab = 'August SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in August')
test3_August = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==4],
                      chloro.im.August,
                      baseline = August.im,
                      bw =1)
# How about chlorophyll?
plot(test3_August, xlim = quantile(chloro.im.August, probs = c(0.1,0.8)),xlab = 'August Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in August')

# test for departures complete spatial randomness for September J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_September = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5],
                        covariate = September.im,#SST.im, 
                        #baseline = September.J.im, #horvitz = T,
                        subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
                        maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_September, xlim = quantile(September.im, probs=c(0.0,0.9)), xlab = 'September Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in September')
pred_September = predict(test_September)
plot(pred_September) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5])
# Now estimate a kernel density estimate (spatially).
kden_September = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5],5,diggle=T)
plot(kden_September)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_September, 
      density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5],5,diggle=T))
# No - more is at play
plot(eval.im(kden_September - pred_September)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_September = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5],
                         SST.im.September,
                         baseline = September.im,
                         bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_September, xlim = quantile(SST.im.September, probs = c(0.1,0.8)),xlab = 'September SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in September')
test3_September = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==5],
                         chloro.im.September,
                         baseline = September.im,
                         bw =1)
# How about chlorophyll?
plot(test3_September, xlim = quantile(chloro.im.September, probs = c(0.1,0.8)),xlab = 'September Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in September')

# test for departures complete spatial randomness for October J relative to search effort
# first look at the kernel density smooth estimate of the effect of estimated search effort on intensity
# Should be approximately a 1-1 relationship (once effects of covariates removed)
test_October = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6],
                        covariate = October.im,#SST.im, 
                        #baseline = October.J.im, #horvitz = T,
                        subset = as.owin(W = spTransform(COAST_simp, COAST_transformed@proj4string)),
                        maxit = 100)
# Now plot the kernel density estimate of the effect of search effort on mean intensity
# Should be a positive linear relationship. If completely spatially random, conditioned on effort - should be 45 degree line.
plot(test_October, xlim = quantile(October.im, probs=c(0.0,0.9)), xlab = 'October Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in October')
pred_October = predict(test_October)
plot(pred_October) # show predicted intensity given search effort only
plot(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6])
# Now estimate a kernel density estimate (spatially).
kden_October = density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6],5,diggle=T)
plot(kden_October)
# Do the two plots agree? Is there a 45 degree line? This would indicate the search effort alone is enough
pairs(pred_October, 
      density(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6],5,diggle=T))
# No - more is at play
plot(eval.im(kden_October - pred_October)) # Should be a map of values close to 0 (white noise) if no spatial effects remain.
# Now add SST, conditioned on search effort as an offset.
test2_October = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6],
                         SST.im.October,
                         baseline = October.im,
                         bw =0.5)
# Is there any kernel density estimated effect of SST?
plot(test2_October, xlim = quantile(SST.im.October, probs = c(0.1,0.8)),xlab = 'October SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in October')
test3_October = rhohat(ppp_total[c(ppp_J$marks$MONTH_INLA,ppp_K$marks$MONTH_INLA,ppp_L$marks$MONTH_INLA)==6],
                         chloro.im.October,
                         baseline = October.im,
                       bw =1)
# How about chlorophyll?
plot(test3_October, xlim = quantile(chloro.im.October, probs = c(0.1,0.8)),xlab = 'October Chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chlorophyll effect, conditional on search effort in October')

# Plot them together
par(mfrow=c(3,2))
plot(test_May, xlim = quantile(May.im, probs=c(0.0,0.9)), legend=F, xlab = 'May Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in May')
plot(test_June, xlim = quantile(June.im, probs=c(0.0,0.9)), legend=F, xlab = 'June Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in June')
plot(test_July, xlim = quantile(July.im, probs=c(0.0,0.9)), legend=F, xlab = 'July Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in July')
plot(test_August, xlim = quantile(August.im, probs=c(0.0,0.9)), legend=F, xlab = 'August Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in August')
plot(test_September, xlim = quantile(September.im, probs=c(0.0,0.9)), legend=F, xlab = 'September Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in September')
plot(test_October, xlim = quantile(October.im, probs=c(0.0,0.9)), legend=F, xlab = 'October Search Effort', ylab = 'Nonparametric intensity function estimate',
     main='A plot of the kernel density estimate of intensity function with search effort in October')
par(mfrow=c(1,1))

par(mfrow=c(3,2))
plot(test_May, xlim = quantile(May.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in May', xaxt='n', yaxt='n')
plot(test_June, xlim = quantile(June.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in June', xaxt='n', yaxt='n')
plot(test_July, xlim = quantile(July.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in July', xaxt='n', yaxt='n')
plot(test_August, xlim = quantile(August.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in August', xaxt='n', yaxt='n')
plot(test_September, xlim = quantile(September.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in September', xaxt='n', yaxt='n')
plot(test_October, xlim = quantile(October.im, probs=c(0.0,0.9)), legend=F, xlab = '', ylab = '',
     main='A plot of the kernel density estimate of intensity function with search effort in October', xaxt='n', yaxt='n')
par(mfrow=c(1,1))

par(mfrow=c(3,2))
plot(test2_May, legend=F, xlim = quantile(SST.im.May, probs = c(0,0.8)),xlab = 'May SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in May')
plot(test2_June, legend=F, xlim = quantile(SST.im.June, probs = c(0,0.8)),xlab = 'June SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in June')
plot(test2_July, legend=F, xlim = quantile(SST.im.July, probs = c(0,0.8)),xlab = 'July SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in July')
plot(test2_August, legend=F, xlim = quantile(SST.im.August, probs = c(0,0.8)),xlab = 'August SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in August')
plot(test2_September, legend=F, xlim = quantile(SST.im.September, probs = c(0,0.8)),xlab = 'September SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in September')
plot(test2_October, legend=F, xlim = quantile(SST.im.October, probs = c(0,0.8)),xlab = 'October SST', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the sea surface temperature effect, conditional on search effort in October')
par(mfrow=c(1,1))

par(mfrow=c(3,2))
plot(test3_May, legend=F, xlim = quantile(chloro.im.May, probs = c(0.1,0.8)),xlab = 'May chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in May')
plot(test3_June, legend=F, xlim = quantile(chloro.im.June, probs = c(0.1,0.8)),xlab = 'June chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in June')
plot(test3_July, legend=F, xlim = quantile(chloro.im.July, probs = c(0.1,0.8)),xlab = 'July chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in July')
plot(test3_August, legend=F, xlim = quantile(chloro.im.August, probs = c(0.1,0.8)),xlab = 'August chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in August')
plot(test3_September, legend=F, xlim = quantile(chloro.im.September, probs = c(0.1,0.8)),xlab = 'September chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in September')
plot(test3_October, legend=F, xlim = quantile(chloro.im.October, probs = c(0.1,0.8)),xlab = 'October chloropyll-A concentration', ylab = 'Nonparametric intensity function estimate',
     main='Kernel density estimate of the chloropyll-A effect, conditional on search effort in October')
par(mfrow=c(1,1))

# K function
ppp_K = Kest(ppp_total)
plot(ppp_K)
plot(envelope(ppp_total))

total.im = May.im
total.im$v = (1/6)*(May.im$v + June.im$v + July.im$v + August.im$v + September.im$v + October.im$v)
# Now conditioned on search effort (mean across months)
ppp_Kinhom = Kinhom(ppp_total, lambda = total.im, nlarge = 10000)
plot(ppp_Kinhom)
#plot(envelope(ppp_total))

# Save files needed for MCMC model fitting
MCMC_files = list(covariates_pp_J = covariates_pp_J,
                  covariates_pp_K = covariates_pp_K,
                  covariates_pp_L = covariates_pp_L,
                  covariates_site_J = covariates_site_J,
                  covariates_site_K = covariates_site_K,
                  covariates_site_L = covariates_site_L,
                  SI_Brian_monthly_corrected_J = SI_Brian_monthly_J,
                  SI_Brian_monthly_corrected_K = SI_Brian_monthly_K,
                  SI_Brian_monthly_corrected_L = SI_Brian_monthly_L,
                  y.pp_J = y.pp_J,
                  y.pp_K = y.pp_K,
                  y.pp_L = y.pp_L,
                  A.pp_J = A.pp_J,
                  A.pp_K = A.pp_K,
                  A.pp_L = A.pp_L,
                  s_index = s_index,
                  s_index_L = s_index_L,
                  max_effort = max_effort,
                  w = w)
saveRDS(MCMC_files, 'MCMC_files_dup_removed.rds')

# how many sightings per month per pod for model assessment
observed_sightings_monthly2 = as.numeric(rbind(table(total_sightings_J$MONTH), table(total_sightings_K$MONTH), table(total_sightings_L$MONTH)))
saveRDS(observed_sightings_monthly2,'observed_sightings_monthly_dup_removed.rds')

# covariates_PP_join_J$Intercept = 1
# covariates_PP_join_J$Search_Effort_SD_scaled = covariates_PP_join_J$Search_Effort_SD / max_effort
# PP_stack_J = inla.stack(data=list(y=y.pp_J, e=(e.pp_J/max(e.pp_J))*(covariates_PP_join_J$Search_Effort/max_effort) ),
#                       A=list(1,A.pp_constr_J), tag='pp',
#                       effects=list(covariates_PP_join_J, 
#                                    s_index)) # i are the indices for the SPDE 
# print('stack pp J complete')
# 
# covariates_PP_join_barrier_J$Intercept = 1
# covariates_PP_join_barrier_J$Search_Effort_SD_scaled = covariates_PP_join_barrier_J$Search_Effort_SD / max_effort_barrier
# PP_stack_barrier_J = inla.stack(data=list(y=y.pp_barrier_J, e=(e.pp_barrier_J/max(e.pp_barrier_J))*(covariates_PP_join_barrier_J$Search_Effort/max_effort_barrier) ),
#                         A=list(1,A.pp_constr_barrier_J), tag='pp',
#                         effects=list(covariates_PP_join_barrier_J, 
#                                      s_index_barrier)) # i are the indices for the SPDE 
# print('stack pp barrier J complete')
# 
# covariates_PP_join_K$Intercept = 1
# covariates_PP_join_K$Search_Effort_SD_scaled = covariates_PP_join_K$Search_Effort_SD / max_effort
# PP_stack_K = inla.stack(data=list(y=y.pp_K, e=(e.pp_K/max(e.pp_K))*(covariates_PP_join_K$Search_Effort/max_effort) ),
#                         A=list(1,A.pp_constr_K), tag='pp',
#                         effects=list(covariates_PP_join_K, 
#                                      s_index)) # i are the indices for the SPDE 
# print('stack pp K complete')
# 
# covariates_PP_join_barrier_K$Intercept = 1
# covariates_PP_join_barrier_K$Search_Effort_SD_scaled = covariates_PP_join_barrier_K$Search_Effort_SD / max_effort_barrier
# PP_stack_barrier_K = inla.stack(data=list(y=y.pp_barrier_K, e=(e.pp_barrier_K/max(e.pp_barrier_K))*(covariates_PP_join_barrier_K$Search_Effort/max_effort_barrier) ),
#                                 A=list(1,A.pp_constr_barrier_K), tag='pp',
#                                 effects=list(covariates_PP_join_barrier_K, 
#                                              s_index_barrier)) # i are the indices for the SPDE 
# print('stack pp barrier K complete')
# 
# covariates_PP_join_L$Intercept = 1
# covariates_PP_join_L$Search_Effort_SD_scaled = covariates_PP_join_L$Search_Effort_SD / max_effort
# PP_stack_L = inla.stack(data=list(y=y.pp_L, e=(e.pp_L/max(e.pp_L))*(covariates_PP_join_L$Search_Effort/max_effort) ),
#                         A=list(1,A.pp_constr_L), tag='pp',
#                         effects=list(covariates_PP_join_L, 
#                                      s_index)) # i are the indices for the SPDE 
# print('stack pp L complete')
# 
# covariates_PP_join_barrier_L$Intercept = 1
# covariates_PP_join_barrier_L$Search_Effort_SD_scaled = covariates_PP_join_barrier_L$Search_Effort_SD / max_effort_barrier
# PP_stack_barrier_L = inla.stack(data=list(y=y.pp_barrier_L, e=(e.pp_barrier_L/max(e.pp_barrier_L))*(covariates_PP_join_barrier_L$Search_Effort/max_effort_barrier) ),
#                                 A=list(1,A.pp_constr_barrier_L), tag='pp',
#                                 effects=list(covariates_PP_join_barrier_L, 
#                                              s_index_barrier)) # i are the indices for the SPDE 
# print('stack pp barrier L complete')

# save.image('SI_WS4.RData')
# 
# model_fitting_files = list(PP_stack_J = PP_stack_J,
#                            PP_stack_K = PP_stack_K,
#                            PP_stack_L = PP_stack_L,
#                            PP_stack_barrier_J = PP_stack_barrier_J,
#                            PP_stack_barrier_K = PP_stack_barrier_K,
#                            PP_stack_barrier_L = PP_stack_barrier_L,
#                            barrier.model = barrier.model,
#                            simple.model = simple.model)
# 
# 
# saveRDS(model_fitting_files, 'Model_fitting_files_scaled.rds')
# 
# 
# ## Need to recreate depth variable due to error
# PP_stack_J$effects$data$depth = PP_stack_J$effects$data$topo
# #PP_stack_J$effects$data$depth[PP_stack_J$effects$data$depth>=51] = 51
# PP_stack_J$effects$data$depth = log(abs(PP_stack_J$effects$data$depth - 52))
# hist(PP_stack_J$effects$data$depth[PP_stack_J$effects$data$depth >0] )
# PP_stack_J$effects$data$depth = scale(PP_stack_J$effects$data$depth)
# hist(PP_stack_J$effects$data$depth)
# 
# PP_stack_K$effects$data$depth = PP_stack_K$effects$data$topo
# #PP_stack_K$effects$data$depth[PP_stack_K$effects$data$depth>=51] = 51
# PP_stack_K$effects$data$depth = log(abs(PP_stack_K$effects$data$depth - 52))
# hist(PP_stack_K$effects$data$depth[PP_stack_K$effects$data$depth >0] )
# PP_stack_K$effects$data$depth = scale(PP_stack_K$effects$data$depth)
# 
# PP_stack_L$effects$data$depth = PP_stack_L$effects$data$topo
# #PP_stack_L$effects$data$depth[PP_stack_L$effects$data$depth>=51] = 51
# PP_stack_L$effects$data$depth = log(abs(PP_stack_L$effects$data$depth - 52))
# hist(PP_stack_L$effects$data$depth[PP_stack_L$effects$data$depth >0] )
# PP_stack_L$effects$data$depth = scale(PP_stack_L$effects$data$depth)
# 
# PP_stack_barrier_J$effects$data$depth = PP_stack_barrier_J$effects$data$topo
# #PP_stack_barrier_J$effects$data$depth[PP_stack_barrier_J$effects$data$depth>=51] = 51
# PP_stack_barrier_J$effects$data$depth = log(abs(PP_stack_barrier_J$effects$data$depth - 52))
# hist(PP_stack_barrier_J$effects$data$depth[PP_stack_barrier_J$effects$data$depth >0] )
# PP_stack_barrier_J$effects$data$depth = scale(PP_stack_barrier_J$effects$data$depth)
# 
# PP_stack_barrier_K$effects$data$depth = PP_stack_barrier_K$effects$data$topo
# #PP_stack_barrier_K$effects$data$depth[PP_stack_barrier_K$effects$data$depth>=51] = 51
# PP_stack_barrier_K$effects$data$depth = log(abs(PP_stack_barrier_K$effects$data$depth - 52))
# hist(PP_stack_barrier_K$effects$data$depth[PP_stack_barrier_K$effects$data$depth >0] )
# PP_stack_barrier_K$effects$data$depth = scale(PP_stack_barrier_K$effects$data$depth)
# 
# PP_stack_barrier_L$effects$data$depth = PP_stack_barrier_L$effects$data$topo
# #PP_stack_barrier_L$effects$data$depth[PP_stack_barrier_L$effects$data$depth>=51] = 51
# PP_stack_barrier_L$effects$data$depth = log(abs(PP_stack_barrier_L$effects$data$depth - 52))
# hist(PP_stack_barrier_L$effects$data$depth[PP_stack_barrier_L$effects$data$depth >0] )
# PP_stack_barrier_L$effects$data$depth = scale(PP_stack_barrier_L$effects$data$depth)
# 
# model_fitting_files_new = list(PP_stack_J = PP_stack_J,
#                            PP_stack_K = PP_stack_K,
#                            PP_stack_L = PP_stack_L,
#                            PP_stack_barrier_J = PP_stack_barrier_J,
#                            PP_stack_barrier_K = PP_stack_barrier_K,
#                            PP_stack_barrier_L = PP_stack_barrier_L,
#                            barrier.model = barrier.model,
#                            simple.model = simple.model)
# 
# 
# saveRDS(model_fitting_files_new, 'Model_fitting_files_scaled_new.rds')
# 
# # Save the files needed for simulating and plotting
# covariates_pp_J$depth = covariates_pp_J$topo
# #covariates_pp_J$depth[covariates_pp_J$depth>=51] = 51
# covariates_pp_J$depth = log(abs(covariates_pp_J$depth - 52))
# covariates_pp_J$depth = scale(covariates_pp_J$depth)
# 
# hist(covariates_pp_J$depth[covariates_pp_J$depth >0] )
# 
# saveRDS(covariates_pp_J,'covariates_pp_new.rds') # same for inla mesh

####--------


# To change the WW search effort, go to line 3733


#####-------


########### Code for INLAbru below not used


# Below is code for INLAbru - note we do not use this due to lack of barrier model fit

# cmp_simple = coordinates ~ mySmooth(map = coordinates,
#                                     model = simple.model) +
#                            Intercept
# 
# inla.setOption("pardiso.license", "~/licenses/pardiso.lic")
# fit_simple = lgcp(cmp_simple, total_sightings, samplers = COAST_simp,
#                   domain = list(coordinates = mesh),
#                   options = list(#control.fixed = list(expand.factor.strategy = "inla"),
#                     control.compute = list(openmp.strategy = "pardiso.parallel", config=TRUE, dic = T, cpo = T, waic = T, smtp = 'pardiso'), 
#                     control.mode = list(theta = c(10.3932, 1.3800), restart=TRUE ), 
#                     control.inla = list(int.strategy = 'eb'),
#                     verbose = T,
#                     run = T)) # 30 seconds with pardiso parallel

# predict based on 100 posterior samples - change by adding n.samples=///
# To increase the resolution of the plot - add nx and ny arguments to pixels function 
#lambda <- predict(fit_simple, pixels(mesh, mask = COAST_simp, nx = 300, ny = 300), ~ exp(mySmooth + Intercept), n.samples = 1000)

#loglambda <- predict(fit_simple, pixels(mesh, mask = COAST_simp, nx = 300, ny = 300), ~ mySmooth + Intercept, n.samples = 1000)

# Make a shortcut to a nicer colour scale:

# pl1 <- ggplot() + 
#   gg(lambda) + 
#   gg(COAST_simp) + 
#   ggtitle("LGCP fit to Points", subtitle = "(Response Scale)") + 
#   coord_fixed() +
#   colsc(lambda$median)
# 
# pl2 <- ggplot() + 
#   gg(loglambda) + 
#   gg(COAST_simp) + 
#   ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
#   coord_fixed() +
#   colsc(loglambda$median)
# 
# multiplot(pl1, pl2, cols = 2)
# 
# int.plot <- plot(fit_simple, "Intercept")
# spde.range <- spde.posterior(fit_simple, "mySmooth", what = "range")
# spde.logvar <- spde.posterior(fit_simple, "mySmooth", what = "log.variance")
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.logvar)
# 
# multiplot(range.plot, var.plot, int.plot)
