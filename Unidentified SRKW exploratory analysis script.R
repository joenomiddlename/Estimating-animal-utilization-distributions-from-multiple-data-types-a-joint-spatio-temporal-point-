# Script to identify unidentified SRKW pods.

WW$DateTimes = as.POSIXct(3600*24 * (WW$Matdate - min(WW$Matdate)), 
                     tz = 'UTC', origin = '2009-05-04 14:15', 
                     format = "%Y-%m-%d %H:%M")
WW$Date = as.Date(WW$DateTimes)
SRKW_days = unique(WW@data$Date[WW@data$Pod == 'SRKW'])

View(WW@data[WW@data$Date %in% SRKW_days,c('Date','DateTimes','Pod')])

hist(hour(WW@data$DateTimes))
hist(hour(WW@data$DateTimes[WW@data$Pod == 'SRKW']), add=T)
# It looks like SRKW is usually found earlier. Agrees with Hauser et al 2006.

table(by(WW@data[WW@data$Date %in% SRKW_days,], WW@data$Date[WW@data$Date %in% SRKW_days], FUN = function(x){
  length(x$Pod[x$Pod != 'SRKW'])
}, simplify = T)) # 78 never-identified sightings/pods

table(by(WW@data[WW@data$Date %in% SRKW_days,], WW@data$Date[WW@data$Date %in% SRKW_days], FUN = function(x){
  min(hour(x$DateTimes[length(x$Pod[x$Pod != 'SRKW']) == 0]))
}, simplify = T)) # some sightings made at strange times of day.

table(by(WW@data[WW@data$Date %in% SRKW_days,], WW@data$Date[WW@data$Date %in% SRKW_days], FUN = function(x){
  min(hour(x$DateTimes))
}, simplify = T)) # some sightings made at strange times of day - could these be BC Ferries?.

length(unique(WW@data$Date))

# return only the observations that were ultimately not verified
temp = by(WW@data[WW@data$Date %in% SRKW_days,], WW@data$Date[WW@data$Date %in% SRKW_days], FUN = function(x){
  if(length(x$Pod[x$Pod != 'SRKW']) == 0){
    return(x)
  }
  }, simplify = F)
temp2 = data.frame()
for(i in 1:length(temp))
{
  if(!is.null(temp[[i]]))
  {
    temp2 = rbind(temp2, temp[[i]])
  }
}
temp2 = SpatialPointsDataFrame(coords = cbind(temp2$lon, temp2$lat),
                               data = data.frame(month = temp2$Month,
                                                 Date = temp2$Date),
                               proj4string = COAST_simp@proj4string)
ggplot() + gg(COAST_plotting) + gg(temp2) + ggtitle('Locations of never-identified sightings')

never_identified_dates = unique(temp2@data$Date)
saveRDS(never_identified_dates,'never_identified_dates.rds')

temp2@data$month