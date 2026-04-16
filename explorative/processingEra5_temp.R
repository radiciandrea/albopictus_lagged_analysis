setwd("/media/ddare/Elements/europe_tas/")
library(terra)
library(sf)
library(dplyr)
library(lubridate)
# library(tmap) #"World" shapefile

#get argentina polygon to crop the temp data
# data("World")
# sf::st_crs(World) #EPSG:4326
# arg<-World[World$name %in% "Argentina", ]
# plot(arg)

#load temperatures ncdf files
b<-list.files(recursive=TRUE, full.names = T, pattern = "tas_era5Land")

weeklyT_list <- lapply(b, rast)
weeklyT_r<- do.call(c, weeklyT_list)

# Calculate weekly average temperature for each year
temp_weekly <- terra::tapp(weeklyT_r, index = "yearweek", fun = "median", cores = 5)
plot(temp_weekly$yw_201801)

# Convert temperature from Kelvin to Celsius
temp_rast <- terra::app(temp_weekly, fun = function(x) { x - 273.15 }, cores = 5)
# terra::time(temp_rast) <- temp_dates
plot(temp_rast$yw_202325)

# Save cropped temperature data as a new NetCDF file
temp_output <- "tas_week_2018_2023.nc"
writeCDF(temp_rast, temp_output, overwrite = TRUE, 
         varname = "t2m", 
         longname = "2 metre median temperature", 
         unit = "°C")


#----- old code
weeklyT_list <- list()
for(i in 1:length(b)){
  # i=1
  tmp<-terra::rast(b[i])
  message("Processing month ", i, " of ", length(b), ": ", unique(lubridate::month(terra::time(tmp))), "/", unique(lubridate::year(terra::time(tmp))))
  cheatTab <-data.frame(date=as.Date(terra::time(tmp)))
  # cheatTab <- transform(cheatTab, id=as.numeric(lubridate::week(date)))
  cheatTab <- transform(cheatTab, id=as.numeric(factor(date)))
  tmp <- terra::tapp(tmp, index=cheatTab$id, fun = median, cores=4)
  # tmp <- terra::mask(tmp, arg)
  tmp <- terra::app(tmp, fun = function(x){x-273.15}, cores=4)
  terra::time(tmp)<-unique(cheatTab$date)
  weeklyT_list[i]<-tmp
  rm(tmp, cheatTab)
}
weeklyT_list<-do.call(c, weeklyT_list)

cheatTab <-data.frame(date=as.Date(terra::time(weeklyT_list)))
cheatTab <- transform(cheatTab, id=as.numeric(lubridate::week(date)))
cheatTab$id <- paste0(lubridate::year(cheatTab$date), "_", cheatTab$id)
weeklyT_list <- terra::tapp(weeklyT_list, index=cheatTab$id, fun = median, cores=4)
names(weeklyT_list) <- paste0("w_", unique(cheatTab$id))

plot(weeklyT_list$w_2023_53)
outname<-paste0("eu_weeklymedianT_2008_2023", ".RDS" )
# outname<-paste0("eu_weeklymedianT_2023", ".RDS" )
saveRDS(weeklyT_list, outname)

