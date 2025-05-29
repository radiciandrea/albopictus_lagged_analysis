setwd("/media/ddare/Elements/europe_tas/")
library(terra)
library(sf)
library(dplyr)
library(lubridate)

b<-list.files( recursive=TRUE, full.names = T, pattern = "prec_era5Land")

# Note: ERA5-Land "hourly total precipitation" refers to the accumulated total precipitation per hour.
# To get the daily cumulative precipitation for day i, use only the total at midnight (hour == 0).
# See more information here:
# - https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview
# - https://forum.ecmwf.int/t/total-precipitation-values-not-in-range-era5-land-hourly-data-from-1981-to-present/1128/7

weeklyP_list <- list()
for(i in 1:length(b)){ #
  # i=1
  # Load precipitation data from the NetCDF file
  cat(b[i])
  precip_rast <- terra::rast(b[i])
  
  # Extract only the layers corresponding to midnight observations (hour == 0)
  time_info <- data.frame(id = 1:terra::nlyr(precip_rast),
                          date = terra::time(precip_rast),
                          hour = lubridate::hour(terra::time(precip_rast)))
  
  midnight_info <- subset(time_info, hour == 0)
  precip_rast <- precip_rast[[midnight_info$id]]

    # Convert precipitation from meters to millimeters
  precip_rast <- terra::app(precip_rast, fun = function(x) { x * 1000 }, cores = 5)
  terra::time(precip_rast) <- midnight_info$date
  terra::units(precip_rast) <- "mm"
  weeklyP_list[[i]]<-precip_rast
}

weeklyP_r<-do.call(c, weeklyP_list)

# Calculate weekly total precipitation for each year
precip_weekly <- terra::tapp(weeklyP_r, index = "yearweek", fun = "sum", cores = 5)

plot(precip_weekly$yw_201801)
plot(precip_weekly$yw_202301)

# Save the processed precipitation data as a new NetCDF file
precip_output <- "tp_week_2018_2023.nc"
writeCDF(precip_weekly, precip_output, overwrite = TRUE, 
         varname = "tp", 
         longname = "Total precipitation", 
         unit = "mm")



#----- old code
weeklyP_list <- list()
for(i in 1:length(b)){ #
  # i=1
  tmp<-terra::rast(b[i])
  message("Processing month ", i, " of ", length(b), ": ", unique(lubridate::month(terra::time(tmp))), "/", unique(lubridate::year(terra::time(tmp))))
  cheatTab <-data.frame(date=as.Date(terra::time(tmp)))
  # cheatTab <- transform(cheatTab, id=as.numeric(lubridate::week(date)))
  #convert from meters to mm
  tmp<-tmp*1000
  terra::units(tmp)<-"mm"
  #get daily sum
  cheatTab <- transform(cheatTab, id=as.numeric(factor(date)))
  tmp <- terra::tapp(tmp, index=cheatTab$id, fun = sum, cores=4)
  # tmp <- terra::mask(tmp, arg)
  # tmp <- terra::app(tmp, fun = function(x){x-273.15}, cores=4)
  terra::time(tmp)<-unique(cheatTab$date)
  weeklyP_list[[i]]<-tmp
  rm(tmp, cheatTab)
}
weeklyP_list<-do.call(c, weeklyP_list)





cheatTab <- data.frame(date=as.Date(terra::time(weeklyP_list)))
cheatTab <- transform(cheatTab, id=sprintf("%02d", as.numeric(lubridate::week(cheatTab$date))))
cheatTab$id <- paste0(lubridate::year(cheatTab$date), "_", cheatTab$id)
pippo <- terra::tapp(weeklyP_list, index=cheatTab$id, fun = sum, cores=4)
names(pippo) <- paste0("w_", unique(cheatTab$id))

plot(pippo$w_2023_53)
outname<-paste0("eu_weeklymedianPrec_2008_2023", ".RDS" )
# outname<-paste0("eu_weeklymedianPrec_2023", ".RDS" )
saveRDS(pippo, outname)


myP<-readRDS("eu_weeklymedianPrec_2008_2022.RDS")
plot(myP$w_2008_1)
