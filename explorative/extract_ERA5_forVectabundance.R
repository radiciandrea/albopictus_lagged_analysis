setwd("/media/dared/Elements/europe_tas/")

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)
library(mapview)

#load vectabundance
vabb <- read_xlsx("Vectabundace_v015.xlsx")
names(vabb)

#creat join ID
vabb$joinID <- paste0(vabb$ID, "_", vabb$year, "_", vabb$week)

#get spatial location
vabb.sf <- vabb %>% 
  select(ID, longitude, latitude) %>%
  distinct() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs=4326)

#load tas and tp
tas <- terra::rast("tas_week_2018_2023.nc")
tp <- terra::rast("tp_week_2018_2023.nc")

#interpolate to avoid NAs located in the sea
w <- matrix(1, 3, 3)
tas.i <- terra::focal(tas, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)
tp.i <- terra::focal(tp, w, mean, na.rm=TRUE, NAonly=TRUE, pad=TRUE)

mapview(vabb.sf) + mapview()

#create dates
# Create a sequence of weekly dates
dates <- seq(ymd("2008-01-01"), ymd("2023-12-31"), by = "week")

# Extract year and ISO week
year_week <- dates %>%
  tibble(date = .) %>%
  mutate(
    year = isoyear(date),
    week = isoweek(date),
    year_week = sprintf("%d_%02d", year, week)
  ) %>%
  pull(year_week)

#extract tas
tas.db <- terra::extract(tas.i, vabb.sf, ID=TRUE)
tas.db$ID<-vabb.sf$ID
names(tas.db) <- c("ID", year_week)
tas.db <- tas.db %>% 
  pivot_longer(-c(ID)) %>% 
  mutate(joinID = paste0(ID, "_", name)) %>% 
  rename(tas=value)

#extract tp
tp.db <- terra::extract(tp.i, vabb.sf, ID=TRUE)
tp.db$ID<-vabb.sf$ID
names(tp.db) <- c("ID", year_week)
tp.db <- tp.db %>% 
  pivot_longer(-c(ID)) %>% 
  mutate(joinID = paste0(ID, "_", name)) %>% 
  rename(tp=value)
dim(tas.db); dim(tp.db)

# join the databases
tas.db$tp <- tp.db$tp
head(tas.db)

vabb <- vabb %>% 
  left_join(tas.db %>%  select(joinID, tas, tp), by="joinID")

saveRDS(vabb, "VectAbundance015_ER5.RDS")
