setwd("/media/dared/Elements/europe_tas/")

library(terra)
library(dplyr)
library(tidyr)
library(lubridate)

#load tas and tp
tas <- terra::rast("tas_week_2008_2023.nc")
tp <- terra::rast("tp_week_2008_2023.nc")

#create dates
# Create a sequence of weekly dates
dates <- seq(ymd("2008-01-01"), ymd("2023-12-31"), by = "week")


#extract tas
tas.db <- as.data.frame(tas, xy=TRUE)
names(tas.db) <- c("x", "y", as.Date(dates))
tas.db <- tas.db %>% 
  pivot_longer(-c(x, y)) %>% 
  rename(tas=value, date=name)%>%
  mutate(date = as.Date(as.integer(date), origin = "1970-01-01"), 
         year = year(date),
         month=  month(date), 
         week = week(date)
         )

#extract tp
tp.db <- as.data.frame(tas, xy=TRUE)
names(tp.db) <- c("x", "y", as.Date(dates))
tp.db <- tp.db %>% 
  pivot_longer(-c(x, y)) %>% 
  rename(tp=value)


# join the databases
tas.db$tp <- tp.db$tp
head(tas.db)
saveRDS(tas.db, "envMatrix_2008_2023.RDS")

# reduce the matrix
m <- readRDS("envMatrix_2008_2023.RDS")
m <- m %>%
  filter (year >=2018)
saveRDS(m, "envMatrix_2018_2023.RDS")
