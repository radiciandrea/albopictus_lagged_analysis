# INTRO ----

#Code to fit and plot the lagged effects of predictors of the best crossvalidated model

# Libraries ----

library(tidyverse)
library(lubridate)
library(skimr)
library(Epi) 
library(dlnm) 
library(tsModel)
library(splines) 
library(mgcv)
library(gnm) 
library(pbs) 
library(plyr)
library(ggplot2)
library(pracma)
library(ggpubr)
library(scales)

# Load data frame ----

bio.matrix_sel <- readRDS("VectAbundance015_ER5.RDS") %>% 
  select(ID, Country, Region, year, latitude, longitude, week, date, value, 
         tas, tp) %>% 
  dplyr::rename(eggs=value, medianTweek= tas, cumPrecweek=tp, long= longitude, lat = latitude) %>% 
  dplyr::mutate(month=month(date), 
                eggs=floor(eggs)) %>% 
  as_tibble()

### add nweek and weekn

bio.matrix_sel <- bio.matrix_sel%>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(nweek = 1:n()) %>%
  dplyr::mutate(weekn = as.numeric(week)) %>%
  dplyr::ungroup()

# Splines and crossbasis definition ----

# tested models

# T.seas, T.P.seas, T.P.seas.trend (anche T.seas.trend, se servisse)

# natural spline over seasonality (df = 4)
bseas <- ns(bio.matrix_sel$weekn, df=4)

# df trend + natural spline
dftrend <- 4 * length(unique(bio.matrix_sel$year))
btrend <- ns(bio.matrix_sel$nweek, knots=equalknots(bio.matrix_sel$nweek, dftrend-1))

# crossbasis Temperature
tknots<-unname(quantile(bio.matrix_sel$medianTweek, probs = c(0.25, 0.75), na.rm = TRUE))
tlag = 12
cbTemp <- crossbasis(bio.matrix_sel$medianTweek, lag=tlag,
                     argvar=list(fun="bs", degree=2, knots=tknots),
                     arglag=list(knots=c(1,4)),
                     group=bio.matrix_sel$Region)

# crossbasis Precipitation
pknots<-unname(quantile(bio.matrix_sel$cumPrecweek, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
cbPrec <- crossbasis(bio.matrix_sel$cumPrecweek, lag=tlag,
                     argvar=list(fun="bs", degree=2, knots=pknots),
                     arglag=list(knots=c(1,4)),
                     group=bio.matrix_sel$Region)

# 1. T.seas ----

# T.seas and year 

# model definition 
mod_T.seas <- gnm(eggs ~ cbTemp + bseas + year,
                    eliminate=as.factor(ID), #grouping factor
                    data= bio.matrix_sel,
                    family=quasipoisson(link="log"),
                    na.action= "na.exclude")

## Effects ----

## Overall (Temperature) ----
#centering point at 15°C
myCen <- 15
cptmean <- crosspred(cbTemp, mod_T.seas, cen=myCen, by=1, at=6:30)

plot(cptmean, ptype="overall", xlab="Temperature (°C)", ylab="Increment ratio", col=2,
     main="Temperature: overall cumulative exposure-response")

## Marginal effects of lags at different temperatures ----
layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, var=6, xlab="Lag", ylab="Increment (ratio)", col=2, 
     main="Temperature = 6 °C", ylim = c(0, 3))

plot(cptmean, var=14, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 14 °C", ylim = c(0, 3))

plot(cptmean, var=22, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 22 °C", ylim = c(0, 3))

plot(cptmean, var=30, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 30 °C", ylim = c(0, 3))

layout(1)

## Marginal effects of temperatues at different lags ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, lag=0, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 0", ylim = c(0, 3))

plot(cptmean, lag=2, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 2", ylim = c(0, 3))

plot(cptmean, lag=4, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 4", ylim = c(0, 3))

plot(cptmean, lag=8, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 8", ylim = c(0, 3))

layout(1)

## Es. plot ----

#da replicare anche negli altri chunk

ID_es = 8649
year_es = 2022

eggs_mod_id = mod_T.seas$model %>%
  mutate(id_x = (year == year_es)*(`as.factor(ID)` == ID_es))

eggs_mod = mod_T.P.seas$fitted.values[which(eggs_mod_id$id_x ==1)]

eggs_obs = mod_T.P.seas$model %>%
  filter(eggs_mod_id$id_x ==1) %>%
  pull(eggs)

date_obs = bio.matrix_sel %>%
  filter(year == year_es) %>%
  filter(ID == ID_es) %>%
  filter(!is.na(eggs)) %>%
  pull(date)

plot(date_obs, eggs_obs)
points(date_obs, eggs_mod, col = 'blue')


# 2. T.P.seas ----

# T.P.seas and year 

# model definition 
mod_T.P.seas <- gnm(eggs ~ cbTemp + cbPrec + bseas + year,
                  eliminate=as.factor(ID), #grouping factor
                  data= bio.matrix_sel,
                  family=quasipoisson(link="log"),
                  na.action= "na.exclude")

## Effects ----

## Overall (Temperature) ----
#centering point at 15°C
myCen <- 15
cptmean <- crosspred(cbTemp, mod_T.P.seas, cen=myCen, by=1, at=6:30)

plot(cptmean, ptype="overall", xlab="Temperature (°C)", ylab="Increment ratio", col=2,
     main="Temperature: overall cumulative exposure-response")


## Marginal effects of lags at different temperatures ----
layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, var=6, xlab="Lag", ylab="Increment (ratio)", col=2, 
     main="Temperature = 6 °C", ylim = c(0, 3))

plot(cptmean, var=14, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 14 °C", ylim = c(0, 3))

plot(cptmean, var=22, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 22 °C", ylim = c(0, 3))

plot(cptmean, var=30, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 30 °C", ylim = c(0, 3))

layout(1)

## Marginal effects of temperature at different lags ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, lag=0, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 0", ylim = c(0, 3))

plot(cptmean, lag=2, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 2", ylim = c(0, 3))

plot(cptmean, lag=4, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 4", ylim = c(0, 3))

plot(cptmean, lag=8, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 8", ylim = c(0, 3))

layout(1)

## Overall (Precipitation) ----
#centering point at 20 mm
myCen <- 3
cppmean <- crosspred(cbPrec, mod_T.P.seas, cen=myCen, by=1, at=0:40)

plot(cppmean, ptype="overall", xlab="Precipitation", ylab="Increment ratio", col=2,
     main="Precipitation: overall cumulative exposure-response", ylim = c(0.5,3))

## Marginal effects of lags at different cumulative precipitation ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cppmean, var=0, xlab="Lag", ylab="Increment (ratio)", col=2, 
     main="Precipitation = 0 mm", ylim = c(0.75, 1.25))

plot(cppmean, var=10, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 10 mm", ylim = c(0.75, 1.25))

plot(cppmean, var=20, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 20 mm", ylim = c(0.75, 1.25))

plot(cppmean, var=30, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 30 mm", ylim = c(0.75, 1.25))

layout(1)

## Marginal effects of precipitation at different lags ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cppmean, lag=0, xlab="P (mm)", ylab="Increment (ratio)", col=2, 
     main="lag = 0", ylim = c(0,2))

plot(cppmean, lag=2, xlab="P (mm)", ylab="Increment (ratio)", col=2, 
     main="lag = 2", ylim = c(0,2))

plot(cppmean, lag=4, xlab="P (mm)", ylab="Increment (ratio)", col=2,
     main="lag = 4", ylim = c(0,2))

plot(cppmean, lag=8, xlab="P (mm)", ylab="Increment (ratio)", col=2,
     main="lag = 8", ylim = c(0,2))

layout(1)

# 3. T.P.seas.trend ----

# model definition 
mod_T.P.seas.trend <- gnm(eggs ~ cbTemp + cbPrec + btrend,
                    eliminate=as.factor(ID), #grouping factor
                    data= bio.matrix_sel,
                    family=quasipoisson(link="log"),
                    na.action= "na.exclude")

## Effects ----

## Overall (Temperature) ----
#centering point at 15°C
myCen <- 15
cptmean <- crosspred(cbTemp, mod_T.P.seas.trend, cen=myCen, by=1, at=6:30)

plot(cptmean, ptype="overall", xlab="Temperature (°C)", ylab="Increment ratio", col=2,
     main="Temperature: overall cumulative exposure-response")

## Marginal effects of lags at different temperatures ----
layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, var=6, xlab="Lag", ylab="Increment (ratio)", col=2, 
     main="Temperature = 6 °C", ylim = c(0, 3))

plot(cptmean, var=14, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 14 °C", ylim = c(0, 3))

plot(cptmean, var=22, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 22 °C", ylim = c(0, 3))

plot(cptmean, var=30, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Temperature = 30 °C", ylim = c(0, 3))

layout(1)

## Marginal effects of temperature at different lags ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, lag=0, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 0", ylim = c(0, 3))

plot(cptmean, lag=2, xlab="T (°C)", ylab="Increment (ratio)", col=2, 
     main="lag = 2", ylim = c(0, 3))

plot(cptmean, lag=4, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 4", ylim = c(0, 3))

plot(cptmean, lag=8, xlab="T (°C)", ylab="Increment (ratio)", col=2,
     main="lag = 8", ylim = c(0, 3))

layout(1)

#centering point at 20 mm

myCen <- 20
cppmean <- crosspred(cbPrec, mod_T.P.seas.trend, cen=myCen, by=1, at=0:275)

## Overall (Precipitation) ----
plot(cppmean, ptype="overall", xlab="Precipitation", ylab="Increment ratio", col=2,
     main="Precipitation: overall cumulative exposure-response", ylim = c(0,5))

## Marginal effects of lags at different cumulative precipitation ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cppmean, var=0, xlab="Lag", ylab="Increment (ratio)", col=2, 
     main="Precipitation = 0 mm", ylim = c(0,2))

plot(cppmean, var=50, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 50 mm", ylim = c(0,2))

plot(cppmean, var=100, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 100 mm", ylim = c(0,2))

plot(cppmean, var=200, xlab="Lag", ylab="Increment (ratio)", col=2,
     main="Precipitation = 200 mm", ylim = c(0,2))

layout(1)

## Marginal effects of precipitation at different lags ----

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cppmean, lag=0, xlab="P (mm)", ylab="Increment (ratio)", col=2, 
     main="lag = 0", ylim = c(0,2))

plot(cppmean, lag=2, xlab="P (mm)", ylab="Increment (ratio)", col=2, 
     main="lag = 2", ylim = c(0,2))

plot(cppmean, lag=4, xlab="P (mm)", ylab="Increment (ratio)", col=2,
     main="lag = 4", ylim = c(0,2))

plot(cppmean, lag=8, xlab="P (mm)", ylab="Increment (ratio)", col=2,
     main="lag = 8", ylim = c(0,2))

layout(1)
