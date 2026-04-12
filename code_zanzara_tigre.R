
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
library(terra)
library(sf)
library(mapview)


bio.matrix_sel <- readRDS("VectAbundance015_ER5.RDS") %>% 
  select(ID, Country, Region, year, latitude, longitude, week, date, value, 
         tas, tp) %>% 
  rename(eggs=value, medianTweek= tas, cumPrecweek=tp, long= longitude, lat = latitude) %>% 
  mutate(month=month(date), 
         eggs=floor(eggs)) %>% 
  as_tibble()
bio.matrix_sel

#library(MASS)

plot(bio.matrix_sel$date)
length(unique(bio.matrix_sel$ID))

bio.matrix_sel %>% 
  select(long, lat) %>% 
  distinct() %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  mapview()

# other methods should be used for validation: e.g. leave-one-out or cross-validation

ID_val = "6776"

mydf <- bio.matrix_sel %>% 
  filter(ID != ID_val) 

# plot validation plot location
bio.matrix_sel %>% 
  filter(ID == ID_val) %>%  
  select(long, lat) %>% 
  distinct() %>% 
  st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
  mapview()

plot(eggs~date,  data=mydf, cex=0.05 ,  xlab="Date", ylab="Number of eggs", bty="l")

summary(mydf$eggs)
# I Dati mancanti "NA" sono diversi dagli "0"?

# summary by month
tapply(mydf$eggs, mydf$month,
       function(x) format(summary(x), scientific = TRUE))

# I mesi informativi sembrano essere da maggio ad ottobre: meglio un'analisi stagionale?

par(mfrow=c(1,2))
hist(mydf$eggs, main = "Histogram of Number of eggs")

ID_subset = "9347"

dsub <- subset(mydf, ID==ID_subset) 

par(mar=c(2,3,1,1))
layout(matrix(1:3, nrow=3), heights=c(2,2,2))
plot(eggs~date, data=dsub, type="s", lty=2, bty="l", xlab="",
     ylab="Number of eggs", mgp=c(2.2,0.7,0), lab=c(5,3,7))
plot(medianTweek~date, data=dsub, type="h", lty=1, bty="l", col=2, xlab="",
     ylab="Temperature  (°C)")
plot(cumPrecweek~date, data=dsub, type="h", lty=1, bty="l", col=4, xlab="Date",
     ylab="Cumulative precipitation (mm)")
hist(log10(mydf$eggs + 1 ), main = "Histogram of log10(Number of eggs +1)")

layout(1)
par(mar=c(5,4,4,1)+0.1)

# trasformazione logaritmica non mi convince molto, meglio quasipoisson o binomiale negativa?

mydf$eggsL10<- log10(mydf$eggs+1)

plot(unique(mydf$lat), ylim=c(40,50), xaxt='n', xlab="", ylab="Latitude",  lty=1, bty="n")

sum(is.na(bio.matrix_sel$date))

mydf$weekn<-as.numeric(mydf$week)
mydf <-mydf %>% group_by(ID) %>%  mutate(nweek = 1:n())

# long interannual temporal trend
# un pò strano il modo di definire il numero di knots, sceglierei un numero di knots per anno (e.g. 4) moltiplicati per il numeri degli anni 5
#dftrend <- round(as.numeric(diff(range(mydf$date))/365.25 * length(unique(mydf$year)))) 
dftrend <- 4 * length(unique(mydf$year))

# nei mesi invernali ci sono molti 0 sarebbe meglio dare più flessibilità ai mesi estivi
# visto che sono dati settimanale forse conviene lavorare direttamente con la variabile week
#btrend <- ns(mydf$date, knots=equalknots(mydf$date, dftrend-1))
btrend <- ns(mydf$nweek, knots=equalknots(mydf$nweek, dftrend-1))


# seasonal intraannual trend
#bseas <- ns(mydf$month, df=4) #week
bseas <- ns(mydf$weekn, df=4) #week
# i comandi precedenti definiscono una struttura che cattura sia il trend che la stagionalità, 
# perchè viene introdotto un termine anche per la stagionalità?

tknots<-unname(quantile(mydf$medianTweek, probs = c(0.25, 0.75), na.rm = TRUE))
#original

cbTemp <- crossbasis(mydf$medianTweek, lag=8,
                     argvar=list(fun="bs", degree=2, knots=tknots),
                     arglag=list(knots=c(1,4)),
                     group=mydf$Region)

# Fit mod_Tel to analyse coefficients
# mod_T <- gnm(eggsL10 ~ cbTemp,
#              eliminate=as.factor(ID), #grouping factor
#              data=mydf,
#              family=gaussian(),
#              na.action= "na.exclude")

# modello quasipoisson
mod_T <- gnm(eggs ~ cbTemp,
             eliminate=as.factor(ID), #grouping factor
             data=mydf,
             family=quasipoisson(link="log"),
             na.action= "na.exclude")

summary(mod_T)

# un esempio con binomiale negativa

library(MASS)

# mod_T <- glm.nb(eggs ~ cbTemp,
#              #eliminate=as.factor(ID), #grouping factor
#              data=mydf,
#              na.action= "na.exclude")

# modello con stagionalità + trend, si potrebbe anche inserire un termine d'interazione con l'anno come variabile categorica

# mod_T.seas <- gnm(eggsL10 ~ cbTemp+ bseas+year, # + btrend,
#                   eliminate=as.factor(ID), #grouping factor
#                   data=mydf,
#                   family=gaussian(),
#                   na.action= "na.exclude")

mod_T.seas <- gnm(eggs ~ cbTemp+ bseas+year, # + btrend,
                  eliminate=as.factor(ID), #grouping factor
                  data=mydf,
                  family=quasipoisson(link="log"),
                  na.action= "na.exclude")

summary(mod_T.seas)

# init<-mod_T.seas$coefficients
# 
# # un esempio con binomiale cbTemp+ bseas+year
# mod_T.seas <- glm.nb(eggs ~ cbTemp+ bseas+year,
#                 #eliminate=stratum, #grouping factor
#                 data=mydf,
#                 na.action= "na.exclude")

#modello con spline trend & seasonality

# mod_T.seas.trend <- gnm(eggsL10 ~ cbTemp+ + btrend,
#                   eliminate=as.factor(ID), #grouping factor
#                   data=mydf,
#                   family=gaussian(),
#                   na.action= "na.exclude")

mod_T.seas.trend <- gnm(eggs ~ cbTemp+ + btrend,
                        eliminate=as.factor(ID), #grouping factor
                        data=mydf,
                        family=quasipoisson(link="log"),
                        na.action= "na.exclude")

summary(mod_T.seas.trend)

# define predictive model
cbTemp_df <- as.data.frame(cbTemp)
colnames(cbTemp_df) <- paste0("cb", seq_len(ncol(cbTemp_df)))

mydf_ext <- cbind(mydf, cbTemp_df)

# Build formula with those columns
#form <- as.formula(paste("eggsL10 ~", paste(colnames(cbTemp_df), collapse = " + ")))
form <- as.formula(paste("eggs ~", paste(colnames(cbTemp_df), collapse = " + ")))

# Fit mod_Tel to predict on new data
# mod_T_pred <- gnm(form,
#                   data = mydf_ext,
#                   family = gaussian(),
#                   na.action = "na.exclude")

mod_T_pred <- gnm(form,
                  eliminate=as.factor(ID),
                  data = mydf_ext,
                  family = quasipoisson(link="log"),
                  na.action = "na.exclude")

# mod_T_pred <- glm.nb(form,
#                   data = mydf_ext,
#                   na.action = "na.exclude")


bseas_df <- as.data.frame(bseas)
colnames(bseas_df) <- paste0("bt", seq_len(ncol(bseas_df)))
mydf_ext.seas <- cbind(mydf, cbTemp_df, bseas_df) #btrend

# Build formula with those columns
form <- as.formula(paste("eggs ~", paste(colnames( cbind(cbTemp_df,  bseas_df,year=mydf$year)), collapse = " + "))) # btrend_df,

# mod_T.seas_pred <- gnm(form,
#                        data = mydf_ext.seas,
#                        family = gaussian(),
#                        na.action = "na.exclude")

mod_T.seas_pred <- gnm(form,
                       eliminate=as.factor(ID),
                       data = mydf_ext.seas,
                       family = quasipoisson(link="log"),
                       na.action = "na.exclude")

#centering point at 15°C
myCen <- 15
cptmean <- crosspred(cbTemp, mod_T, cen=myCen, by=1, at=0:30)
cptmean.seas <- crosspred(cbTemp, mod_T.seas, cen=myCen, by=1, at=0:30)


# Overall cumulative
plot(cptmean, ptype="overall", xlab="Temperature (°C)", ylab="Increment ratio", col=2,
     main="Temperature: overall cumulative exposure-response")

plot(cptmean.seas, ptype="overall", xlab="Temperature (°C)", ylab="Increment ratio", col=2,
     main="Temperature + Seasonality: overall cumulative exposure-response")

# the data have a nested structure they should be analysed with a two-stage design

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, var=5, xlab="Lag", ylab="Increment (Log10)", col=2, 
     main="Temperature at 5 °C")

plot(cptmean, var=20, xlab="Lag", ylab="Increment (Log10)", col=2,
     main="Temperature at 20 °C")

plot(cptmean, var=25, xlab="Lag", ylab="Increment (Log10)", col=2,
     main="Temperature at 25 °C")

plot(cptmean, var=30, xlab="Lag", ylab="Increment (Log10)", col=2,
     main="Temperature at 30 °C")

layout(1)

layout(matrix(1:4, ncol=2, nrow = 2))
plot(cptmean, lag=0, xlab="T", ylab="Increment (Log10)", col=2, 
     main="lag=0")

plot(cptmean, lag=2, xlab="T", ylab="Increment (Log10)", col=2, 
     main="lag=2")

plot(cptmean, lag=4, xlab="T", ylab="Increment (Log10)", col=2,
     main="lag=4")

plot(cptmean, lag=8, xlab="T", ylab="Increment (Log10)", col=2,
     main="lag=8")

# fitted values

mydf$fitted <- fitted(mod_T_pred)
mydf$fitted.seas <- fitted(mod_T.seas_pred)

mydf_temp <- mydf %>%
  filter(ID == ID_subset)

# temperature only
plot(mydf_temp$date, mydf_temp$eggs, type="p" ,col="blue", ylab="Y", main= "Temperature only model")
lines(mydf_temp$date, mydf_temp$fitted, col="darkgreen", lty=2)
legend("topright", legend=c("Observed", "Fitted T"), col=c("blue", "darkgreen"), lty=c(1,2))
text(x=mydf_temp$date[100],y=500, labels=paste("Pearson's r:", round(cor.test(mydf_temp$eggs, mydf_temp$fitted)$estimate,3)))

# with season and trend
plot(mydf_temp$date, mydf_temp$eggs, type="p" ,col="blue", ylab="Y", main= "Temperature + seas model")
lines(mydf_temp$date, mydf_temp$fitted.seas, col="darkgreen", lty=2)
legend("topright", legend=c("Observed", "Fitted T + Seas"), col=c("blue", "darkgreen"), lty=c(1,2))
text(x=mydf_temp$date[100],y=500, labels=paste("Pearson's r:", round(cor.test(mydf_temp$eggs, mydf_temp$fitted.seas)$estimate,3)))

# plot in the validation sample
mydf_val <- bio.matrix_sel %>% 
  filter(ID == ID_val)

regions <- unique(mydf$Region)

# Build newdf with same shape as original (if needed)
newdf <- data.frame(date = mydf_val$date,
  medianTweek = mydf_val$medianTweek, ID=rep(3339,nrow(mydf_val)))

# Build new crossbasis
cbTemp_val <- crossbasis(newdf$medianTweek, lag=8,
                         argvar=list(fun="bs", degree=2, knots=tknots),
                         arglag=list(knots=c(1,4)),
                         group=newdf$Region)

cbTemp_val_df <- as.data.frame(cbTemp_val)
colnames(cbTemp_val_df) <- paste0("cb", seq_len(ncol(cbTemp_val_df)))
newdf_ext <- cbind(newdf, cbTemp_val_df)

# Predict
log10eggsPred <- predict(mod_T_pred, newdata = newdf_ext, type = "response")
mydf_val$eggsPred  = (log10eggsPred)

#plot
plot(mydf_val$date, mydf_val$eggs, col = "blue")
lines(newdf$date, mydf_val$eggsPred , col = "darkgreen", lty=2)
# plot(newdf$date, predY, col = "darkgreen")
legend("topright", legend=c("observed eggs", "predicted eggs"), col=c("blue", "darkgreen"), lty=c(1,2))
text(x=mydf_val$date[100],y=500, labels=paste("Pearson's r:", round(cor.test(mydf_val$eggs, mydf_val$eggsPred)$estimate,3)))

# long interannual temporal trend
dftrend <- 4 * length(unique(mydf_val$year))
mydf_val <-mydf_val %>% group_by(ID) %>%  mutate(nweek = 1:n())
btrend <- ns(mydf_val$nweek, knots=equalknots(mydf_val$nweek, dftrend-1))

# seasonal intraannual trend
mydf_val$weekn<-as.numeric(mydf_val$week)
bseas_pred <- ns(mydf_val$weekn, df=4) #week

bseas_pred_df <- as.data.frame(bseas_pred)
colnames(bseas_pred_df) <- paste0("bt", seq_len(ncol(bseas_pred_df)))

mydf_ext.seas <- cbind(newdf, cbTemp_val_df, bseas_pred_df, year=mydf_val$year)

# Predict
log10eggsPred <- predict(mod_T.seas_pred, newdata = mydf_ext.seas, type = "response")
mydf_val$eggsPred_seas  = log10eggsPred
mydf_val$eggsPred_seas[is.na(mydf_val$eggs)]<-NA

# il confronto fra valori predetti e ossrvati andrebbe fatto sui dati non mancanti

# plot
plot(mydf_val$date, mydf_val$eggs, type="p" ,col="blue", ylab="Y", main= "Temperature + seas model")
lines(mydf_val$date, mydf_val$eggsPred_seas , col="darkgreen", lty=2)
legend("topright", legend=c("Observed", "Fitted T + Seas"), col=c("blue", "darkgreen"), lty=c(1,2))
text(x=mydf_val$date[100],y=500, labels=paste("Pearson's r:", round(cor.test(mydf_val$eggs, mydf_val$eggsPred_seas )$estimate,3)))
