# INTRO ----

#Code to run a model selection based by a cross-validation upon non linear lagged model to study the relation between weather and oviposition

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

# Define cross-valitation settings ----

n_ids = nrow(bio.matrix_sel)
reps = 20
prop_fit = 0.7

#settings lag
tlag = 8

# add n row + and new column for stratified sampling

bio.matrix_sel <- bio.matrix_sel%>%
  dplyr::mutate(IDdata = 1:n_ids) %>%
  dplyr::mutate(IDy = paste0(ID, "_", year))

# Performances df ----
#we define a data.frame with 10 rows, one for random repetitions, and n+1 columns (n = models)

# RMSE: root mean square error
RMSE_val_df = data.frame(rep = 1:reps,
                         seas = NA,
                         seas.trend = NA,
                         Temp = NA,
                         T.seas = NA,
                         T.P.seas = NA,
                         T.seas.trend = NA,
                         T.P.seas.trend = NA)

RMSE_fit_df = RMSE_val_df

# MAE: mean absolute error
MAE_val_df = RMSE_val_df
MAE_fit_df = RMSE_val_df

# RMSLE: root mean square logarithmic error
RMSLE_val_df = RMSE_val_df
RMSLE_fit_df = RMSE_val_df

# Residual deviance
ResDev_df = RMSE_val_df

# Overdispersion
Dispersion_df = RMSE_val_df

# Cycle to fit and validate the model ----
tic()
for(i in 1:reps){
  
  # for stratified sampling:
  
  # I must leave the same IDs and the same years in both subsets
  # I should not remove weather from any subsets (but only observations)

  # extract random batch from fit and validation  ----
  ids_i = bio.matrix_sel %>%
    dplyr::group_by(IDy)%>%
    sample_n(size = 52*prop_fit) %>% pull(IDdata)
    
  bio.matrix_fit_i = bio.matrix_sel
  bio.matrix_fit_i$eggs[-ids_i] = NA
  
  bio.matrix_val_i = bio.matrix_sel
  bio.matrix_val_i$eggs[ids_i] = NA
  
  # Extract obs_val for validation
  obs_val = bio.matrix_val_i$eggs
  obs_val_na = which(is.na(obs_val))
  
  # 1. Seas and year ----
  
  # natural spline over seasonality (df = 4)
  bseas <- ns(bio.matrix_fit_i$weekn, df=4)
  
  # predictive formulation 
  bseas_df <- as.data.frame(bseas)
  colnames(bseas_df) <- paste0("bs", seq_len(ncol(bseas_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_seas <- cbind(bio.matrix_fit_i, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(bseas_df), "year"), collapse = " + ")))
  
  # model definition 
  # with only spline trend & seasonality
  mod_seas_i <- gnm(form,
                  eliminate=as.factor(ID), #grouping factor
                  data=bio.matrix_fit_i_ext_seas,
                  family=quasipoisson(link="log"),
                  na.action= "na.exclude")
  
  # assess fitting error
  
  pred_seas_i = mod_seas_i$fitted
  obs_fit_i = mod_seas_i$fitted + mod_seas_i$residuals
  
  RMSE_fit_df$seas[i] = sqrt(mean((pred_seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$seas[i]  = sqrt(mean((log(pred_seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$seas[i]  = mean(abs(pred_seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # natural spline over trend 
  bseas_val <- ns(bio.matrix_val_i$weekn, df=4)
  
  bseas_val_df <- as.data.frame(bseas_val)
  colnames(bseas_val_df) <- paste0("bs", seq_len(ncol(bseas_val_df)))
  
  bio.matrix_val_i_seas <- cbind(bio.matrix_val_i, bseas_val_df)
  
  # predict
  pred_seas <- predict(mod_seas_i, newdata = bio.matrix_val_i_seas, type = "response")
  
  RMSE_val_df$seas[i] = sqrt(mean((pred_seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$seas[i]  = sqrt(mean((log(pred_seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$seas[i]  = mean(abs(pred_seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 2. Seas.trend ----
  
  # df trend
  dftrend <- 4 * length(unique(bio.matrix_fit_i$year))
  
  # natural spline over trend
  btrend <- ns(bio.matrix_fit_i$nweek, knots=equalknots(bio.matrix_fit_i$nweek, dftrend-1))
  
  # predictive formulation 
  btrend_df <- as.data.frame(btrend)
  colnames(btrend_df) <- paste0("bt", seq_len(ncol(btrend_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_seas.trend <- cbind(bio.matrix_fit_i, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(btrend_df), collapse = " + ")))
  
  # model definition 
  # with only spline trend & seasonality
  mod_seas.trend_i <- gnm(form,
                        eliminate=as.factor(ID), #grouping factor
                        data=bio.matrix_fit_i_ext_seas.trend,
                        family=quasipoisson(link="log"),
                        na.action= "na.exclude")
  
  # assess fitting error
  
  pred_seas.trend_i = mod_seas.trend_i$fitted
  obs_fit_i = mod_seas.trend_i$fitted + mod_seas.trend_i$residuals
  
  RMSE_fit_df$seas.trend[i] = sqrt(mean((pred_seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$seas.trend[i]  = sqrt(mean((log(pred_seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$seas.trend[i]  = mean(abs(pred_seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # natural spline over trend (same dftrend as before)
  btrend_val <- ns(bio.matrix_val_i$nweek, knots=equalknots(bio.matrix_val_i$nweek, dftrend-1))
  
  btrend_val_df <- as.data.frame(btrend_val)
  colnames(btrend_val_df) <- paste0("bt", seq_len(ncol(btrend_val_df)))
  
  bio.matrix_val_i_seas.trend <- cbind(bio.matrix_val_i, btrend_val_df)
  
  # predict
  pred_seas.trend <- predict(mod_seas.trend_i, newdata = bio.matrix_val_i_seas.trend, type = "response")
  
  RMSE_val_df$seas.trend[i] = sqrt(mean((pred_seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$seas.trend[i]  = sqrt(mean((log(pred_seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$seas.trend[i]  = mean(abs(pred_seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 3 Temp  ----
  
  # knots for temperature
  tknots<-unname(quantile(bio.matrix_fit_i$medianTweek, probs = c(0.25, 0.75), na.rm = TRUE))
  
  # cross basis over temperature
  cbTemp <- crossbasis(bio.matrix_fit_i$medianTweek, lag=tlag,
                       argvar=list(fun="bs", degree=2, knots=tknots),
                       arglag=list(knots=c(1,4)),
                       group=bio.matrix_fit_i$Region)
  
  # define model (preditive formulation)
  cbTemp_df <- as.data.frame(cbTemp)
  colnames(cbTemp_df) <- paste0("cbt", seq_len(ncol(cbTemp_df)))
  
  # and cbind
  bio.matrix_fit_i_ext_T <- cbind(bio.matrix_fit_i, cbTemp_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbTemp_df), collapse = " + ")))
  
  # model definition 
  # with temperature, spline trend & seasonality
  mod_T_i <- gnm(form,
               eliminate=as.factor(ID), #grouping factor
               data=bio.matrix_fit_i_ext_T,
               family=quasipoisson(link="log"),
               na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T_i = mod_T_i$fitted
  obs_fit_i = mod_T_i$fitted + mod_T_i$residuals
 
  RMSE_fit_df$Temp[i] = sqrt(mean((pred_T_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$Temp[i]  = sqrt(mean((log(pred_T_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$Temp[i]  = mean(abs(pred_T_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # Build new crossbasis
  cbTemp_val <- crossbasis(bio.matrix_val_i$medianTweek, lag=tlag,
                           argvar=list(fun="bs", degree=2, knots=tknots),
                           arglag=list(knots=c(1,4)),
                           group=bio.matrix_val_i$Region)
  
  cbTemp_val_df <- as.data.frame(cbTemp_val)
  colnames(cbTemp_val_df) <- paste0("cbt", seq_len(ncol(cbTemp_val_df)))
  
  bio.matrix_val_i_T <- cbind(bio.matrix_val_i, cbTemp_val_df)
  
  # predict
  pred_T<- predict(mod_T_i, newdata = bio.matrix_val_i_T, type = "response")
  
  RMSE_val_df$Temp[i] = sqrt(mean((pred_T[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$Temp[i]  = sqrt(mean((log(pred_T[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$Temp[i]  = mean(abs(pred_T[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 4 T.seas  and year ----
  
  # cbind
  bio.matrix_fit_i_ext_T.seas <- cbind(bio.matrix_fit_i, cbTemp_df, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(cbind(cbTemp_df, bseas_df)), "year"), collapse = " + ")))
  
  # model definition 
  # with temperature, spline seasonality
  mod_T.seas_i <- gnm(form,
                    eliminate=as.factor(ID), #grouping factor
                    data=bio.matrix_fit_i_ext_T.seas,
                    family=quasipoisson(link="log"),
                    na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T.seas_i = mod_T.seas_i$fitted
  obs_fit_i = mod_T.seas_i$fitted + mod_T.seas_i$residuals
  
  RMSE_fit_df$T.seas[i] = sqrt(mean((pred_T.seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.seas[i]  = sqrt(mean((log(pred_T.seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.seas[i]  = mean(abs(pred_T.seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.seas <- cbind(bio.matrix_val_i, cbTemp_val_df, bseas_val_df)
  
  # predict
  pred_T.seas<- predict(mod_T.seas_i, newdata = bio.matrix_val_i_T.seas, type = "response")
  
  RMSE_val_df$T.seas[i] = sqrt(mean((pred_T.seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.seas[i]  = sqrt(mean((log(pred_T.seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.seas[i]  = mean(abs(pred_T.seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 5. T.P.seas and year ----
  
  # knots for precipitation
  pknots<-unname(quantile(bio.matrix_fit_i$cumPrecweek, probs = c(0.25, 0.75, 0.9), na.rm = TRUE))
  
  cbPrec <- crossbasis(bio.matrix_fit_i$cumPrecweek, lag=tlag,
                       argvar=list(fun="bs", degree=2, knots=pknots),
                       arglag=list(knots=c(1,4)),
                       group=bio.matrix_fit_i$Region)
  
  # define model (preditive formulation)
  cbPrec_df <- as.data.frame(cbPrec)
  colnames(cbPrec_df) <- paste0("cbp", seq_len(ncol(cbPrec_df)))
  
  # cbind
  bio.matrix_fit_i_ext_T.P.seas <- cbind(bio.matrix_fit_i, cbTemp_df, cbPrec_df, bseas_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(c(colnames(cbind(cbTemp_df, cbPrec_df,bseas_df)), "year"), collapse = " + ")))
  
  # model definition 
  # with temperature, spline seasonality
  mod_T.P.seas_i <- gnm(form,
                      eliminate=as.factor(ID), #grouping factor
                      data=bio.matrix_fit_i_ext_T.P.seas,
                      family=quasipoisson(link="log"),
                      na.action= "na.exclude")
  
  # assess fitting error
  
  pred_T.P.seas_i = mod_T.P.seas_i$fitted
  obs_fit_i = mod_T.P.seas_i$fitted + mod_T.P.seas_i$residuals
  
  RMSE_fit_df$T.P.seas[i] = sqrt(mean((pred_T.P.seas_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.P.seas[i]  = sqrt(mean((log(pred_T.P.seas_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.P.seas[i]  = mean(abs(pred_T.P.seas_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  # Build new crossbasis for Prec
  cbPrec_val <- crossbasis(bio.matrix_val_i$cumPrecweek, lag=tlag,
                           argvar=list(fun="bs", degree=2, knots=pknots),
                           arglag=list(knots=c(1,4)),
                           group=bio.matrix_val_i$Region)
  
  cbPrec_val_df <- as.data.frame(cbPrec_val)
  colnames(cbPrec_val_df) <- paste0("cbp", seq_len(ncol(cbPrec_val_df)))
  
  bio.matrix_val_i_T.P.seas <- cbind(bio.matrix_val_i, cbTemp_val_df, cbPrec_val_df, bseas_val_df)
  
  # predict
  pred_T.P.seas<- predict(mod_T.P.seas_i, newdata = bio.matrix_val_i_T.P.seas, type = "response")
  
  RMSE_val_df$T.P.seas[i] = sqrt(mean((pred_T.P.seas[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.P.seas[i]  = sqrt(mean((log(pred_T.P.seas[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.P.seas[i]  = mean(abs(pred_T.P.seas[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)

  # 6. T.seas.trend ----
  
  # and cbind
  bio.matrix_fit_i_ext_T.seas.trend <- cbind(bio.matrix_fit_i, cbTemp_df, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbind(cbTemp_df, btrend_df)), collapse = " + ")))
  
  
  # model definition 
  # with temperature, spline trend & seasonality
  mod_T.seas.trend_i <- gnm(form,
                          eliminate=as.factor(ID), #grouping factor
                          data=bio.matrix_fit_i_ext_T.seas.trend,
                          family=quasipoisson(link="log"),
                          na.action= "na.exclude")
  
  # assess fitting error
  
  pred_mod_T.seas.trend_i = mod_T.seas.trend_i$fitted
  obs_fit_i = mod_T.seas.trend_i$fitted + mod_T.seas.trend_i$residuals
  
  RMSE_fit_df$T.seas.trend[i] = sqrt(mean((pred_mod_T.seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.seas.trend[i]  = sqrt(mean((log(pred_mod_T.seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.seas.trend[i]  = mean(abs(pred_mod_T.seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.seas.trend <- cbind(bio.matrix_val_i, cbTemp_val_df, btrend_val_df)
  
  # predict
  pred_T.seas.trend <- predict(mod_T.seas.trend_i, newdata = bio.matrix_val_i_T.seas.trend, type = "response")
  
  RMSE_val_df$T.seas.trend[i] = sqrt(mean((pred_T.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.seas.trend[i]  = sqrt(mean((log(pred_T.seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.seas.trend[i]  = mean(abs(pred_T.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  # 7. T.P.seas.trend ----
  
  # define cbind
  bio.matrix_fit_i_ext_T.P.seas.trend <- cbind(bio.matrix_fit_i, cbTemp_df, cbPrec_df, btrend_df) 
  
  # Build formula with those columns
  form <- as.formula(paste("eggs ~",
                           paste(colnames(cbind(cbTemp_df, cbPrec_df, btrend_df)), collapse = " + ")))
  
  # model definition 
  # with temperature, precipitaiton, with spline trend & seasonality
  mod_T.P.seas.trend_i <- gnm(form,
                            eliminate=as.factor(ID), #grouping factor
                            data=bio.matrix_fit_i_ext_T.P.seas.trend,
                            family=quasipoisson(link="log"),
                            na.action= "na.exclude")
  
  # assess fitting error
  
  pred_mod_T.P.seas.trend_i = mod_T.P.seas.trend_i$fitted
  obs_fit_i = mod_T.P.seas.trend_i$fitted + mod_T.P.seas.trend_i$residuals
  
  RMSE_fit_df$T.P.seas.trend[i] = sqrt(mean((pred_mod_T.P.seas.trend_i - obs_fit_i)^2, na.rm = T))
  RMSLE_fit_df$T.P.seas.trend[i]  = sqrt(mean((log(pred_mod_T.P.seas.trend_i+1) - log(obs_fit_i+1))^2, na.rm = T))
  MAE_fit_df$T.P.seas.trend[i]  = mean(abs(pred_mod_T.P.seas.trend_i - obs_fit_i), na.rm = T)
  
  ## assess model over validation sample ----
  
  bio.matrix_val_i_T.P.seas.trend <- cbind(bio.matrix_val_i, cbTemp_val_df, cbPrec_val_df, btrend_val_df)
  
  # predict
  pred_T.P.seas.trend <- predict(mod_T.P.seas.trend_i, newdata = bio.matrix_val_i_T.P.seas.trend, type = "response")
  
  RMSE_val_df$T.P.seas.trend[i] = sqrt(mean((pred_T.P.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)])^2, na.rm = T))
  RMSLE_val_df$T.P.seas.trend[i]  = sqrt(mean((log(pred_T.P.seas.trend[-(obs_val_na)]+1) - log(obs_val[-(obs_val_na)]+1))^2, na.rm = T))
  MAE_val_df$T.P.seas.trend[i]  = mean(abs(pred_T.P.seas.trend[-(obs_val_na)] - obs_val[-(obs_val_na)]), na.rm = T)
  
  
  ## Other metrics ----
  # Overdispersion parameters and residual deviance
  
  Dispersion_df$seas[i] = summary(mod_seas_i)$dispersion
  Dispersion_df$seas.trend[i] = summary(mod_seas.trend_i)$dispersion
  Dispersion_df$Temp[i] = summary(mod_T_i)$dispersion
  Dispersion_df$T.seas[i] = summary(mod_T.seas_i)$dispersion
  Dispersion_df$T.P.seas[i] = summary(mod_T.P.seas_i)$dispersion
  Dispersion_df$T.seas.trend[i] = summary(mod_T.seas.trend_i)$dispersion
  Dispersion_df$T.P.seas.trend[i] = summary(mod_T.P.seas.trend_i)$dispersion
  
  ResDev_df$seas[i] = summary(mod_seas_i)$deviance
  ResDev_df$seas.trend[i] = summary(mod_seas.trend_i)$deviance
  ResDev_df$Temp[i] = summary(mod_T_i)$deviance
  ResDev_df$T.seas[i] = summary(mod_T.seas_i)$deviance
  ResDev_df$T.P.seas[i] = summary(mod_T.P.seas_i)$deviance
  ResDev_df$T.seas.trend[i] = summary(mod_T.seas.trend_i)$deviance
  ResDev_df$T.P.seas.trend[i] = summary(mod_T.P.seas.trend_i)$deviance
}
toc()

# Plots ----
# fitting vs + validation, different metrics
# think of asdding wilcoson test

## RMSE ----

#tidy R
# pivot_longer(RMSE_val_df, !rep)

RMSE_fit_ldf <- reshape2::melt(RMSE_fit_df,
                               id.vars = "rep",
                               variable.name = "model",
                               value.name = "RMSE")%>%
  mutate(dataset = "Fitting")

RMSE_val_ldf <- reshape2::melt(RMSE_val_df,
                               id.vars = "rep",
                               variable.name = "model",
                               value.name = "RMSE") %>%
  mutate(dataset = "Validation")

RMSE_fit <- ddply(RMSE_fit_ldf, "model", summarise, 
                  mean=mean(RMSE),
                  median=median(RMSE),
                  P10=quantile(x = RMSE, probs = 0.1),
                  P90=quantile(x = RMSE, probs = 0.9))

RMSE_val <- ddply(RMSE_val_ldf, "model", summarise, 
                 mean=mean(RMSE),
                 median=median(RMSE),
                 P10=quantile(x = RMSE, probs = 0.1),
                 P90=quantile(x = RMSE, probs = 0.9))

p_RMSE_fit <- ggplot(RMSE_fit_ldf, aes(x = reorder(model, RMSE, FUN = median), y = RMSE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Fitting")+
  theme_test() +
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )
 
p_RMSE_val <- ggplot(RMSE_val_ldf, aes(x = reorder(model, RMSE, FUN = median), y = RMSE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Validation")+
  theme_test()+
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

ggarrange(p_RMSE_fit + scale_y_continuous(trans='log10'),
          p_RMSE_val,
          ncol = 1, nrow = 2)

## RMSLE ----

RMSLE_fit_ldf <- reshape2::melt(RMSLE_fit_df,
                               id.vars = "rep",
                               variable.name = "model",
                               value.name = "RMSLE")%>%
  mutate(dataset = "Fitting")

RMSLE_val_ldf <- reshape2::melt(RMSLE_val_df,
                               id.vars = "rep",
                               variable.name = "model",
                               value.name = "RMSLE") %>%
  mutate(dataset = "Validation")

RMSLE_fit <- ddply(RMSLE_fit_ldf, "model", summarise, 
                  mean=mean(RMSLE),
                  median=median(RMSLE),
                  P10=quantile(x = RMSLE, probs = 0.1),
                  P90=quantile(x = RMSLE, probs = 0.9))

RMSLE_val <- ddply(RMSLE_val_ldf, "model", summarise, 
                  mean=mean(RMSLE),
                  median=median(RMSLE),
                  P10=quantile(x = RMSLE, probs = 0.1),
                  P90=quantile(x = RMSLE, probs = 0.9))

p_RMSLE_fit <- ggplot(RMSLE_fit_ldf, aes(x = reorder(model, RMSLE, FUN = median), y = RMSLE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Fitting")+
  theme_test() +
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

p_RMSLE_val <- ggplot(RMSLE_val_ldf, aes(x = reorder(model, RMSLE, FUN = median), y = RMSLE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Validation")+
  theme_test()+
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

ggarrange(p_RMSLE_fit, p_RMSLE_val, ncol = 1, nrow = 2)

## MAE ----

MAE_fit_ldf <- reshape2::melt(MAE_fit_df,
                                id.vars = "rep",
                                variable.name = "model",
                                value.name = "MAE")%>%
  mutate(dataset = "Fitting")

MAE_val_ldf <- reshape2::melt(MAE_val_df,
                                id.vars = "rep",
                                variable.name = "model",
                                value.name = "MAE") %>%
  mutate(dataset = "Validation")

MAE_fit <- ddply(MAE_fit_ldf, "model", summarise, 
                   mean=mean(MAE),
                   median=median(MAE),
                   P10=quantile(x = MAE, probs = 0.1),
                   P90=quantile(x = MAE, probs = 0.9))

MAE_val <- ddply(MAE_val_ldf, "model", summarise, 
                   mean=mean(MAE),
                   median=median(MAE),
                   P10=quantile(x = MAE, probs = 0.1),
                   P90=quantile(x = MAE, probs = 0.9))

p_MAE_fit <- ggplot(MAE_fit_ldf, aes(x = reorder(model, MAE, FUN = median), y = MAE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Fitting")+
  theme_test() +
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

p_MAE_val <- ggplot(MAE_val_ldf, aes(x = reorder(model, MAE, FUN = median), y = MAE, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Validation")+
  theme_test()+
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

ggarrange(p_MAE_fit + coord_transform(y='log10'), p_MAE_val, ncol = 1, nrow = 2)

## Residual deviance and dispersion ----

ResDev_ldf <- reshape2::melt(ResDev_df,
                              id.vars = "rep",
                              variable.name = "model",
                              value.name = "ResDev")

Dispersion_ldf <- reshape2::melt(Dispersion_df,
                             id.vars = "rep",
                             variable.name = "model",
                             value.name = "Dispersion")

ResDev <- ddply(ResDev_ldf, "model", summarise, 
                 mean=mean(ResDev),
                 median=median(ResDev),
                 P10=quantile(x = ResDev, probs = 0.1),
                 P90=quantile(x = ResDev, probs = 0.9))

Dispersion <- ddply(Dispersion_ldf, "model", summarise, 
                 mean=mean(Dispersion),
                 median=median(Dispersion),
                 P10=quantile(x = Dispersion, probs = 0.1),
                 P90=quantile(x = Dispersion, probs = 0.9))


p_ResDev <- ggplot(ResDev_ldf, aes(x = reorder(model, ResDev, FUN = median), y = ResDev, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("ResDev")+
  theme_test() +
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

p_Dispersion <- ggplot(Dispersion_ldf, aes(x = reorder(model, Dispersion, FUN = median), y = Dispersion, fill = model))+
  geom_boxplot(alpha=0.5, position="identity", width=0.15)+
  ggtitle("Dispersion")+
  theme_test()+
  theme(axis.text.x = element_text(angle = +90))+
  stat_summary(
    geom = "text",
    fun = function(y) quantile(y, probs = c(0.5), na.rm = TRUE),
    aes(label = formatC(..y.., format = "e", digits = 2)),
    position = position_nudge(x = 0.35),
    size = 2.5
  )

ggarrange(p_ResDev + scale_y_continuous(trans='log10'),
          p_Dispersion + scale_y_continuous(trans='log10'),
          ncol = 1, nrow = 2)

