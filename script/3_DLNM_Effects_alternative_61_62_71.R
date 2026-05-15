## Alternative to 6.1 -6.2  ----

# 6.1 rebuilt model for prediction (test years) ----

cbTemp_df <- as.data.frame(cbTemp)
colnames(cbTemp_df) <- paste0("cbt", seq_len(ncol(cbTemp_df)))

cbPrec_df <- as.data.frame(cbPrec)
colnames(cbPrec_df) <- paste0("cbp", seq_len(ncol(cbPrec_df)))

bseas <- ns(bio.matrix_sel$weekn, df = 4)
bseas_df <- as.data.frame(bseas)
colnames(bseas_df) <- paste0("bs", seq_len(ncol(bseas_df)))

years_df <- data.frame("year" = bio.matrix_sel$year)

ID_df <- data.frame("ID" = bio.matrix_sel$ID)

matrix_sel_pred <- cbind(cbTemp_df, cbPrec_df, bseas_df, years_df) #without eggs

form <- as.formula(paste("eggs ~", paste(colnames(matrix_sel_pred), collapse = " + ")))

bio.matrix_sel_pred <- cbind(bio.matrix_sel$eggs, matrix_sel_pred, ID_df)
names(bio.matrix_sel_pred) <- c("eggs", names(matrix_sel_pred), "ID")

mod_T.P.seas_pred <- gnm(form,
                         eliminate=as.factor(ID),
                         data = bio.matrix_sel_pred,
                         family = quasipoisson(link="log"),
                         na.action = "na.exclude")

# Prediction function  
predict_mod <- function(newdata, mod, 
                        tlag = 8) {
  cbT_new <- crossbasis(newdata$medianTweek, lag = tlag,
                        argvar = list(fun = "bs", degree = 2, knots = tknots),
                        arglag = list(knots = c(1, 4)),
                        group  = newdata$Region)
  
  cbP_new <- crossbasis(newdata$cumPrecweek, lag = tlag,
                        argvar = list(fun = "bs", degree = 2, knots = pknots),
                        arglag = list(knots = c(1, 4)),
                        group  = newdata$Region)
  
  #### as Data Frame
  
  cbTemp_df <- as.data.frame(cbT_new)
  colnames(cbTemp_df) <- paste0("cbt", seq_len(ncol(cbTemp_df)))
  
  cbPrec_df <- as.data.frame(cbP_new )
  colnames(cbPrec_df) <- paste0("cbp", seq_len(ncol(cbPrec_df)))
  
  bseas <- ns(newdata$weekn, df=4)
  
  bseas_df <- as.data.frame(bseas)
  colnames(bseas_df) <- paste0("bs", seq_len(ncol(bseas_df)))
  
  if(sum(names(newdata)=="ID")==0){
    newdata$ID =NA
  }
  
  yearsID_df <- data.frame("year" = newdata$year,
                           "ID" = newdata$ID)
  
  bio.matrix_val = cbind(cbTemp_df, cbPrec_df, bseas_df, yearsID_df)
  
  prediction = predict(mod,
                       newdata = bio.matrix_val,
                       type = "response") 
}

pred_years_known <- predict_mod(newdata = bio.matrix_test_years_known, 
                                            mod = mod_T.P.seas_pred)

# 6.2 rebuilt model for prediction (test years): no ID ----

mod_T.P.seas_pred_noID <- gnm(form,
                              data = bio.matrix_sel_pred,
                              family = quasipoisson(link="log"),
                              na.action = "na.exclude")

pred_regions <- predict_mod(newdata = bio.matrix_test_regions, 
                                mod = mod_T.P.seas_pred_noID)

## 7.1 Predict on spatial grid ----

# env.matrix has no trap IDs — use avg_intercept for all locations
# Region column needed for crossbasis group argument — check if present
cat("Region column present:", "Region" %in% names(env.matrix), "\n")

# If no Region column, create a dummy — crossbasis group just breaks lag continuity
# at region boundaries; a single group is valid for a spatial grid
if (!"Region" %in% names(env.matrix)) {
  env.matrix <- env.matrix %>% dplyr::mutate(Region = "grid")
  cat("Region column added as dummy\n")
}

env.matrix$eggsPred <- predict_mod(newdata = env.matrix, 
                                   mod = mod_T.P.seas_pred_noID)
