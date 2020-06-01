# Regression of treelet modelled food group intake on periodontal disease measures from NHANES

library(tidyverse)
library(betareg) # For Beta regression
library(Formula) # For multipart formulas used in Beta regression
library(lmtest) # For likelihood ratio tests of Beta regressions
library(quantreg) # For quantile regression
library(lqr) # For logistic quantile regression
library(cdfquantreg) # For logistic quantile regression (more recent package)
library(splines)
library(broom)

cat("\n Beta regression of CAL by food groups")


fg1 <- betareg(formula = prop_CAL_sites3mm_beta ~ RIDAGEYR + RIAGENDR + 
                 SMQ020 + diabetes + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + 
                 TC7 + TC8 | RIDAGEYR + RIAGENDR + SMQ020 + diabetes, data = nh_grps)

# Age-TC interaction (selected TCs)
fg2 <- update(fg1, .~. + RIDAGEYR:TC1 + RIDAGEYR:TC7 |., data = nh_grps)

# TCs in precision model
fg3 <- update(fg1, .~.|. + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + TC7 + TC8, data = nh_grps)

# Quadratic fit with TCs in both mean and variance models
fg4 <- update(fg3, .~. + I(TC1^2) + I(TC2^2) + I(TC3^2) + I(TC4^2) + I(TC5^2) + I(TC6^2) + I(TC7^2) + I(TC8^2)  |. + I(TC1^2) + I(TC2^2) + I(TC3^2) + I(TC4^2) + I(TC5^2) + I(TC6^2) + I(TC7^2) + I(TC8^2), data = nh_grps)

# Deciles for mean model
fg5 <- betareg(formula = prop_CAL_sites3mm_beta ~ RIDAGEYR + RIAGENDR + 
                 SMQ020 + diabetes + TC1_dec + TC2_dec + TC3_dec + TC4_dec + TC5_dec + TC6_dec + 
                 TC7_dec + TC8_dec | RIDAGEYR + RIAGENDR + SMQ020 + diabetes, data = nh_grps)

# Deciles for mean and precision model
fg6 <- update(fg5, .~.|. + TC1_dec + TC2_dec + TC3_dec + TC4_dec + TC5_dec + TC6_dec + 
                TC7_dec + TC8_dec, data = nh_grps)

# Create prediction frame
# Take quantiles and stack
tc_qtx <- StackQuantiles(food_groups_tc_scores, probs = seq(0.02, 0.98, by = 0.01))

# Prediction for average aged, male, non-smoker, non-diabetic, average energy intake
cal_pred_frame_grps <- cbind(RIDAGEYR = mean(nh_grps$RIDAGEYR), 
                             RIAGENDR = 1, 
                             SMQ020 = 2, 
                             diabetes = F,
                             KCAL = mean(nh_grps$KCAL),
                             tc_qtx) %>% 
  as_tibble()


# Logistic regression
lg1 <- glm(prop_CAL_sites3mm ~ RIDAGEYR + RIAGENDR + 
             SMQ020 + diabetes + TC1_dec + TC2_dec + TC3_dec + TC4_dec + TC5_dec + TC6_dec + TC7_dec + TC8_dec, family = binomial(link = "logit"), weights = CAL_sites_assessed,  data = nh_grps)

# Actual and predicted means by TC1 decile
TC1_means <- nh_grps %>% 
  cbind(predicted_logistic = predict(lg1, type = "response"),
        predicted_beta = predict(fg6)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual mean" = mean(prop_CAL_sites3mm_beta),
            "Predicted mean: beta regression (decile)" = mean(predicted_beta),
            "Predicted mean: logistic regression (decile)" = mean(predicted_logistic)) %>% 
  gather(Measurement, value, -TC_decile)


cat("\n Robust quantile regression by food groups")
# Specify model for the median 
rob_mod <- ~ 1 + RIDAGEYR + RIAGENDR + 
  SMQ020 + diabetes + TC1_dec + TC2_dec + TC3_dec + TC4_dec + TC5_dec + TC6_dec + 
  TC7_dec + TC8_dec

# Model the median
rob1 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                x = model.matrix(rob_mod, data = nh_grps),
                p = 0.5)

# Model the lower quartile
rob1_lower <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod, data = nh_grps),
                      p = 0.25) 

# Model the upper quartile
rob1_upper <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod, data = nh_grps),
                      p = 0.75) 

#model.matrix(rob_mod, data = nh)[2,] %*% rob1$beta
#invlogit(rob1$fitted.values[2])
#invlogit(as.numeric(rob_pred_frame) %*% rob1$beta)

# Predict at the mean age, male, non-smoker, non-diabetic, all TCs at decile 5
rob_pred_frame <- model.matrix(rob_mod, data = nh_grps)[1,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(RIDAGEYR = 0,
         RIAGENDR = 1,
         SMQ020 = 2, 
         diabetesTRUE = 0) %>% 
  mutate_at(vars(matches("_dec")), function(.){0}) %>% 
  mutate_at(vars(matches("_dec5")), function(.){1}) 

# Median predictions for TC1 from robust regression
TC1_medians <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob1$fitted.values),
        predicted_upper = invlogit(rob1_upper$fitted.values),
        predicted_lower = invlogit(rob1_lower$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)

# Extract coefficients for median model
rob1_coef <- bind_cols(Term = colnames(model.matrix(rob_mod, data = nh_grps)),
                       as_tibble(rob1$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)#,
         # Why are these values so high?
         #"Predicted proportion" = round(invlogit(ifelse(Term == "(Intercept)", Estimate[1], Estimate[1] + Estimate)), 2)
  ) %>% 
  rename(Sig = V1)


cat("\n Logistic CDF quantile regression")
lrq1 <- cdfquantreg(prop_CAL_sites3mm ~ RIDAGEYR + RIAGENDR + SMQ020 + diabetes + KCAL 
                    + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + TC7 + TC8  
                    | RIDAGEYR + RIAGENDR + SMQ020 + diabetes,
                    fd = "logit", sd = "logistic", data = as.data.frame(nh_grps))

# # Decile model slow but similar distribution
# lrq1 <- cdfquantreg(prop_CAL_sites3mm ~ RIDAGEYR + RIAGENDR + SMQ020 + diabetes
#                       + TC1_dec + TC2_dec + TC3_dec + TC4_dec + TC5_dec + TC6_dec + TC7_dec + TC8_dec 
#                       | RIDAGEYR + RIAGENDR + SMQ020 + diabetes,
#        fd = "logit", sd = "logistic", data = as.data.frame(nh_grps))

# Median predictions for TC1 from CDF regression
TC1_cdf_medians <- nh_grps %>% 
  cbind(predicted_cdf = predict(lrq1),
        predicted_robust = invlogit(rob1$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted median: CDF regression (decile)" = median(predicted_cdf),
            "Predicted median: robust regression (decile)" = median(predicted_robust)
  ) %>% 
  gather(Measurement, value, -TC_decile)


# Median predictions for TC7 from robust regression
TC7_medians <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob1$fitted.values),
        predicted_upper = invlogit(rob1_upper$fitted.values),
        predicted_lower = invlogit(rob1_lower$fitted.values)) %>% 
  group_by(TC_decile = TC7_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)


## Predictions from robust regression ##

# TC1 at decile 6
rob_pred_frame %>% 
  mutate(TC1_dec5 = 0,
         TC1_dec6 = 1) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

# TC1 at decile 1
tc1_dec1 <- rob_pred_frame %>% 
  mutate(TC1_dec5 = 0) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100


# TC1 at decile 10
tc1_dec10 <- rob_pred_frame %>% 
  mutate(TC1_dec5 = 0,
         TC1_dec10 = 1) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100 

# TC7 at decile 1
tc7_dec1 <- rob_pred_frame %>% 
  mutate(TC7_dec5 = 0) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

# TC7 at decile 5
rob_pred_frame %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

# TC1 at decile 10
tc7_dec10 <- rob_pred_frame %>% 
  mutate(TC7_dec5 = 0,
         TC7_dec10 = 1) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

# Age at mean
invlogit(as.numeric(rob_pred_frame) %*% rob1$beta)*100

# Age at mean + xx years
rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == round(mean(nhanes$RIDAGEYR)) + 17])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

# Age at given years
age47 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 47])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

age46 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 46])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

age45 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 45])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100


age61 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 61])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100


# Age at 5th percentile
age_5pc <- rob_pred_frame %>% 
  mutate(RIDAGEYR = quantile(nh_grps$RIDAGEYR, 0.05)) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()

# Age at 95th percentile
age_95pc <- rob_pred_frame %>% 
  mutate(RIDAGEYR = quantile(nh_grps$RIDAGEYR, 0.95)) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()


# Diabetes effect
rob_pred_frame %>% 
  mutate(diabetesTRUE = 1) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()

# Smoking effect
rob_pred_frame %>% 
  mutate(SMQ020 = 1) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()

# Gender effect
rob_pred_frame %>% 
  mutate(RIAGENDR = 2) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()



## Displaying daily intake of food groups by TC score ##

# Join the daily intake values to the individual characteristics (including the TC scores and deciles)
daily_intake <- nh_grps %>%
  inner_join(food_grps_grms_per_day %>% 
               gather(grp_code, grms, -SEQN), by = "SEQN")

# Function to extract daily intake by food groups and deciles for a given TC
# Args:
# daily_intake = tibble with daily intake in grams
# TC = number of TC to extract
# decile_column = name of column with deciles to group by
ExtractTCIntake <- function(daily_intake, TC, decile_column){
  dec_col <- enquo(decile_column)
  daily_intake %>% 
    # Select only those food groups with a non zero score for the TC of interest
    inner_join(food_groups_loadings %>% 
                 filter(Component == paste0("TC", TC), Value != 0) %>%  
                 mutate(grp_code = as.character(Variable)),
               by = "grp_code") %>% 
    # Calculate summary statistics for each decile and food group
    group_by(!!dec_col, Variable) %>% 
    # Proportion with any intake 
    summarise(any_intake = round(sum(grms>0)/length(grms), digits = 2),
              # Median intake among those with any intake
              grms_intake = round(median(grms[grms>0]))) %>% 
    ungroup()
}
# # TC 1
# ExtractTCIntake(daily_intake, 1, TC1_dec) %>% 
#   transmute(TC1_dec, `Food group` = Variable, intake = paste0(any_intake, " (", grms_intake, ")")) %>% 
#   spread(TC1_dec, intake) %>% 
#   knitr::kable(caption = "Daily intake of food groups by TC1 decile. Proportion with any intake (median among those with any intake, g/day.")


# Select variables correlated with TC1
# Kendall's tau because skewed TC1 score
cor_tc1 <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC1, ., method = "kendall"), 3)) 

cor_tc1_pearson <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC1, ., method = "pearson"), 3)) 


# Select variables correlated with TC7
cor_tc7 <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC7, ., method = "kendall"), 3))

cor_tc7_pearson <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC7, ., method = "pearson"), 3)) 

ggplot(nh_grps, aes(x = THEO, y = TC1)) + geom_density2d()


# Intake of carbohydrates from food and drink by TC7
food_drink_intake <- nh_grps %>% 
  select(SEQN, TC7, TC7_dec) %>% 
  inner_join(
    food_grps_per_day %>% 
      mutate(Intake_Carbs = if_else(grp_code %in% (bevs_included), "Drink", "Food")) %>% 
      group_by(SEQN, Intake_Carbs) %>% 
      summarise(CARB = sum(CARB)) %>% 
      spread(Intake_Carbs, CARB, sep = "_", fill = 0),
    by = "SEQN")


# Robust regression including carbohydrate intake from food and drink
rob_mod2 <- update(rob_mod, ~. + Intake_Carbs_Drink + Intake_Carbs_Food)
# rob_mod2 <- update(rob_mod, ~. + Intake_Carbs_Drink + Intake_Carbs_Food + I(Intake_Carbs_Drink^2) + I(Intake_Carbs_Food^2))

# Setup model matrix
rob2_matrix <- food_drink_intake %>% 
  select(SEQN, Intake_Carbs_Drink, Intake_Carbs_Food) %>% 
  mutate_at(vars(matches("Intake")), scale) %>% 
  inner_join(nh_grps, by = "SEQN") %>% 
  model.matrix(rob_mod2, data = .)

# Fit model
rob2 <- Log.lqr(rob2_matrix, y = nh_grps$prop_CAL_sites3mm, p = 0.5)

# Extract coefficients 
rob2_coef <- bind_cols(Term = colnames(rob2_matrix),
                       as_tibble(rob2$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)


# Median predictions for TC1 from robust regression
TC1_medians_rob2 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob2$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted median: robust regression (decile)" = median(predicted_robust)) %>% 
  gather(Measurement, value, -TC_decile)

TC7_medians_rob2 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob2$fitted.values)) %>% 
  group_by(TC_decile = TC7_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted median: robust regression (decile)" = median(predicted_robust)) %>% 
  gather(Measurement, value, -TC_decile)

cat("\n Treelet analysis on food groups completed")


# Calculate other periodontal measures by TC1 and TC7 deciles #

TC1_perio <- nh_grps %>% 
  group_by(TC1_dec) %>% 
  summarise(CAL_mouth = mean(CAL_mouth_mean),
            CAL_mouth_se =
              sd(CAL_mouth_mean)/sqrt(length(CAL_mouth_mean)),
            PD_mouth = mean(PD_mouth_mean),
            PD_mouth_se = sd(PD_mouth_mean)/sqrt(length(CAL_mouth_mean))) %>% 
  mutate_at(vars(CAL_mouth, PD_mouth), formatC, digits = 1, format = "f") %>% mutate_at(vars(CAL_mouth_se, PD_mouth_se), formatC, digits = 2, format = "f") %>% 
  transmute(Variable = paste("TC1 Decile", TC1_dec),
        `Mean CAL (mm)` = paste0(CAL_mouth, " (", CAL_mouth_se, ")"),
         `Mean PPD (mm)` = paste0(PD_mouth, " (", PD_mouth_se, ")"))

TC7_perio <- nh_grps %>% 
  group_by(TC7_dec) %>% 
  summarise(CAL_mouth = mean(CAL_mouth_mean),
            CAL_mouth_se =
              sd(CAL_mouth_mean)/sqrt(length(CAL_mouth_mean)),
            PD_mouth = mean(PD_mouth_mean),
            PD_mouth_se = sd(PD_mouth_mean)/sqrt(length(CAL_mouth_mean))) %>% 
  mutate_at(vars(CAL_mouth, PD_mouth), formatC, digits = 1, format = "f") %>% mutate_at(vars(CAL_mouth_se, PD_mouth_se), formatC, digits = 2, format = "f") %>% 
  transmute(Variable = paste("TC7 Decile", TC7_dec),
            `Mean CAL (mm)` = paste0(CAL_mouth, " (", CAL_mouth_se, ")"),
            `Mean PPD (mm)` = paste0(PD_mouth, " (", PD_mouth_se, ")"))



# Robust regression excluding those with low tooth count

# Model the median
rob3 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm[nh_grps$tooth_count>=20], 
                x = model.matrix(rob_mod, data = nh_grps[nh_grps$tooth_count>=20,]),
                p = 0.5)

# Extract coefficients for median model
rob3_coef <- bind_cols(Term = colnames(model.matrix(rob_mod, data = nh_grps[nh_grps$tooth_count>=20,])),
                       as_tibble(rob3$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)



## Robust regression adding total energy intake ##
# Multivariate nutrient density model from Willett et al., 1997.

rob4 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                x = model.matrix(update(rob_mod, ~.+KCAL_raw), data = nh_grps),
                p = 0.5)

rob4_coef <- bind_cols(Term = colnames(model.matrix(update(rob_mod, ~.+KCAL_raw), data = nh_grps)),
                       as_tibble(rob4$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)



## Robust regression with three category smoking status ##

rob5 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                x = model.matrix(update(rob_mod, ~.-SMQ020 + smoking), data = nh_grps),
                p = 0.5)

rob5_lower = Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                     x = model.matrix(update(rob_mod, ~.-SMQ020 + smoking), data = nh_grps),
                     p = 0.25)

rob5_upper = Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                     x = model.matrix(update(rob_mod, ~.-SMQ020 + smoking), data = nh_grps),
                     p = 0.75)

rob5_coef <- bind_cols(Term = colnames(model.matrix(update(rob_mod, ~.-SMQ020 + smoking), data = nh_grps)),
                       as_tibble(rob5$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)


# Median predictions for TC1 and TC7 from robust regression

TC1_medians_rob5 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob5$fitted.values),
        predicted_upper = invlogit(rob5_upper$fitted.values),
        predicted_lower = invlogit(rob5_lower$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)

TC7_medians_rob5 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob5$fitted.values),
        predicted_upper = invlogit(rob5_upper$fitted.values),
        predicted_lower = invlogit(rob5_lower$fitted.values)) %>% 
  group_by(TC_decile = TC7_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)



# Format all the TC estimates 
rob5_summary <- nh_grps %>% 
  select(SEQN, matches("^TC[0-9]{1}_dec"),
         CAL_mouth_mean, PD_mouth_mean) %>% 
  cbind(predicted_robust = invlogit(rob5$fitted.values)) %>% 
  pivot_longer(cols = matches("^TC[0-9]{1}_dec"), 
               names_to = "TC", 
               values_to = "TC_decile") %>% 
  group_by(TC, TC_decile) %>% 
  summarise(`Median sites with CAL` = perc(median(predicted_robust)),
            CAL_mouth = mean(CAL_mouth_mean),
            CAL_mouth_se =
              sd(CAL_mouth_mean)/sqrt(length(CAL_mouth_mean)),
            PD_mouth = mean(PD_mouth_mean),
            PD_mouth_se = sd(PD_mouth_mean)/sqrt(length(CAL_mouth_mean))) %>%
  mutate_at(vars(CAL_mouth, PD_mouth), formatC, digits = 1, format = "f") %>% 
  mutate_at(vars(CAL_mouth_se, PD_mouth_se), formatC, digits = 2, format = "f") %>% 
  mutate(`Mean CAL (mm)` = paste0(CAL_mouth, " (", CAL_mouth_se, ")"),
         `Mean PPD (mm)` = paste0(PD_mouth, " (", PD_mouth_se, ")")) %>% 
  unite(col = Term, TC, TC_decile, sep = "", remove = F) %>% 
  left_join(rob5_coef %>% 
  mutate(`Odds Ratio` = FormatOddsRatio(Estimate, `Std. Error`),
         P = format.pval(`Pr(>|z|)`, digits = 1, eps = 0.001)),
  by = "Term") %>% 
  ungroup() %>% 
  mutate_at(vars(`Median sites with CAL`, P), ~if_else(is.na(.), "", .)) %>% 
  mutate(`Odds Ratio` = if_else(is.na(`Odds Ratio`), "1.00", `Odds Ratio`)) %>% 
  mutate(TC = str_extract(TC, "^TC[0-9]{1}")) %>% 
  
  select(TC, TC_decile, `Odds Ratio`, P, `Median sites with CAL`, `Mean CAL (mm)`, `Mean PPD (mm)`) %>% 
  rename(`Median % sites with CAL $\\geq$ 3mm` = `Median sites with CAL`,
         `Mean CAL (± SE) mm` = `Mean CAL (mm)`,
         `Mean PPD  (± SE) mm` = `Mean PPD (mm)`) %>% 
  # rename(`Odds Ratio (median sites with CAL)` = `Odds Ratio`) %>% 
  pivot_longer(cols = -one_of("TC", "TC_decile"), names_to = "Metric") %>% 
  pivot_wider(names_from = TC_decile)


nh_grps %>% 
  select(SEQN, matches("^TC[0-9]{1}_dec"),
         CAL_mouth_mean, PD_mouth_mean) %>% 
  cbind(predicted_robust = invlogit(rob5$fitted.values)) %>% 
  pivot_longer(cols = matches("^TC[0-9]{1}_dec"), 
               names_to = "TC", 
               values_to = "TC_decile") %>% 
  group_by(TC, TC_decile) %>% 
  summarise(median_sites_with_CAL = median(predicted_robust)) %>% 
  group_by(TC) %>% 
  summarise(min_sites = min(median_sites_with_CAL),
            max_sites = max(median_sites_with_CAL),
            range_sites = max_sites - min_sites)
  
# Predictions for various groups

rob5_pred_frame <-  model.matrix(update(rob_mod, ~.-SMQ020 + smoking), data = nh_grps)[1,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(RIDAGEYR = 0,
         RIAGENDR = 1,
         smokingFormer = 0,
         smokingCurrent = 0, 
         diabetesTRUE = 0) %>% 
  mutate_at(vars(matches("_dec")), function(.){0}) %>% 
  mutate_at(vars(matches("_dec5")), function(.){1}) 

# Base prediction
rob5_pred_frame %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100

# Age at given years
rob5_age46 <- rob5_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 46])) %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100

rob5_age61 <- rob5_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 61])) %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100


# Diabetes effect
rob5_pred_frame %>% 
  mutate(diabetesTRUE = 1) %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100

# Smoking effect
rob5_pred_frame %>% 
  mutate(smokingCurrent = 1) %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100

# Gender effect
rob5_pred_frame %>% 
  mutate(RIAGENDR = 2) %>% 
  as.numeric() %*% rob5$beta %>% 
  invlogit()*100










## Robust regression excluding current smokers ##

rob6 <- Log.lqr(y = filter(nh_grps, smoking != "Current")$prop_CAL_sites3mm, 
                x = model.matrix(rob_mod, data = filter(nh_grps, smoking != "Current")),
                p = 0.5)

# Model the lower quartile
rob6_lower <- Log.lqr(y = filter(nh_grps, smoking != "Current")$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod, data = filter(nh_grps, smoking != "Current")),
                      p = 0.25) 

# Model the upper quartile
rob6_upper <- Log.lqr(y = filter(nh_grps, smoking != "Current")$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod, data = filter(nh_grps, smoking != "Current")),
                      p = 0.75) 

# Extract coefficients for median model
rob6_coef <- bind_cols(Term = colnames(model.matrix(rob_mod, data = filter(nh_grps, smoking != "Current"))),
                       as_tibble(rob6$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)


TC1_medians_rob6 <- filter(nh_grps, smoking != "Current") %>% 
  cbind(predicted_robust = invlogit(rob6$fitted.values),
        predicted_upper = invlogit(rob6_upper$fitted.values),
        predicted_lower = invlogit(rob6_lower$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)

TC7_medians_rob6 <- filter(nh_grps, smoking != "Current") %>% 
  cbind(predicted_robust = invlogit(rob6$fitted.values),
        predicted_upper = invlogit(rob6_upper$fitted.values),
        predicted_lower = invlogit(rob6_lower$fitted.values)) %>% 
  group_by(TC_decile = TC7_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -TC_decile)


## Robust regression with three category smoking status and food/drink carbohydrate correction ##

rob_mod7 <- update(rob_mod2, ~. -SMQ020 + smoking)

# Setup model matrix
rob7_matrix <- food_drink_intake %>% 
  select(SEQN, Intake_Carbs_Drink, Intake_Carbs_Food) %>% 
  mutate_at(vars(matches("Intake")), scale) %>% 
  inner_join(nh_grps, by = "SEQN") %>% 
  model.matrix(rob_mod7, data = .)

# Fit model
rob7 <- Log.lqr(rob7_matrix, y = nh_grps$prop_CAL_sites3mm, p = 0.5)

# Extract coefficients 
rob7_coef <- bind_cols(Term = colnames(rob7_matrix),
                       as_tibble(rob7$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)


# Median predictions for TC1 from robust regression
TC1_medians_rob7 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob7$fitted.values)) %>% 
  group_by(TC_decile = TC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted median: robust regression (decile)" = median(predicted_robust)) %>% 
  gather(Measurement, value, -TC_decile)

TC7_medians_rob7 <- nh_grps %>% 
  cbind(predicted_robust = invlogit(rob7$fitted.values)) %>% 
  group_by(TC_decile = TC7_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted median: robust regression (decile)" = median(predicted_robust)) %>% 
  gather(Measurement, value, -TC_decile)


## Robust regression with three category smoking status and linear functions of TCs ##

rob_mod8 <- ~ 1 + RIDAGEYR + RIAGENDR + 
  smoking + diabetes + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + 
  TC7 + TC8

# Model the median
rob8 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                x = model.matrix(rob_mod8, data = nh_grps),
                p = 0.5)

rob8_coef <- bind_cols(Term = colnames(model.matrix(rob_mod8, data = nh_grps)),
                       as_tibble(rob8$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)

## Robust regression with three category smoking status and quadratic functions of TCs ##

rob_mod9 <- ~ 1 + RIDAGEYR + RIAGENDR + 
  smoking + diabetes + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + 
  TC7 + TC8 + I(TC1^2) + I(TC2^2) + I(TC3^2) + I(TC4^2) + I(TC5^2) + I(TC6^2) + I(TC7^2) + I(TC8^2)

# Model the median
rob9 <- Log.lqr(y = nh_grps$prop_CAL_sites3mm, 
                x = model.matrix(rob_mod9, data = nh_grps),
                p = 0.5)

rob9_coef <- bind_cols(Term = colnames(model.matrix(rob_mod9, data = nh_grps)),
                       as_tibble(rob9$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)) %>% 
  rename(Sig = V1)


# Compare the AIC of the three models with three category smoking status
rob5$AIC
rob8$AIC
rob9$AIC


### Analysis of tooth count ###


# Model the median
tooth1 <- Log.lqr(y = nh_grps$tooth_count, 
                  a =1, b = 32,
                x = model.matrix(rob_mod, data = nh_grps),
                p = 0.5)

tooth1_coef <- bind_cols(Term = colnames(model.matrix(rob_mod, data = nh_grps)),
                       as_tibble(tooth1$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)#,
         # Why are these values so high?
         #"Predicted proportion" = round(invlogit(ifelse(Term == "(Intercept)", Estimate[1], Estimate[1] + Estimate)), 2)
  ) %>% 
  rename(Sig = V1)



# Truncate tooth count at 28
nh_grps$tooth_count_trunc <- if_else(nh_grps$tooth_count > 28, as.integer(28), nh_grps$tooth_count)

# Model the median
tooth2 <- Log.lqr(y = nh_grps$tooth_count_trunc, 
                  a =1, b = 28,
                  x = model.matrix(rob_mod, data = nh_grps),
                  p = 0.5)

tooth2_coef <- bind_cols(Term = colnames(model.matrix(rob_mod, data = nh_grps)),
                         as_tibble(tooth1$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)#,
         # Why are these values so high?
         #"Predicted proportion" = round(invlogit(ifelse(Term == "(Intercept)", Estimate[1], Estimate[1] + Estimate)), 2)
  ) %>% 
  rename(Sig = V1)

