# Analysis of food group intake and periodontal disease measures from NHANES

require(tidyverse)
require(treelet) # For treelet transform
require(betareg) # For Beta regression
require(Formula) # For multipart formulas used in Beta regression
require(lmtest) # For likelihood ratio tests of Beta regressions
require(quantreg) # For quantile regression
require(lqr) # For logistic quantile regression
require(cdfquantreg) # For logistic quantile regression (more recent package)
require(ggdendro)
require(splines)
require(broom)

cat("\n Treelet analysis on food groups\n")

# Define food groups
food_groups_nhanes <- food_grps_grms_per_day %>% 
  
  select(SEQN, 
         matches("BEV(0|21|22|231|232|233|241|242)$"),
         EGG0,
         matches("FAT(1|2)$"),
         matches("FRUIT(11|2|31|32|33|34|35)$"),
         matches("GRAIN(1|21|22|23|3|4|5|6)$"),
         matches("MEAT(1|2|3|4|5|6|7|8)$"),
         matches("MILK(11|2|3|4)$"),
         matches("SUGAR(1|2)$"),
         matches("VEG(1|2|3|4|5|6|7|8)$"),
         WATER1)

# List of groups of beverages included
bevs_included <- names(select(food_groups_nhanes, matches("BEV|WATER")))

# Scale by overall energy intake
food_groups_scl <- food_groups_nhanes %>% 
  mutate_at(vars(-SEQN), ~./nhanes$KCAL) %>% 
  select(-SEQN) %>%
  scale()

food_groups_cor <- cor(food_groups_scl)

# Generate maximum height tree
# Save basis vectors at all cut levels so the most useful can be selected later
food_groups_tree_full <- Run_JTree(food_groups_cor, maxlev = ncol(food_groups_cor)-1, whichsave = 1:(ncol(food_groups_cor)-1))
# Extract the treelet components and associated variances
food_groups_tc_full <- TTVariances(food_groups_tree_full, food_groups_cor)

# Dendrogram of maximum height tree
food_groups_dendro <- dendro_data(ConvertTTDendro(food_groups_tree_full))


# Initial analysis - data driven selection of number of components and cut points

## Cross validation ##
# Data driven cross validation to choose a cut level that can describe the 
# data with few components (TCs)
# Set the desired number of components (m) based on scree plot and then find the best cut level
# Cross validation scores for each cut level
m_grps <- 8
cvs_grps <- replicate(5, CrossValidateTT(food_groups_scl, m = m_grps))


# Fit reduced treelet
# Project selected components to get scores for each individual using the selected cut level
# Cut level selected based on  cross validation and inspection of dendrogram
# Original
food_groups_tc_reduced <- TTVariances(food_groups_tree_full, cor(food_groups_scl), cut_level = 28, components = m_grps)
food_groups_tc_scores <- food_groups_scl %*% food_groups_tc_reduced$tc

cat("\n Calculate correlation of TC scores ")
food_groups_tc_cor <- cor(food_groups_tc_scores, method = "kendall") %>% 
  formatC(digits = 2, format = "f")
food_groups_tc_cor[upper.tri(food_groups_tc_cor, diag = T)] <- ""


### Add TC scores to main NHANES dataset
# The resulting dataset is distinct from the one used in the nutrient analysis ('nh')
nh_grps <- bind_cols(nhanes, as_tibble(food_groups_tc_scores)) %>% 
  # Classify TC scores into deciles
  mutate_at(vars(matches("^TC[0-9]{1,2}$")), list(~as.factor(ntile(., n = 10)))) %>% 
  rename_at(vars(matches("^TC[0-9]{1,2}$")), ~paste0(., "_dec")) %>% 
  # Raw TC scores
  bind_cols(as_tibble(food_groups_tc_scores)) %>% 
  inner_join(dietary, by = "SEQN") %>% 
  
  # Standardise age and KCAL to make predictions easier
  mutate(RIDAGEYR = as.numeric(scale(RIDAGEYR)),
         KCAL = as.numeric(scale(KCAL)))

#food_groups_out <- bind_cols(diet, data.frame(food_groups_tc_scores))

# Extract loadings for TCs
food_groups_loadings <- food_groups_tc_reduced$tc %>% 
  as_tibble(rownames = "Variable") %>% 
  gather(Component, Value, -Variable) %>% 
  # Order by leaf labelling order
  mutate(Variable = factor(Variable, levels = as.character(food_groups_dendro$labels$label)))


cat("\n Beta regression of CAL by food groups")

nh_grps$prop_CAL_sites3mm_beta <- PropTransform(nh_grps$prop_CAL_sites3mm)

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

age45 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 45])) %>% 
  as.numeric() %*% rob1$beta %>% 
  invlogit()*100

age59 <- rob_pred_frame %>% 
  mutate(RIDAGEYR = unique(nh_grps$RIDAGEYR[nhanes$RIDAGEYR == 59])) %>% 
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
cor_tc1 <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC1, ., method = "kendall"), 3)) 


# Select variables correlated with TC7
cor_tc7 <- nh_grps %>% 
  select(one_of(pull(diet_names, `Variable Name`))) %>% 
  map_dbl(~round(cor(nh_grps$TC7, ., method = "kendall"), 3))

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
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)#,
         # Why are these values so high?
         #"Predicted proportion" = round(invlogit(ifelse(Term == "(Intercept)", Estimate[1], Estimate[1] + Estimate)), 2)
  ) %>% 
  rename(Sig = V1)

cat("\n Treelet analysis on food groups completed")


