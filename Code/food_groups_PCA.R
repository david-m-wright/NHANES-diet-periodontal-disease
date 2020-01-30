# Principal components analysis on food groups

library(tidyverse)
library(lqr) # For logistic quantile regression

# Peform PCA
# Note that because the PCA is on the scaled food group dataset this is equivalent to using the correlation matrix
pca_food_groups <- prcomp(food_groups_scl, center = FALSE, scale. = FALSE)

# Reverse the signs on the eigenvectors (loadings) and scores so PC loadings are directly comparable with TC loadings
pca_food_groups$rotation <- -pca_food_groups$rotation
pca_food_groups$x <- -pca_food_groups$x



# PC scores
food_groups_pc_scores <- as_tibble(pca_food_groups$x) %>% 
  rename_all(~str_replace(., "Comp.([0-9]{1,2})", "PC\\1")) %>% 
  select(1:8)

## Add PC scores for first 8 components to main NHANES dataset
# The resulting dataset is distinct from the one used in the nutrient analysis ('nh')
nh_grps_pca <- bind_cols(nhanes, food_groups_pc_scores) %>% 
                         
  # Classify PC scores into deciles
  mutate_at(vars(matches("^PC[0-9]{1,2}$")), list(~as.factor(ntile(., n = 10)))) %>% 
  rename_at(vars(matches("^PC[0-9]{1,2}$")), ~paste0(., "_dec")) %>% 
  # Raw TC scores
  bind_cols(food_groups_pc_scores) %>% 
  inner_join(dietary, by = "SEQN") %>% 
  
  # Standardise age and KCAL to make predictions easier
  mutate(RIDAGEYR = as.numeric(scale(RIDAGEYR)),
         KCAL = as.numeric(scale(KCAL)))


# Extract loadings for PCs
food_groups_pc_loadings <- pca_food_groups$rotation %>% 
  as.data.frame() %>% 
  as_tibble(rownames="Variable") %>% 
  rename_all(~str_replace(., "Comp.([0-9]{1,2})", "PC\\1")) %>% 
  select(Variable, PC1:PC8) %>% 
  gather(Component, Value, -Variable) %>% 
  # Order by leaf labelling order
  mutate(Variable = factor(Variable, levels = food_groups_dendro_order))




cat("\n Robust quantile regression by food groups")

# Specify model for the median 
rob_mod_pca <- ~ 1 + RIDAGEYR + RIAGENDR + 
  SMQ020 + diabetes + PC1_dec + PC2_dec + PC3_dec + PC4_dec + PC5_dec + PC6_dec + 
  PC7_dec + PC8_dec

# Model the median
rob1_pca <- Log.lqr(y = nh_grps_pca$prop_CAL_sites3mm, 
                x = model.matrix(rob_mod_pca, data = nh_grps_pca),
                p = 0.5)

# Model the lower quartile
rob1_lower_pca <- Log.lqr(y = nh_grps_pca$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod_pca, data = nh_grps_pca),
                      p = 0.25) 

# Model the upper quartile
rob1_upper_pca <- Log.lqr(y = nh_grps_pca$prop_CAL_sites3mm, 
                      x = model.matrix(rob_mod_pca, data = nh_grps_pca),
                      p = 0.75) 

#model.matrix(rob_mod, data = nh)[2,] %*% rob1$beta
#invlogit(rob1$fitted.values[2])
#invlogit(as.numeric(rob_pred_frame) %*% rob1$beta)

# Predict at the mean age, male, non-smoker, non-diabetic, all TCs at decile 5
rob_pred_frame_pca <- model.matrix(rob_mod_pca, data = nh_grps_pca)[1,] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(RIDAGEYR = 0,
         RIAGENDR = 1,
         SMQ020 = 2, 
         diabetesTRUE = 0) %>% 
  mutate_at(vars(matches("_dec")), function(.){0}) %>% 
  mutate_at(vars(matches("_dec5")), function(.){1}) 

# Extract coefficients for median model
rob1_coef_pca <- bind_cols(Term = colnames(model.matrix(rob_mod_pca, data = nh_grps_pca)),
                       as_tibble(rob1_pca$table, .name_repair = "minimal")) %>% 
  mutate(OR = formatC(exp(Estimate), format = "f", digits = 2)#,
         # Why are these values so high?
         #"Predicted proportion" = round(invlogit(ifelse(Term == "(Intercept)", Estimate[1], Estimate[1] + Estimate)), 2)
  ) %>% 
  rename(Sig = V1)


# Median predictions for PC1 from robust regression
PC1_medians <- nh_grps_pca %>% 
  cbind(predicted_robust = invlogit(rob1_pca$fitted.values),
        predicted_upper = invlogit(rob1_upper_pca$fitted.values),
        predicted_lower = invlogit(rob1_lower_pca$fitted.values)) %>% 
  group_by(PC_decile = PC1_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -PC_decile)

# Median predictions for PC1 from robust regression
PC6_medians <- nh_grps_pca %>% 
  cbind(predicted_robust = invlogit(rob1_pca$fitted.values),
        predicted_upper = invlogit(rob1_upper_pca$fitted.values),
        predicted_lower = invlogit(rob1_lower_pca$fitted.values)) %>% 
  group_by(PC_decile = PC6_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -PC_decile)

# Median predictions for PC1 from robust regression
PC8_medians <- nh_grps_pca %>% 
  cbind(predicted_robust = invlogit(rob1_pca$fitted.values),
        predicted_upper = invlogit(rob1_upper_pca$fitted.values),
        predicted_lower = invlogit(rob1_lower_pca$fitted.values)) %>% 
  group_by(PC_decile = PC8_dec) %>% 
  summarise("Actual median" = median(prop_CAL_sites3mm),
            "Predicted 25th percentile" = median(predicted_lower),
            "Predicted median: robust regression (decile)" = median(predicted_robust),
            "Predicted 75th percentile" = median(predicted_upper)) %>% 
  gather(Measurement, value, -PC_decile)


