# Analysis of nutrient intake and periodontal disease measures from NHANES
require(tidyverse)
require(treelet) # For treelet transform
require(betareg) # For Beta regression
require(Formula) # For multipart formulas used in Beta regression
require(lmtest) # For likelihood ratio tests of Beta regressions
require(ggdendro)
require(splines)
require(broom)

cat("\n Running treelet analysis on nutrients\n")

# Scale data so can intepret treelets on correlation scale
diet_scl <- dietary %>% 
  select(-SEQN) %>% 
  scale()

# Generate maximum height tree
# Save basis vectors at all cut levels so the most useful can be selected later
diet_tree_full <- Run_JTree(cor(diet_scl), maxlev = ncol(diet_scl)-1, whichsave = 1:(ncol(diet_scl)-1))
# Extract the treelet components and associated variances
diet_tc_full <- TTVariances(diet_tree_full, cor(diet_scl))

# Dendrogram of maximum height tree
diet_dendro <- dendro_data(ConvertTTDendro(diet_tree_full))


## Initial analysis - data driven selection of number of components and cut points

## Cross validation ##
# Data driven cross validation to choose a cut level that can describe the 
# data with few components (TCs)
# Set the desired number of components (m) and then find the best cut level
# Cross validation scores for each cut level
m = 7
cvs <- replicate(5, CrossValidateTT(diet_scl, m = m))

# Project selected components to get scores for each individual using the selected cut level
# Original
diet_tc_reduced <- TTVariances(diet_tree_full, cor(diet_scl), cut_level = 50, components = m)

diet_tc_scores <- diet_scl %*% diet_tc_reduced$tc

cat("\n Calculating TC correlations")
diet_tc_cor <- cor(diet_tc_scores) %>% 
  formatC(digits = 2, format = "f")
diet_tc_cor[upper.tri(diet_tc_cor, diag = T)] <- ""


### Add TC scores to main NHANES dataset
nh <- bind_cols(nhanes, as_tibble(diet_tc_scores)) %>% 
  # Classify TC scores into deciles
  mutate_at(vars(matches("^TC[0-9]{1,2}$")), list(~as.factor(ntile(., n = 10)))) %>% 
  rename_at(vars(matches("^TC[0-9]{1,2}$")), ~paste0(., "_dec")) %>% 
  # Raw TC scores
  bind_cols(as_tibble(diet_tc_scores)) %>% 
  inner_join(dietary, by = "SEQN")

# Extract loadings for TCs
diet_loadings <- diet_tc_reduced$tc %>% 
  as_tibble(rownames = "Variable") %>% 
  gather(Component, Value, -Variable) %>% 
  # Order by leaf labelling order
  mutate(Variable = factor(Variable, levels = as.character(diet_dendro$labels$label)))


cat("\n Beta regression analysis of CAL")

# Variables to adjust for
adj_vars <- c("RIDAGEYR", "RIAGENDR", "SMQ020", "diabetes", "KCAL")

# Transform the proportion of CAL sites so that all values lie within the interval (0,1) for beta regression
nh$prop_CAL_sites3mm_beta <- PropTransform(nh$prop_CAL_sites3mm)

# Beta regression model with fixed dispersion
cal1_form <- as.Formula(paste("prop_CAL_sites3mm_beta ~", paste(adj_vars, collapse = "+")))
cal1 <- betareg(cal1_form, data = nh)
# Beta regression with variable dispersion model
cal2_form <- as.Formula(paste("prop_CAL_sites3mm_beta ~", paste(adj_vars, collapse = "+"), "|", paste(adj_vars, collapse = "+")))
cal2 <- betareg(cal2_form, data = nh)
# Dropping KCAL from the precision model (not significant)
cal3 <- update(cal2, .~. |.-KCAL, data = nh)

# Adding TCs to the precision model
cal4 <- update(cal3, .~. + TC1 + TC2 + TC3 + TC4 + TC5 + TC6 + TC7 |., data = nh)

# CAL TC polynomial models
# Cubic spline fit
cal5 <- update(cal4, .~. - TC1 - TC2 - TC3 - TC4 - TC5 - TC6 - TC7 + bs(TC1) + bs(TC2) + bs(TC3) + bs(TC4) + bs(TC5) + bs(TC6) + bs(TC7)|., data = nh)

# Quadratic fit
cal6 <- update(cal4, .~. + I(TC1^2) + I(TC2^2) + I(TC3^2) + I(TC4^2) + I(TC5^2) + I(TC6^2) + I(TC7^2) |., data = nh)

# Note that bias correction or bias reduction made very little difference to estimates
cal4bc <- update(cal4, type = "BC")
cal4br <- update(cal4, type = "BR")
#format(exp(cbind(coef(cal4), coef(cal4bc), coef(cal4br))), digits = 4)

# Model to predict from
mpred <- cal4

# Take TC quantiles and stack
tc_qtx <- StackQuantiles(diet_tc_scores, probs = seq(0.02, 0.98, by = 0.01))

# Prediction for average aged, male, non-smoker, non-diabetic, average energy intake
cal_pred_frame <- cbind(RIDAGEYR = mean(nh$RIDAGEYR), 
                        RIAGENDR = 1, 
                        SMQ020 = 2, 
                        diabetes = F,
                        KCAL = mean(nh$KCAL),
                        tc_qtx) %>% 
  as_tibble()


## Varying the cut level of the cluster tree ##

# Second treelet analysis, varying the number of components and cut levels
diet_tc_reduced2 <- TTVariances(diet_tree_full, cor(diet_scl), cut_level = 40, components = 9)
diet_tc_scores2 <- diet_scl %*% diet_tc_reduced2$tc

# Correlation of TC scores
diet_tc_cor2 <- cor(diet_tc_scores2) %>% 
  formatC(digits = 2, format = "f")
diet_tc_cor2[upper.tri(diet_tc_cor2, diag = T)] <- ""

# Cross validation
cvs2 <- replicate(5, CrossValidateTT(diet_scl, m = 9))

nh2 <- bind_cols(nhanes, as_tibble(diet_tc_scores2)) %>% 
  inner_join(dietary, by = "SEQN") %>% 
  mutate(prop_CAL_sites3mm_beta = PropTransform(prop_CAL_sites3mm))

# Extract loadings for TCs
diet_loadings2 <- diet_tc_reduced2$tc %>% 
  as_tibble(rownames = "Variable") %>% 
  gather(Component, Value, -Variable) %>% 
  # Order by leaf labelling order
  mutate(Variable = factor(Variable, levels = as.character(diet_dendro$labels$label)))


# Second analysis CAL TC models
cal7 <- update(cal4, .~. + TC8 + TC9 |., data = nh2)
cal8 <- update(cal6, .~. + TC8 + TC9 + I(TC8^2) + I(TC9^2)|., data = nh2)

# Prediction for average aged, male, non-smoker, non-diabetic, average energy intake
cal_pred_frame2 <- cbind(RIDAGEYR = mean(nh$RIDAGEYR), 
                         RIAGENDR = 1, 
                         SMQ020 = 2, 
                         diabetes = F,
                         KCAL = mean(nh$KCAL),
                         StackQuantiles(diet_tc_scores2, probs = seq(0.02, 0.98, by = 0.01))) %>% 
  as_tibble()

cat("\n Nutrient analysis completed")
