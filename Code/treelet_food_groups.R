# Treelet analysis of food group data
library(tidyverse)
library(treelet) # For treelet transform
library(ggdendro)
library(here)

source(here("Code", "treelet_functions.R"))

cat("\n Treelet analysis on food groups\n")

# Generate maximum height tree
# Save basis vectors at all cut levels so the most useful can be selected later
food_groups_tree_full <- Run_JTree(food_groups_cor, maxlev = ncol(food_groups_cor)-1, whichsave = 1:(ncol(food_groups_cor)-1))
# Extract the treelet components and associated variances
food_groups_tc_full <- TTVariances(food_groups_tree_full, food_groups_cor)

# Dendrogram of maximum height tree
food_groups_dendro <- dendro_data(ConvertTTDendro(food_groups_tree_full))
# Plotting order
food_groups_dendro_order <- as.character(food_groups_dendro$labels$label)

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
         KCAL_raw = KCAL,
         KCAL = as.numeric(scale(KCAL))) %>% 
  
  # Set up outcome variable for Beta regression
  mutate(prop_CAL_sites3mm_beta = PropTransform(prop_CAL_sites3mm)) 

#food_groups_out <- bind_cols(diet, data.frame(food_groups_tc_scores))

WrapLabel <- function(x){
  xwrap <- str_wrap(x, width = 30)
  # str_pad(
  if_else(str_detect(xwrap, "\n"), 
          as.character(xwrap), 
          paste0("\n", as.character(xwrap))) %>% 
    as_factor()
  # width = 30, side = "right")
}

# Extract loadings for TCs
food_groups_loadings <- food_groups_tc_reduced$tc %>% 
  as_tibble(rownames = "Variable") %>% 
  gather(Component, Value, -Variable) %>% 
  # Add food group descriptions
  inner_join(fgrp %>% 
               distinct(grp_code, grp_description),
             by = c("Variable" = "grp_code")) %>% 
  # Order by leaf labelling order
  mutate(Variable = factor(Variable, levels = food_groups_dendro_order)) %>% 
    arrange(Variable) %>% 
    mutate(grp_description = as_factor(grp_description),
           grp_padded = WrapLabel(grp_description)) 

