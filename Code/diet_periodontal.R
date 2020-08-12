# Analysis of diet and periodontal disease measures from NHANES
# Joining of NHANES datasets
# Using average of two 24hour recall periods
library(here)
library(tidyverse)

source(here("Code", "prepare_periodontal_data.R"))
source(here("Code", "prepare_dietary_data.R"))
source(here("Code", "treelet_functions.R"))

cat("\n Preparing analysis cohort")

## Smoking status questionaire ##
smoking <- list.files(path = here("NHANES"), pattern = "SMQ", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  transmute(SEQN, SMQ020, SMQ040,
            smoking = factor(case_when(SMQ020 %in% c(2, 7, 9) ~ "Never",
                                SMQ020 == 1 & SMQ040 %in% c(3, 7, 9) ~ "Former",
                                SMQ020 == 1 & is.na(SMQ040) ~ "Former",
                                SMQ020 == 1 & SMQ040 %in% c(1,2) ~ "Current"), 
            levels = c("Never", "Former", "Current")))

## Demographic data ##

demographic <- list.files(path = here("NHANES"), pattern = "DEMO", recursive = T, full.names = T) %>% 
  map(~ cbind(read_xpt(.), 
              # Including wave membership
              WAVE = str_extract(., "[A-Z]{1}(?=.XPT)"),
              stringsAsFactors = FALSE)) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(WAVE = if_else(WAVE == "O", "A", WAVE)) %>% 
  # Set refused or don't know for edcuational attainment to missing
  mutate_at(vars(DMDEDUC2), ~if_else(. %in% c(7, 9), as.numeric(NA), .)) %>% 
  mutate_at(vars(DMDEDUC2, RIDRETH1, INDFMIN2), ~fct_explicit_na(as.factor(.))) %>% 
  select(SEQN, WAVE, SDDSRVYR, RIDSTATR, RIAGENDR, RIDAGEYR, RIDRETH1, INDFMIN2, DMDEDUC2)

demographic %>% count(DMDEDUC2)

## Diabetes data ##
# HbA1c
hba1c <- list.files(path = here("NHANES"), pattern = "GHB", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  select(SEQN, LBXGH)
# Questionnaire
diabetes <- list.files(path = here("NHANES"), pattern = "DIQ", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  select(SEQN, DIQ010) %>% 
  full_join(hba1c, by = "SEQN") %>% 
  # Composite measure using HbA1c or questionnaire
  mutate(diabetes = if_else((!is.na(LBXGH) & LBXGH > 6.5 )|(!is.na(DIQ010) & DIQ010 == 1), T, F))
  


### Construct the cohort ###
# Join all datasets except dietary #

nhanes_all <- demographic %>% 
  left_join(smoking, by = "SEQN") %>% 
  left_join(tooth_count, by = "SEQN") %>% 
  left_join(energy, by = "SEQN") %>% 
  left_join(select(perio_raw, SEQN, EXCLU, PDSTS), by = "SEQN") %>% 
  left_join(perio_cal, by = "SEQN") %>% 
  left_join(perio_pocket, by = "SEQN") %>% 
  left_join(diabetes, by = "SEQN") %>% 
  left_join(transmute(dietary, SEQN, dietary = T), by = "SEQN")

    
## Apply exclusion criteria
nhanes <- nhanes_all %>% 
  # Include only the selected waves
  filter(between(SDDSRVYR, 6, 8)) %>% 
  # Exclude aged < 30
  filter(RIDAGEYR >= 30) %>% 
  # Exclude with no medical exam
  filter(RIDSTATR == 2) %>% 
  # Exclude edentulous people
  filter(tooth_count > 0)  %>% 
  # Exclude people with medical exclusion from periodontal exam
  filter(EXCLU != 1) %>% 
  # Exclude people with no completed periodontal exam
  filter(PDSTS == 1) %>% 
  #7 missing from these exclusions compared with previous paper
  # Remove those where a dental exam was not carried out
  #  filter(OHDEXSTS != 3) #%>% 
    
  # Exclude people with no sites assessed for CAL
  filter(CAL_sites_assessed > 0) %>% 
  # Exclude people with no sites assessed for pocket depth
  #filter(PD_sites_assessed > 0) %>% 
  
  # Exclude people with no in person dietary questionnaire
  # Check that this is not flagged elsewhere 
  filter(dietary == T) %>% 

  # Exclude people with no diabetes flag (none)
  filter(!is.na(diabetes)) %>% 
  
  # Classify into quartiles for descriptive tables
  mutate(CAL_sites_quartile = ntile(prop_CAL_sites3mm, 4))

# Extract dietary data for just those in the cohort
dietary <- dietary %>% 
  inner_join(select(nhanes, SEQN), by = "SEQN") 

food_grps_grms_per_day <- food_grps_grms_per_day %>% 
  # Only records for cohort members
  inner_join(select(nhanes, SEQN), by = "SEQN")
  
food_total_grms_per_day <- food_total_grms_per_day %>% 
  # Only records for cohort members
  inner_join(select(nhanes, SEQN), by = "SEQN")

food_grps_per_day <- food_grps_per_day %>% 
  inner_join(select(nhanes, SEQN), by = "SEQN") 

# food_grps_carbs_per_day <- food_grps_carbs_per_day %>% 
#   inner_join(select(nhanes, SEQN), by = "SEQN")
# 
# food_grps_sugars_per_day <- food_grps_sugars_per_day %>% 
#   inner_join(select(nhanes, SEQN), by = "SEQN")

# Extract site level periodontal measurements for just those in the cohort
perio_cal_site <- perio_raw %>% 
  select(SEQN, matches("^[0-9]{2}LA")) %>% 
  gather(variable, value, -SEQN) %>% 
  mutate(tooth = str_extract(variable, "^[0-9]{2}"),
         site = str_extract(variable, "[A-Z]{3}$")) %>% 
  inner_join(select(nhanes, SEQN, RIDAGEYR), by = "SEQN") %>% 
  filter(!is.na(value))

perio_pocket_site <- perio_raw %>% 
  select(SEQN, matches("^[0-9]{2}PC")) %>% 
  gather(variable, value, -SEQN) %>% 
  mutate(tooth = str_extract(variable, "^[0-9]{2}"),
         site = str_extract(variable, "[A-Z]{3}$")) %>% 
  inner_join(select(nhanes, SEQN, RIDAGEYR), by = "SEQN") %>% 
  filter(!is.na(value))

cat("\n Cohort prepared")

# Comparison of total food weight and total energy intake 
food_energy_weight <- food_total_grms_per_day %>% 
  inner_join(nhanes, by = "SEQN") %>% 
  select(KCAL, GRMS)

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

# Scale by overall energy intake (the nutrient density method)
food_groups_scl <- food_groups_nhanes %>% 
  mutate_at(vars(-SEQN), ~./nhanes$KCAL) %>% 
  select(-SEQN) %>%
  scale()

food_groups_cor <- cor(food_groups_scl)

# Scale by overall energy intake (the residual method)
food_groups_cor_resid <- food_groups_nhanes %>% 
  select(-SEQN) %>% 
  modify(~resid(lm(.~ nhanes$KCAL, data = food_groups_nhanes))) %>% cor()
#food_groups_cor <- food_groups_cor_resid


cat("\n Food groups prepared")
