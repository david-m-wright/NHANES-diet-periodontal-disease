# Analysis of diet and periodontal disease measures from NHANES
# Joining of NHANES datasets
# Using average of two 24hour recall periods
require(here)
require(tidyverse)

source(here("Code", "prepare_periodontal_data.R"))
source(here("Code", "prepare_dietary_data.R"))
source(here("Code", "treelet_functions.R"))

## Smoking status questionaire ##
smoking <- list.files(path = here("NHANES"), pattern = "SMQ", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  select(SEQN, SMQ020)

## Demographic data ##
demographic <- list.files(path = here("NHANES"), pattern = "DEMO", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  select(SEQN, SDDSRVYR, RIDSTATR, RIAGENDR, RIDAGEYR, INDFMPIR, DMDEDUC2)


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

nhanes <- demographic %>% 
  left_join(smoking, by = "SEQN") %>% 
  left_join(tooth_count, by = "SEQN") %>% 
  left_join(energy, by = "SEQN") %>% 
  left_join(select(perio_raw, SEQN, EXCLU, PDSTS), by = "SEQN") %>% 
  left_join(perio_cal, by = "SEQN") %>% 
  left_join(perio_pocket, by = "SEQN") %>% 
  left_join(diabetes, by = "SEQN") %>% 
  left_join(transmute(dietary, SEQN, dietary = T), by = "SEQN") %>% 
  
## Exclusion criteria
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
  filter(PD_sites_assessed > 0) %>% 
  
  # Exclude people with no in person dietary questionnaire
  # Check that this is not flagged elsewhere 
  filter(dietary == T) %>% 

  # Exclude people with no diabetes flag (none)
  filter(!is.na(diabetes))

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

cat("\n Cohort prepared")