# Analysis of diet and periodontal disease measures from NHANES
# Joining of NHANES datasets
# Using average of two 24hour recall periods
require(here)


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
  select(SEQN, SDDSRVYR, RIDSTATR, RIAGENDR, RIDAGEYR)


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



# # Lookup table for non-dietary variables
# names(nhanes)
# library(nhanesA)
# 
# nhanesTableVars()
# #DEMO_F, SMQ_F,
# nm2 <- lapply(list("DEMO_F", "SMQ_F", "OHXPER_F"), nhanesTranslate, colnames = names(nhanes)[str_detect(names(nhanes), "^[A-Z]+$")])
# nhanesTranslate("SMQ_F", colnames = names(smoking))
# nm3 <- nhanesTranslate("DEMO_F", colnames = names(demographic))
# nhanesTranslate("DEMO_F", "SDDSRVYR")
# demographic2 <- 
# demographic3 <- nhanesTranslate("DEMO_F", "RIDSTATR")
# demographic4 <- nhanesTranslate("DEMO_F", colnames = c("SEQN", "SDDSRVYR", "RIDSTATR", "RIAGENDR", "RIDAGEYR"), data = nhanes("DEMO_F"))
# prp <- nhanesTranslate("OHXPER_F", names)



#rm(demographic, smoking, tooth_count, energy, perio_raw, perio_cal, perio_pocket, diabetes)

### Treelet analysis ###


# Stability of components (sign test)


# ### Models ###
#   
#   invlogit <- 
#     function (x) 
#     {
#       1/(1 + exp(-x))
#     }  
#   
# 
# 
# ### Clinical attachment loss ###
# 
# exp(coef(cal3)) %>% enframe() %>% print(n = Inf)
# invlogit(coef(cal3)) %>% enframe() %>% print(n = Inf)
# nh %>% count(sites_assessed)
# 
# TC3 -ve
# TC4 ++ve
# TC5 +ve
# 
# # # Random forests
# library(randomForest)
# nhrf <- select(nh, c("CAL_sites", adj_vars, diet_vars))
# train <- sample(1:nrow(nhrf), round(nrow(nhrf)*0.8))
# # Slow
# cal0_rf <- randomForest(formula = CAL_sites ~ ., data = nhrf, subset = train)



