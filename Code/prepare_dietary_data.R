# Analysis of diet and periodontal disease measures from NHANES
# Preparation of NHANES dietary data
# Averaging over two 24 hour recall periods


### Dietary data ###

# Day 1 recall questionnaire was in person
# Day 2 was by phone



## Total nutrient intake files ##

# Read in data from both days for each person
# 1.5 mins
dietary <- list.files(path = here("NHANES"), pattern = "DR(1|2)TOT_[A-Z]{1}.XPT", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  select(SEQN, matches("^DR(1|2)T"), -matches("NUMF|WS$")) %>% 
  # Take average values across days, use existing value if only one present
  gather(key, value, -SEQN) %>% 
  mutate(variable = str_replace(key, "DR[0-9]{1}T", "")) %>% 
  group_by(SEQN, variable) %>% 
  summarise(value = mean(value, na.rm = T)) %>% 
  mutate(value = if_else(is.nan(value), as.double(NA), value)) %>% 
  spread(key = variable, value = value) %>% 
  ungroup() %>% 
  # Drop all rows with missing values
  na.omit()
  
# Variable names
# Note that in one wave (H) DR1TNUMF covers both foods and beverages
diet_names <- list.files(path = here("NHANES"), pattern = "DR1TOT_[A-Z]{1}_variable_names.csv", recursive = T, full.names = T) %>% 
  lapply(read_csv) %>% 
  bind_rows() %>% 
  distinct() %>% 
  mutate(`Variable Name` = str_replace(`Day1 Name`, "DR[0-9]{1}T", "")) %>% 
  transmute(`Variable Name`, Description = `Variable Label`) %>% 
  filter(`Variable Name` %in% setdiff(names(dietary), "SEQN")) %>% 
  filter(Description %in% c("Total folate (mcg)", "Number of foods/beverages reported") == F)


# Remove energy intake from dietary dataset and join separately
energy <- dietary %>% 
  select(SEQN, KCAL)

dietary <- dietary %>% 
  select(-KCAL)



