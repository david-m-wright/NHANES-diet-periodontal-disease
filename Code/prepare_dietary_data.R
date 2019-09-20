# Analysis of diet and periodontal disease measures from NHANES
# Preparation of NHANES dietary data
# Averaging over two 24 hour recall periods

require(here)
require(tidyverse)
require(haven)
require(fuzzyjoin)


### Helper function ###

# Function to check whether x falls within an interval defined in y
# Args: x = numeric value
# y = interval specified in string format "start end" where start <= end
# Value: logical indicating whether x falls within y
CheckWithinInterval <- function(x, y){
  tibble(xcol = x,
         interval_start = as.numeric(str_extract(y, ".+(?= )")),
         interval_end = as.numeric(str_extract(y, "(?<= ).+")),
         isWithin = if_else(xcol >= interval_start & xcol <= interval_end, T, F)) %>% 
    select(isWithin)
}
# CheckWithinInterval(2, "1 10")
# CheckWithinInterval(0, "1 10")


cat("\n Preparing dietary data")

### Dietary data ###

# Day 1 recall questionnaire was in person
# Day 2 was by phone


## Total nutrient intake files ##

# Read in data from both days for each person

# 1.5 mins
dietary <-  list.files(path = here("NHANES"), pattern = "DR(1|2)TOT_[A-Z]{1}.XPT", recursive = T, full.names = T) %>% 
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



## Individual food files ##
# iff1 <- read_xpt(here("NHANES/2009_2010/DR1IFF_F.XPT"))

# First day intake
food1 <- list.files(path = here("NHANES"), pattern = "DR1IFF_[A-Z]{1}.XPT", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  transmute(SEQN, FDCD = DR1IFDCD, 
            DR1ISUGR, 
            DR1ICARB,
            DR1IGRMS, DRDINT) %>% 
  # Remove rows with no grams (human milk)
  filter(!is.na(DR1IGRMS)) %>% 
  # Sum total intake on day 1
  group_by(SEQN, FDCD, DRDINT) %>% 
  summarise(DR1IGRMS = sum(DR1IGRMS),
            DR1ICARB = sum(DR1ICARB),
            DR1ISUGR = sum(DR1ISUGR))


# Second day intake
food2 <- list.files(path = here("NHANES"), pattern = "DR2IFF_[A-Z]{1}.XPT", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows() %>% 
  transmute(SEQN, FDCD = DR2IFDCD, 
            DR2ISUGR,
            DR2ICARB,
            DR2IGRMS, DRDINT) %>% 
  # Remove rows with no grams (human milk)
  filter(!is.na(DR2IGRMS)) %>% 
  # Sum total intake on day 2
  group_by(SEQN, FDCD, DRDINT) %>% 
  summarise(DR2IGRMS = sum(DR2IGRMS),
            DR2ISUGR = sum(DR2ISUGR),
            DR2ICARB = sum(DR2ICARB))


# # Average intake across days 
# foods_grms_per_day <- food1 %>% 
#   full_join(food2, by = c("SEQN", "FDCD", "DRDINT")) %>% 
#   # Take the mean if an individual was measured on both days, 
#   mutate_at(c("DR1IGRMS", "DR2IGRMS"), list(function(.){if_else(is.na(.), 0, .)})) %>% 
#   # Otherwise take the exising value (the parallel maximum here)
#   mutate(GRMS = if_else(DRDINT == 2, (DR1IGRMS + DR2IGRMS)/2, pmax(DR1IGRMS, DR2IGRMS))) %>% 
#   ungroup() %>% 
#   select(SEQN, FDCD, GRMS) 


### Groups of foods ###
## First day intake ##


# Food code descriptions
# Ordinary codes
fcd <- read_xpt(here("NHANES/2009_2010/DRXFCD_F.XPT")) %>% 
  mutate(DRXFDCD_char = format(DRXFDCD, digits = 8))
# Modification codes
#fcdl <- read_xpt(here("NHANES/2009_2010/DRXMCD_F.XPT"))




# Food groups defined by FRSG for USDA food codes
fgrp <- read_csv(here("NHANES", "2009_2010", "FRSG_food_group_lookup_FNDSS_5.csv")) %>% 
  transmute(grp_code = `Variable name`, grp_description = `Food group`, food_code_raw = `Food Code Number`) %>% 
  mutate(scode = str_split(food_code_raw, "or")) %>% 
  unnest(scode) %>% 
  transmute(grp_code,
            grp_description,
            ss = str_replace_all(scode, " ", ""),
            fdcd_regex = str_replace_all(str_extract(ss, "^.{8}$"), "-", "."),
            start_range = str_replace_all(str_extract(ss, ".{8}(?=thru)"), "-", "0"),
            end_range = str_replace_all(str_extract(ss, "(?<=thru).{8}"), "-", "9"),
            rng = paste(start_range, end_range)) 

#fgrp %>% print(n=Inf)

# Match food codes to food groups
# Long format
fdd <- fcd %>% 
  # All the food codes
  select(DRXFDCD, DRXFDCD_char, DRXFCSD) %>% 
  
  # Regex matches of food codes to food groups
  left_join(fcd %>% 
              transmute(DRXFDCD, DRXFDCD_char, DRXFCSD) %>% 
              regex_inner_join(fgrp %>% 
                                 select(grp_code, fdcd_regex), by = c("DRXFDCD_char" = "fdcd_regex")) %>%
              
              # Matches of food codes to intervals specified in food group lookup
              bind_rows(fcd %>% 
                          transmute(DRXFDCD, DRXFDCD_char, DRXFCSD) %>% 
                          fuzzy_inner_join(fgrp %>%
                                             select(grp_code, rng),
                                           by = c("DRXFDCD" = "rng"),
                                           match_fun = CheckWithinInterval)),
            by = c("DRXFDCD", "DRXFDCD_char", "DRXFCSD")) %>% 
  mutate(grp_code = if_else(is.na(grp_code), "UNGROUPED", grp_code))


# Changes in food groupings across years
# Draw from changes file in FNDDS notes - relatively minor


# Total intake by group
food_grp1 <- food1 %>% 
  inner_join(fdd %>% 
               select(DRXFDCD, grp_code), 
             by = c("FDCD" = "DRXFDCD")) %>% 
  group_by(SEQN, DRDINT,grp_code) %>% 
  summarise(DR1IGRMS = sum(DR1IGRMS),
            DR1ICARB = sum(DR1ICARB),
            DR1ISUGR = sum(DR1ISUGR))

food_grp2 <- food2 %>% 
  inner_join(fdd %>% 
               select(DRXFDCD, grp_code), 
             by = c("FDCD" = "DRXFDCD")) %>% 
  group_by(SEQN, DRDINT, grp_code) %>% 
  summarise(DR2IGRMS = sum(DR2IGRMS),
            DR2ICARB = sum(DR2ICARB),
            DR2ISUGR = sum(DR2ISUGR))


# Average intake across days 
food_grps_per_day <- food_grp1 %>% 
  full_join(food_grp2, by = c("SEQN", "grp_code", "DRDINT")) %>% 
  # Take the mean if an individual was measured on both days, 
  mutate_at(c("DR1IGRMS", "DR2IGRMS", "DR1ISUGR","DR2ISUGR", "DR1ICARB", "DR2ICARB"), function(.){if_else(is.na(.), 0, .)}) %>% 
  # Otherwise take the existing value (the parallel maximum here)
  mutate(GRMS = if_else(DRDINT == 2, (DR1IGRMS + DR2IGRMS)/2, pmax(DR1IGRMS, DR2IGRMS)),
         CARB = if_else(DRDINT == 2, (DR1ICARB + DR2ICARB)/2, pmax(DR1ICARB, DR2ICARB)),
         SUGR = if_else(DRDINT == 2, (DR1ISUGR + DR2ISUGR)/2, pmax(DR1ISUGR, DR2ISUGR))) %>% 
  # # The simple average
  # mutate(GRMS = (DR1IGRMS + DR2IGRMS)/2,
  #        CARB = (DR1ICARB + DR2ICARB)/2,
  #        SUGR = (DR1ISUGR + DR2ISUGR)/2) %>% 
  ungroup()

# Wide format
food_grps_grms_per_day <- food_grps_per_day %>% 
    select(SEQN, grp_code, GRMS) %>% 
  spread(grp_code, GRMS, fill = 0) 

# food_grps_carbs_per_day <- food_grps_per_day %>% 
#   select(SEQN, grp_code, CARB) %>% 
#   spread(grp_code, CARB, fill = 0) 
# 
# food_grps_sugars_per_day <- food_grps_per_day %>% 
#   select(SEQN, grp_code, SUGR) %>% 
#   spread(grp_code, SUGR, fill = 0) 


# Total intake (not grouped)
# Note that these totals contain both dry matter and liquid
food_total_grms_per_day <- food1 %>% 
  group_by(SEQN) %>% 
  summarise(DRDINT = DRDINT[1], DR1IGRMS = sum(DR1IGRMS)) %>% 
  full_join(food2 %>% 
               group_by(SEQN) %>% 
               summarise(DR2IGRMS = sum(DR2IGRMS)),
             by = "SEQN") %>% 
  # Take the mean if an individual was measured on both days, 
  mutate_at(c("DR1IGRMS", "DR2IGRMS"), function(.){if_else(is.na(.), 0, .)}) %>% 
  # Otherwise take the existing value (the parallel maximum here)
  mutate(GRMS = if_else(DRDINT == 2, (DR1IGRMS + DR2IGRMS)/2, pmax(DR1IGRMS, DR2IGRMS)))
  

cat("\n Dietary data prepared")

### End of dietary data preparation script ###


# Treelet analysis on individual food items not feasible 
# does not compute and tree too large anyway.

# foods <- foods_av %>% 
#   # Only records for cohort members
#   inner_join(select(nhanes, SEQN), by = "SEQN") %>% 
#   # Into wide format
#   spread(key = FDCD, value = GRMS, fill = 0) 
#   
# 
#foods_av %>% filter(FDCD == "11111150") %>% slice(1:10)

## Initial treelets
# # Scale data so can intepret treelets on correlation scale
# foods_scl <- foods %>% 
#   select(-SEQN) %>% 
#   scale()
# 
# # 5 mins
# foods_cor <- cor(foods_scl)
# foods_cor <- cor(foods_scl[,1:1000])
# 
# #foods_scl_red <- foods_scl[,1:100]
# # FUll dataaset
# #12:30 - out of memory
# # 12:43 - first 1000 columns
# 
# 
# # Generate maximum height tree
# # Save basis vectors at all cut levels so the most useful can be selected later
# foods_tree_full <- Run_JTree(foods_cor, maxlev = ncol(foods_cor)-1, whichsave = 1:(ncol(foods_cor)-1))
# # Extract the treelet components and associated variances
# foods_tc_full <- TTVariances(foods_tree_full, foods_cor)
# 
# # Dendrogram of maximum height tree
# foods_dendro <- dendro_data(ConvertTTDendro(foods_tree_full))
# 
# ggdendrogram(foods_dendro, rotate = T)

