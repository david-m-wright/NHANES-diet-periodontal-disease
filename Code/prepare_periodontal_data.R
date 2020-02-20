# Analysis of diet and periodontal disease measures from NHANES
# Preparation of NHANES dental examination and periodontal data

library(tidyverse)
library(haven)
library(here)

# Function to drop wave specific 3 letter prefixes from NHANES periodontal variable names
# Arguments:
# dat = data.frame to be renamed
# Value: data.frame with renamed variables
DropPerioPrefix <- function(dat){
    dat %>% rename_at(vars(matches("OHA|OHD|OHX")), list(~str_replace(., "OHA|OHD|OHX", "")))
}


cat("\n Preparing periodontal data")

## Dental examination ##

# Calculate tooth count from dental exam data
dental_raw <- list.files(path = here("NHANES"), pattern = "OHXDEN", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  bind_rows()  

tooth_count_status <- dental_raw %>% 
  # Remove those where a dental exam was not carried out
# Select only the tooth count variables
select(SEQN, matches("^OHX[0-9]{2}TC$")) %>% 
  gather(variable, value, -SEQN) %>% 
  # Remove missing values (variable not assessed at this wave)
  filter(!is.na(value)) %>% 
  mutate(tooth = str_extract(variable, "[0-9]{2}"),
         tooth_present = if_else(value == 1 | value == 2, T, F),
         implant_present = if_else(value == 3, T, F)) 

tooth_count <- tooth_count_status %>% 
  group_by(SEQN) %>% 
  summarise(tooth_count = sum(tooth_present),
            implant_count = sum(implant_present)) 


# tooth_count %>% count(tooth_count) %>% print(n = Inf)
# tooth_count %>% count(implant_count) %>% print(n = Inf)



# Combine periodontal datasets across waves
perio_raw <- list.files(path = here("NHANES"), pattern = "OHXP", recursive = T, full.names = T) %>% 
  lapply(read_xpt) %>% 
  lapply(DropPerioPrefix) %>% 
  bind_rows() %>% 
  # Convert '99' ("Cannot be assessed") measurements to NA
  mutate_at(vars(matches("^[0-9]{2}[A-Z]{3}$")), list(~if_else(.==99, as.numeric(NA), .))) #%>%
  # Count the number of LA measurements taken for each individual
  #mutate(LA_measurements = rowSums(!is.na(select(., matches("^[0-9]{2}"))))) 


### Prepare outcome variables ###

## Clinical attachment loss (CAL) ##

# CAL only reported for 28 teeth as recommended in Holftreter
# Note that sites with implants have an NA value anyway
perio_cal <- perio_raw %>% 
  select(SEQN, matches("^[0-9]{2}LA")) %>% 
  gather(variable, value, -SEQN) %>% 
  mutate(tooth = str_extract(variable, "^[0-9]{2}"),
         site = str_extract(variable, "[A-Z]{3}$"),
         CAL3mm = if_else(value >= 3, T, F),
         CAL5mm = if_else(value >= 5, T, F)) %>% 
  group_by(SEQN) %>% 
  # Number of sites per mouth with CAL
  summarise(CAL_sites_assessed = sum(!is.na(value)),
            CAL_sites3mm = sum(CAL3mm, na.rm = T),
            CAL_sites5mm = sum(CAL5mm, na.rm = T),
            # Proportion of sites per mouth with CAL
            prop_CAL_sites3mm = CAL_sites3mm/CAL_sites_assessed,
            # Mean CAL (mm) per mouth
            CAL_mouth_mean = mean(value, na.rm =T))



## Pocket depth ##

# >= 4mm and >= 6mm are commonly used thresholds
perio_pocket <- perio_raw %>% 
  select(SEQN, matches("^[0-9]{2}PC")) %>% 
  gather(variable, value, -SEQN) %>% 
  mutate(tooth = str_extract(variable, "^[0-9]{2}"),
         site = str_extract(variable, "[A-Z]{3}$"),
         PD4mm = if_else(value >= 4, T, F),
         PD6mm = if_else(value >= 6, T, F)) %>% 
  group_by(SEQN) %>% 
  # Number of sites per mouth with PD over threshold
  summarise(PD_sites_assessed = sum(!is.na(value)),
            PD_sites4mm = sum(PD4mm, na.rm = T),
            PD_sites6mm = sum(PD6mm, na.rm = T),
            # Proportion of sites per mouth with PD
            prop_PD_sites4mm = PD_sites4mm/PD_sites_assessed,
            # Mean PD (mm) per mouth
            PD_mouth_mean = mean(value, na.rm =T))



rm(dental_raw, tooth_count_status)
cat("\n Periodontal data prepared")
