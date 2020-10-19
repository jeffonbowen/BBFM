#########################################
#### Create Site C Survey Covariates ####
#########################################

# This script builds the dataset used to analyze the songbird data.
# The data includes (1) survey-level covariates to account for temporal variation
# in sampling and (2) habitat data. 

library(tidyverse)
library(lubridate)
library(hms)
library(readxl)

# The covariate data is built using the songbird station and survey data as a base 
# and must already be in memory. 
# Run this if songbird data not already in memory.
# source("SongImport.R")

#### TIME AND DATE COVARIATES ####

# Import and create time and date covariates that can be used to correct for 
# survey-level temporal differences. The covariates will be those needed for 
# QPAD. These covariates can be used for non-QPAD analyses as well. 

## Temporal Covariates ----

# This is for Time Since Sunrise (TSSR) covariate.
# Import sunrise data that has been downloaded from NRC for a given location: 
# https://nrc.canada.ca/en/research-development/products-services/software-applications/sun-calculator/
# The NRC data is specific to 2020. Can join data using day of the year (yday) to 
# allow for matches for other years.
sunrise <- read_csv("point_counts/data_raw/sunrisesunsetFSJ.csv", skip = 1) %>% 
  select(SunDate = Date, `Sun rise`) %>% 
  mutate(SunDate = mdy(SunDate), YDAY = yday(SunDate), 
         `Sun rise` = as_hms(`Sun rise`))

# Create a table for survey timing covariates
survey_cov <- surveys %>% 
  select(StationID, Visit, Date, Time, SurveyDuration) %>%
# Year variable. Integer right now. Create factor later when needed.
  mutate(Year = year(Time)) %>%
# For linear effects of year, use monitoring year to avoid scaling problems. 
  mutate(MonYear = Year - 2006) %>% 
# Use day-of-the-year to join the sunrise data and calculate TSSR.
  mutate(YDAY = yday(Time)) %>% 
  left_join(sunrise, by = "YDAY") %>% 
  mutate(TSSR = (as.double(as_hms(Time) - (`Sun rise`)))/60/60/24) %>% 
# Add JDAY (day-of-the-year). This is not needed when DSLS is used.
  mutate(JDAY = YDAY / 365) %>% 
# Continue on with Days Since Local Spring.
# Date of local spring: May 17 for JSJ area. 
# From https://climateatlas.ca/, under agriculture, last frost.
  mutate(DSLS = (YDAY - yday(ymd("2020-05-17")))/365) %>% 
# Clean up intermediate variables
  select(-SunDate, -`Sun rise`, -YDAY)

# Add Easting and then rescale (centre on mean).
survey_cov <- survey_cov %>% 
  left_join(select(stations, StationID, Easting), by = "StationID") %>% 
  mutate(EastWest = scale(Easting)) %>% 
  select(-Easting)    # Clean up intermediate variables.   

# Remove sunrise table since we are done with it. 
rm(sunrise)

#### HABITAT DATA ####

# QPAD requires 
#   (1) tree cover from MODIS VCF 
#   (2) NALCMS for vegetation type (two covariates)
# Best not to wrangle with the spatial data toomuch. Get the RS layers
# clipped and in the right projection from GIS and then develop script for extracting.
# Do this later when needed.

# Other Habitat Data
#   (1) BHC7 and BHC20 from TEM. This was used for the 2019 analyses. 
#   (2) NALCMS data for most dominant class and for continuous values. For further exploration.
# The BHC assignment has to have human validation. Master table outside of R.

# For now, will use the previously generated habitat variables from the 2019 analyses.
# Need to move to progrmatic generation of variables from the spatial layers. 
# The exception is BHC, whcih is generated from a combination of CEM and field data. 
# Wonder if the BHC assignment should be attached to `station` data since it must be 
# assigned and not generated.

## Tree Cover ----

TreeCover <- read_excel("point_counts/data_raw/TREEto2019.xlsx", "Tree") %>% 
  select(StationID, TREE)
survey_cov <- survey_cov %>% 
  left_join(TreeCover, by = "StationID")
# Check for missing values. (2020-07-03) There are two that will be fixed when 
# previous errors fixed. 
filter(survey_cov, is.na(TREE))

## BHC and LCC ---- 
BHC_LCC <- read_excel("point_counts/data_raw/BHCto2019.xlsx", "Habitat") %>%
  select(StationID, BHC20_100, NALCMS100m)
survey_cov <- survey_cov %>% 
  left_join(BHC_LCC, by = "StationID")
filter(survey_cov, is.na(BHC20_100))

# Load the BHC lookup table
BHClookup <- read_excel("point_counts/data_raw/BHClookup.xlsx", "BHClookup")
# Clean up intermediate tables.  
rm(BHC_LCC, TreeCover)

