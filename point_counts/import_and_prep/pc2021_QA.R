# 2021 Data QA

library(tidyverse)
library(readxl)
library(here)
library(lubridate)

stations <- read_excel(here("point_counts", "data_raw", "Site C Songbird Data 2021.xlsx"), sheet = 1)


stations <- read_excel(here("point_counts", "data_raw", "Site C Songbird Data 2021.xlsx"), "Station") %>% 
  rename(StationID = 'Sample Station Label',
         Easting = `Easting Sample Station`,
         Northing = `Northing Sample Station`) %>%
  mutate(Year = str_sub(`Survey Name`,-4)) %>% 
  select(StationID, Year, everything(), 
         -c(`Sample Station Photos`, `Sample Station Comments`,
            `UTM Zone Sample Station`, `Survey Name`, Station))

surveys <- read_excel(here("point_counts", "data_raw", "Site C Songbird Data 2021.xlsx"), 
                      "Surveys") %>% 
  rename(StationID = 'Sample Station Label') %>% 
  select(StationID, Visit, Date, Time, SurveyDuration) %>% 
  mutate(Date = as_date(Date), 
         Time = as_datetime(Time)) 
# Update the date in `Time` with correct date from `Date` 
date(surveys$Time) <- date(surveys$Date)
# Most analysis of the songbird data must account for temporal survey-level effects
# These variables can be created now and to be used later.
# Year, TSSR, YDAY, DSLS

obs <- read_excel(here("point_counts", "data_raw", "Site C Songbird Data 2021.xlsx"), 
                  "Observations",
                  col_types = c("text", "numeric", "text", "numeric", "text", "numeric",
                                "numeric", "numeric", "numeric", "numeric", 
                                "text", "text", "text", "text", "numeric", 
                                "numeric", "text", "guess", "guess", "guess", 
                                "guess"),
                  na = "NA") %>% 
  rename(StationID = 'Sample Station Label', SpCode = Species) %>% 
  select(StationID, Visit, SpCode, `Count 5 min`, `Count 0-3 min`, `Count 3-5 min`, 
         `Count 5-10 min`, Count, `Distance Category`, Flyovers) %>% 
  filter(str_detect(SpCode, "B-U", negate = TRUE)) # Clear out the unknowns

# Bring in the BC Bird list from BC Species and Ecosystem Explorer
BCbirds <- read_excel(here("point_counts", "data_raw", "Songto2019data.xlsx"), 
                      "BC_Bird_List") %>% 
  select(ID, `English Name`, `Scientific Name`, `Species Code`, Order, 
         `BC List`, COSEWIC, SARA) %>% 
  rename(SpCode = "Species Code") %>% 
  mutate(IsSongbird = ifelse(Order %in% c("Columbiformes", "Caprimulgiformes", 
                                          "Piciformes", "Passeriformes"), 
                             TRUE, 
                             FALSE))

# No go run the code under import