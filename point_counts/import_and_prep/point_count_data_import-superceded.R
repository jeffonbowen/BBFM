########################################################
#### Import Site C Songbird Data, Tidy and QA Check ####
########################################################

# This script is for import, tidy and quality assurance checks of the survey data.


library(tidyverse)
library(readxl)
library(lubridate)
library(skimr)

## Read and tidy data ----

# The data imported below is a copy of the master data file that is located
# in the project sp site. Make sure it is up to date.
# Set some nicer names for variables regularly used.
# Get rid of variables that will not be used in any analysis, just for cleaner tables.

stations <- read_excel("point_counts/data_raw/Songto2019data.xlsx", "Stations") %>% 
  rename(StationID = 'Sample Station Label',
         Easting = `Easting Sample Station`,
         Northing = `Northing Sample Station`) %>%
  mutate(Year = str_sub(`Survey Name`,-4)) %>% 
  select(StationID, Year, everything(), 
         -c(`Sample Station Photos`, `Sample Station Comments`,
            `UTM Zone Sample Station`, `Survey Name`, Station))

surveys <- read_excel("point_counts/data_raw/Songto2019data.xlsx", "Surveys") %>% 
  rename(StationID = 'Sample Station Label') %>% 
  select(StationID, Visit, Date, Time, SurveyDuration) %>% 
  mutate(Date = as_date(Date), 
         Time = as_datetime(Time)) 
# Update the date in `Time` with correct date from `Date` 
date(surveys$Time) <- date(surveys$Date)
# Most analysis of the songbird data must account for temporal survey-level effects
# These variables can be created now and to be used later.
# Year, TSSR, YDAY, DSLS
  
obs <- read_excel("point_counts/data_raw/Songto2019data.xlsx", "Observations",
                  col_types = c("text", "numeric", "numeric", "text", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", "text", "text",
                                "text", "text", "numeric", "numeric", "text", 
                                "numeric", "text", "numeric", "numeric"),
                  na = "NA") %>% 
  rename(StationID = 'Sample Station Label', SpCode = Species) %>% 
  select(StationID, Visit, SpCode, `Count 5 min`, `Count 0-3 min`, `Count 3-5 min`, 
         `Count 5-10 min`, Count, `Distance Category`) %>% 
  filter(str_detect(SpCode, "B-U", negate = TRUE)) # Clear out the unknowns

# Bring in the BC Bird list from BC Species and Ecosystem Explorer
BCbirds <- read_excel("point_counts/data_raw/Songto2019data.xlsx", "BC_Bird_List")
BCbirds <- BCbirds %>% 
  select(ID, `English Name`, `Scientific Name`, `Species Code`, Order, 
         `BC List`, COSEWIC, SARA) %>% 
  rename(SpCode = "Species Code") %>% 
  mutate(IsSongbird = ifelse(Order %in% c("Columbiformes", "Caprimulgiformes", 
                                          "Piciformes", "Passeriformes"), TRUE, FALSE))

## Check for Errors ----

# Check for duplicates in stations$StationID 
stations %>% 
  count(StationID) %>% 
  filter(n>1)

# Check for duplicates in surveys$StationID & Visit
# Note: There are some duplicates in 2011 data. Need to fix (2020-06-29).
surveys %>% 
  count(StationID, Visit) %>% 
  filter(n>1)

# Check to see if there is a match for StationID between `stations` and `surveys`.
# This can be done using a full join and then filter. 

# Join the two tables
StationID_match <- stations %>% 
  full_join(surveys, by = "StationID")
# Filter to check for StationID in `surveys` but not in `stations`.
StationID_match %>% filter(is.na(Year))
# Filter to check for StationID in `stations` but not in `surveys`.
StationID_match %>% filter(is.na(Visit))

# Now check for match for StationID between 'surveys' and 'obs'.
# Note that `stations` and `surveys` must be clean before this is useful.
StationID_match <- surveys %>% 
  full_join(obs, by = c("StationID", "Visit"))
# Filter to check for `StationID $ Visit` in `obs` but not in `surveys`.
StationID_match %>% filter(is.na(Date))
# Filter to check for `StationID $ Visit` in `surveys` but not in `obs`.
# NB: This might be because of no birds observations rather than data entry error.
StationID_match %>% filter(is.na(SpCode))

# Validate species codes in obs table.
obs %>% 
  left_join(BCbirds, by = "SpCode") %>% 
  filter(is.na(`English Name`))

# Check  missing values (nulls).
is.null(stations)
is.null(surveys)
is.null(obs)

## Create summary tables to see if any values out of range ----

# Map station coordinates for anything out of range.
stations %>% 
  ggplot(aes(x = Easting, y = Northing, colour = Year)) + 
  geom_point()

# Check number of surveys per station per year
table(surveys$Visit, year(surveys$Date))
# This is a tidyverse way of creating frequency tables. 
# Longer but good for exporting or further analysis. 
surveys %>% 
  group_by(Visit, Year = year(Date)) %>% 
  summarize(Frequency = n()) %>% 
  pivot_wider(Visit, names_from = Year, values_from = Frequency) %>% 
  replace(is.na(.), 0)  # Could also use mutate with ifelse instead.
# The `Visit` numbering is not correct for all stations. 
# In general it doesn't matter for other analyses but it does for the 
# analysis of the number of surveys. The script below does a proper re-numbering. 
surveys <- surveys %>% 
  select(StationID, Date, Time, Visit,SurveyDuration) %>% 
  arrange(StationID, Date) %>% 
  group_by(StationID) %>% 
  mutate(Visit = row_number())

# The above approach just counts visit numbers, which is ok for QA, but not for reporting.
# Below is an approach to displaying by number of visits.
s.tmp <- surveys %>% 
  count(StationID, Year = year(Date), name  = "Visits") %>% 
  group_by(Year, Visits) %>% 
  summarize(Count = n())
# Table for presentation
s.tmp %>% 
  pivot_wider(id_cols = "Visits", names_from = "Year", values_from = "Count") %>% 
  replace(is.na(.), 0)
# Plot for presentation
s.tmp %>% 
  ggplot(aes(x = as.character(Year), y = Count, 
             fill = fct_rev(factor(Visits)))) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette="reds") +
  xlab("Year") +
  labs(fill='# of Visits')

# Survey dates and times.
table(year(surveys$Date))
table(month(surveys$Date))
table(hour(surveys$Time))

# Observations counts
table(obs$Count)
# Distances
table(obs$`Distance Category`)

# Clean up temp tables
rm(StationID_match)
rm(s.tmp)

# Some extra QA on surveys
qa <- surveys %>% 
  select(StationID, Visit) %>% 
  mutate(Visit = factor(as.character(`Visit`)),
         survey = "x") %>% 
  pivot_wider(id_cols = StationID, names_from = Visit, values_from = survey)


## Create wide format ----

# Prepare Data

-   Will use the 5 minute count column.

-   Since some survey years were fixed radius and some unlimited, select only the obs =\< 100.

-   Filter for songbirds only.

-   Then reformat data to wide format for use in vegan: samples by species.

-   To make the wide matrix, it is necessary to join up with the survey data. This is because there were some surveys with no species detections and therefore the surveys are not in the observations.

Filter for songbirds only and then pivot wider:
  
  ```{r message=FALSE}
# Make list of songbirds to filter data on.
song_list <- BCbirds %>%
  filter(IsSongbird == TRUE)
song_list <-song_list$SpCode

# Filter
songbirds_wide <- obs %>% 
  filter(SpCode %in% song_list) %>% 
  filter(`Distance Category` != ">100") %>%
  group_by(StationID, Visit, SpCode) %>% 
  summarize(count = sum(`Count 5 min`)) %>% 
  pivot_wider(id_cols = c("StationID", "Visit"), 
              names_from = "SpCode",
              values_from = "count",
              values_fill = 0) %>% 
  select(sort(names(.))) %>% 
  select(StationID, Visit, everything())
```

Join to survey data to pick up the extra stations that had no observations.

```{r}

temp <- full_join(songbirds_wide, surveys, by = c("StationID", "Visit")) %>% 
  filter(is.na(StationID))
temp




# Will need this later
replace(is.na(.), 0)

```

