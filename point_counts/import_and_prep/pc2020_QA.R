# 2020 Data QA

library(tidyverse)
library(readxl)
library(here)

stations <- read_excel(here("point_counts", "data_raw", "stations_2020_prelim.xlsx"), sheet = 1)

list.files(here("point_counts", "data_processed", "stations_2020_prelim.xlxs"))
