
rm(list = ls()) 

library(tidyverse)

###
BHC20_100 <- read_csv("BHC20-100.csv")
BHC20_100 <- BHC20_100 %>% rename(area = "Area (m)")
BHC20_100 <- BHC20_100 %>% rename(BHC20_100 = "BHC20_1")

BHC20_100sum <- BHC20_100 %>% 
  group_by(StationID) %>%
  filter(area == max(area))

###
BHC20_250 <- read_csv("BHC20-250.csv")
BHC20_250 <- BHC20_250 %>% rename(area = "Area (m)")
BHC20_250 <- BHC20_250 %>% rename(BHC20_250 = "BHC20_1")

BHC20_250sum <- BHC20_250 %>% 
  group_by(StationID) %>%
  filter(area == max(area))

###
BHC33_100 <- read_csv("BHC33-100.csv")
BHC33_100 <- BHC33_100 %>% rename(area = "Area (m)")
BHC33_100 <- BHC33_100 %>% rename(BHC33_100 = "BHC33_1")

BHC33_100sum <- BHC33_100 %>% 
  group_by(StationID) %>%
  filter(area == max(area))

###
BHC33_250 <- read_csv("BHC33-250.csv")
BHC33_250 <- BHC33_250 %>% rename(area = "Area (m)")
BHC33_250 <- BHC33_250 %>% rename(BHC33_250 = "BHC33_1")

BHC33_250sum <- BHC33_250 %>% 
  group_by(StationID) %>%
  filter(area == max(area))

###

out <- data.frame(BHC20_100sum$StationID)
out <- out %>% rename(StationID = "BHC20_100sum.StationID")

out <- left_join(out, BHC20_100sum, by = "StationID")

out <- left_join(out, BHC20_250sum, by = "StationID") 

out <- left_join(out, BHC33_100sum, by = "StationID") 

out <- left_join(out, BHC33_250sum, by = "StationID") 

write_csv(out, "BHC-draft.csv")


