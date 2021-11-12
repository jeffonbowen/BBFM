## Process New Habitat Data

library(readxl)

# Note later. If Shrubland is of interest, then need to adjust for clearcut.

nalcms <- read_xlsx(here("point_counts", "data_raw", "Songto2021_NALCMS_100.xlsx"),
                    sheet = 1) %>% 
  rename(StationID = SampleStationLabel) %>% 
  pivot_wider(c(StationID, CC_flag), names_from = Description, values_from = `Area (m2)`,
              values_fill = 0) %>% 
  rename(Conifer_lc = "Temperate or sub-polar needleleaf forest",
         Decid_lc = "Temperate or sub-polar broadleaf deciduous forest", 
         Mixed_lc = "Mixed forest", 
         Shrubland_lc = "Temperate or sub-polar shrubland",
         Grassland_lc = "Temperate or sub-polar grassland",
         Urban_lc = "Urban and built-up",
         Wetland_lc = Wetland,
         Cropland_lc = Cropland,
         Barren_lc = `Barren land`,
         Water_lc = Water) %>% 
  mutate(across(Conifer_lc:`Barren_lc`, ~.x / 31400)) %>%
  mutate(Conifer_lc = if_else(CC_flag == "y", 0, Conifer_lc),
         Decid_lc = if_else(CC_flag == "y", 0, Decid_lc),
         Mixed_lc = if_else(CC_flag == "y", 0, Mixed_lc))

# The CEM raw data comes from raster>zonal histogram in QGIS. Note that
# this cannot be done with overlapping zones in ArcMap (not easily).

CEM <- read_xlsx(here("point_counts", "data_raw", "CEM_out.xlsx")) %>% 
  rename(NVAC = HISTO_1, Upl_Forest = HISTO_2, Gram_Wetland = HISTO_3,
         Water = HISTO_4, For_Wetland = HISTO_5, Rip_Forest = HISTO_6,
         Other = HISTO_7, Rip_Wetland = HISTO_8, Grassland = HISTO_9,
         Asp_Shrub = HISTO_10, Unk = HISTO_11) %>% 
  mutate(Wetland = Gram_Wetland + For_Wetland,
         Dry_Slopes = Grassland + Asp_Shrub,
         Forest_All = Upl_Forest + Rip_Forest) %>% 
  select(StationID, NVAC:Forest_All) %>% 
  mutate(across(!StationID, ~.x * 25 / 31400)) %>% 
  left_join(nalcms, by = "StationID") %>% 
  mutate(Upl_Forest = if_else(CC_flag == "y", 0, Upl_Forest),
         Rip_Forest = if_else(CC_flag == "y", 0, Rip_Forest),
         Forest_All = if_else(CC_flag == "y", 0, Forest_All)
         )

env_cont <- CEM

write_csv(env_cont, here("point_counts", "data_processed", "env_continuous.csv"))
