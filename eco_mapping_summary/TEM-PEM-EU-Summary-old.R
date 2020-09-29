## Ecosystem Unit Summary Tool ## 

# This tool creates a list of ecosystem units and calculates area mapped 
# across all deciles.

# An ecosystem unit is the combination of site series/map code, seral 
# association and structural stage, though could also include modifiers 
# and other attributes.

# Step 1. Import and organize ----

# Libraries
library(tidyverse)
library(readxl)
library(here)

here("CEM-footprint.csv")

# Read attribute file.
# The attribute file is a csv export the original shapefile or other source. 
eco_attributes <- read_csv("dat/CEM-footprint.csv", 
                  col_types = cols(.default = "c",
                                   'Shape_Area' = col_double(),
                                   'SDEC_1' = col_integer(), 
                                   'SDEC_2' = col_integer(), 
                                   'SDEC_3' = col_integer()))

# Select the attributes of interest.
# Customize list based on dataset. 

eco_attributes <- select(eco_attributes, PAZ, CEM_ECP, ECO_SEC, BGC_ZONE, BGC_SUBZON, BGC_VRT, BGC_PHASE, SDEC_1, 
        SITE_S1, SITEMC_S1, SITE_M1A, SITE_M1B, STRCT_S1, STRCT_M1, SERAL_1, SDEC_2, 
        SITE_S2, SITEMC_S2, SITE_M2A, SITE_M2B, STRCT_S2, STRCT_M2, SERAL_2, SDEC_3, 
        SITE_S3, SITEMC_S3, SITE_M3A, SITE_M3B, STRCT_S3, STRCT_M3, SERAL_3, Shape_Area) 

# Conduct QA of codes and values
table(eco_attributes$BGC_ZONE, eco_attributes$BGC_SUBZON)

table(eco_attributes$SDEC_1)
table(eco_attributes$SITEMC_S1)
table(eco_attributes$SERAL_1)
table(eco_attributes$STRCT_S1)

table(eco_attributes$SDEC_2)
table(eco_attributes$SITEMC_S2)
table(eco_attributes$SERAL_2)
table(eco_attributes$STRCT_S2)

table(eco_attributes$SDEC_3)
table(eco_attributes$SITEMC_S3)
table(eco_attributes$SERAL_3)
table(eco_attributes$STRCT_S3)

# Site C CEM data error found: Recode the "AS" to "as" in SERAL_1
eco_attributes[eco_attributes$SERAL_1 == "AS" & !is.na(eco_attributes$SERAL_1), ]
eco_attributes[eco_attributes$SERAL_1 == "AS" & !is.na(eco_attributes$SERAL_1), ]$SERAL_1 <- "as"
# If there is more than just a few edits, then perhaps the original data needs to be corrected.

# Step 2. Create summary across all deciles ----
 
d1 <- eco_attributes %>% 
  select(PAZ,
         BGC_ZONE,
         BGC_SUBZON,
         BGC_VRT,
         "Map Code" = SITEMC_S1, 
         "Seral Association" = SERAL_1, 
         "Structural Stage" = STRCT_S1,
         "Decile" = SDEC_1, 
         Shape_Area) %>% 
  mutate("Area" = Decile/10*Shape_Area)
d2 <- eco_attributes %>% 
  select(PAZ,
         BGC_ZONE,
         BGC_SUBZON,
         BGC_VRT,
         "Map Code" = SITEMC_S2, 
         "Seral Association" = SERAL_2, 
         "Structural Stage" = STRCT_S2,
         "Decile" = SDEC_2, 
         Shape_Area) %>% 
  mutate("Area" = Decile/10*Shape_Area)
d3 <- eco_attributes %>% 
  select(PAZ,
         BGC_ZONE,
         BGC_SUBZON,
         BGC_VRT,
         "Map Code" = SITEMC_S3, 
         "Seral Association" = SERAL_3, 
         "Structural Stage" = STRCT_S3,
         "Decile" = SDEC_3, 
         Shape_Area) %>% 
  mutate("Area" = Decile/10*Shape_Area)

eu_summary <- bind_rows(d1, d2, d3) %>% 
  mutate("BEC Variant" = str_c(`BGC_ZONE`, `BGC_SUBZON`, `BGC_VRT`)) %>% 
  group_by(`PAZ`, `BEC Variant`, `Map Code`, `Seral Association`, `Structural Stage`) %>% 
  summarise("Total (ha)" = sum(Area)/10000) %>% 
  mutate("Proportion" = round(`Total (ha)`/sum(`Total (ha)`),3)) %>% 
  mutate("Percent of Total" = Proportion * 100) %>% 
  mutate(Ecosystem = ifelse(is.na(`Seral Association`), 
                             `Map Code`,
                             str_c(`Map Code`, ":", `Seral Association`))) %>% 
  select(1:3, `Ecosystem`, everything()) 
  
write_csv(eu_summary, "out/Footprint_CEM_EU.csv")

mc_summary <- eu_summary %>% 
  group_by(`PAZ`, `BEC Variant`, `Ecosystem`) %>%
  summarize("Total (ha)" = sum(`Total (ha)`)) %>%
  mutate("unit" = str_c(`BEC Variant`, Ecosystem)) %>% 
  filter(!is.na(Ecosystem)) %>% 
  pivot_wider(id_cols = c(`BEC Variant`, Ecosystem, unit), names_from = PAZ, values_from = `Total (ha)`) 

mc_summary <- mc_summary %>% 
  left_join(eco_ref, by = "unit") %>% 
  arrange(Sort) %>% 
  select(Sort, "BEC Variant" = `BEC Variant.x`, Ecosystem, `Site Series`, `Site Series Name`,
         Type, Dam, Reservoir, "TL" = `Transmission Line`, `Hwy 29`, Road, Quarry, 
         "Erosion Impact" = `Erosion Impact Area`)

write_csv(mc_summary, "out/Footprint_CEM_MC.csv")

write_csv(mc_summary, "out/Hwy29.csv")


## For later when using BHC

level_order <- c("CSH","CYF","CMF","DSH","DYF", "DMF", "RSH", "RYF", "RMF", "FBS", "FBT", "WGR", "WSH", "WRI", "DSG", "DSS", "CUL", "NVE", "ANT", "WAT")
g <- ggplot(summary, aes(factor(BHC20, level = level_order), Total)) +
      geom_col()
g




# Step X. Create new EU field ----

# It can be useful to create ecosystem unit fields that can be joined back 
# with the orignal spatial data for easier viewing, labelling and query of polygon EUs 
# in GIS, Google Earth or other spatial application.

eco_attributes <- eco_attributes %>% 
  mutate(MC1 = ifelse(is.na(SERAL_1), 
                      SITEMC_S1,
                      str_c(SITEMC_S1, ":", SERAL_1)))
unique(eco_attributes$MC1)


# Read the list of BHCs for names and sorting
BHC <- read_csv('BHC.csv', col_types = cols('HC Sort' = col_integer()))
BHC <- BHC %>% rename(BHC20 = 'Habitat Code', Sort = 'HC Sort')



                                                