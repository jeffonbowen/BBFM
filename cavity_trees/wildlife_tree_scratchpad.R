

temp <- obs %>% 
  pivot_wider(id_cols = c("StationID", "Visit"), 
              names_from = "SpCode", 
              values_from = "Count")

temp %>% 
  group_by("StationID", "Visit") %>% 
  