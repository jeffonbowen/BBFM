

temp <- obs %>% 
  pivot_wider(id_cols = c("StationID", "Visit"), 
              names_from = "SpCode", 
              values_from = "Count")

temp %>% 
  group_by("StationID", "Visit") %>% 
  
  
  
  
x <- tree_dat %>% 
  group_by(Structural_Stage, SMR) %>% 
  summarize(n = n_distinct(HabDat_ID)) %>% 
  pivot_wider(names_from = "SMR", values_from = "n")  
x
