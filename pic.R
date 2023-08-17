loggers_raw %>% 
    mutate(#tm = 1:nrow(.), 
           # D2 = str_replace_all(D, "/2016", ""),
           # cell = case_when(str_detect(D2, "16|17") ~ "right_down", 
           #                     TRUE ~ "left_down"), 
           D = as.numeric(substr(D, 1,2))) %>% 
    left_join(factors, by = c("D", "H")) %>% 
    mutate(H2 = case_when(D == 7 | D == 17 ~ H, TRUE ~ H*2)) %>% View 
    select(-temp_0cm, -hum_abs)
    
