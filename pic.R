library(ggh4x) # dendrograms as well

for_pic <- loggers_raw %>% 
    mutate(D = as.numeric(substr(D, 1,2))) %>% 
    left_join(select(factors, D, H, lux), by = c("D", "H")) %>% 
    mutate(H2 = case_when(D == 7 | D == 16 ~ H, 
                          D == 8 | D == 17 ~ H+24, 
                          TRUE ~ H+48), 
           cell1 = case_when(D < 15 ~ "07-08.08", TRUE ~ "16-17.08"),
           .after = H) %>% 
    select(D, H, H2, cell1, 
           `1. Light level, klx` = lux,
           `2. Temperature, °C` = temp_2m, 
           `3. Humidity, %` = hum_rel
           )

b <- min(for_pic$`3. Humidity, %`)
for_pic <- for_pic %>% 
    mutate(`3. Humidity, %_2` = `3. Humidity, %`-b) 
a <- max(for_pic$`2. Temperature, °C`) / max(for_pic$`3. Humidity, %_2`)
for_pic <- for_pic %>% 
    mutate(`3. Humidity, %_3` = `3. Humidity, %_2`*a)

p_down <- for_pic %>% 
    ggplot(aes(x = H2)) + 
    geom_line(aes(y = `2. Temperature, °C`), color = "red") + 
    geom_point(aes(y = `2. Temperature, °C`), color = "red", size = 0.5) +
    geom_line(aes(y = `3. Humidity, %_3`), color = "blue") +
    geom_point(aes(y = `3. Humidity, %_3`), color = "blue", size = 0.5) +
    scale_y_continuous(
        sec.axis = sec_axis(trans=~./a+b, name="3. Humidity, %")) +
    facet_wrap(~cell1) + 
    labs(x = NULL)

p_top <- for_pic %>% 
    filter(`1. Light level, klx` > 0) %>% 
    ggplot(aes(H2, `1. Light level, klx`)) + 
    geom_line(color = "yellow") + 
    geom_point(color = "yellow") + 
    facet_wrap(~cell1) + 
    scale_y_log10() +
    labs(x = NULL)

gridExtra::grid.arrange(p_top, p_down, ncol = 1)










