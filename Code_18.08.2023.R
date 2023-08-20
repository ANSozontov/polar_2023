library(tidyverse)
theme_set(theme_bw() + theme(legend.position = "bottom"))

# custom functions --------------------------------------------------------
per2 <- function(data, periodicity_freq){
	data <- as.numeric(data)
	# Compute the power spectral density
	psd <- spec.pgram(data, log = "no", detrend = FALSE)

	# Compute the spectral density ratio (SDR)
	denom_freq_range <- c(periodicity_freq*0.9, periodicity_freq*1.1)
	
	denom_power <- mean(psd$spec[psd$freq >= denom_freq_range[1] | psd$freq <= denom_freq_range[2]])
	num_power <- psd$spec[abs(psd$freq - periodicity_freq) == min(abs(psd$freq - periodicity_freq))]
	sdr <- num_power/denom_power

	# Perform an F-test to test the significance of the SDR
	num_df <- 1
	denom_df <- psd$df
	num_mean <- num_power
	denom_mean <- denom_power
	f_stat <- (num_mean/denom_mean)/((num_df/num_mean)+(denom_df/denom_mean))
	p_val <- pf(f_stat, num_df, denom_df, lower.tail = FALSE)
	
	# Results
	list(x = data, total = sum(data), SDR = sdr, 
		F.statistic = f_stat, psd_df = psd$df, p.value = p_val)
}

decompress <- function(a){
    # Export decomposed time series as tibble columns
	var.name <- deparse(substitute(a))
	a |> 
		ts(frequency = 6) |> 
		decompose() |> 
		_[1:4] |> 
		as_tibble() |> 
		rename_with(~paste0(var.name, "_", .x))
}

# data load ---------------------------------------------------------------
loggers_raw <- readxl::read_excel("Data_28.04.2023.xlsx", 
		sheet = "logger_data", skip = 1)
loggers <- loggers_raw |>
	mutate(D = as.numeric(substr(D, 1, 2)), 
		  H = (floor(H/4)+1)*4) |>
	group_by(D, H) |>
	summarise_all(mean) |>
	ungroup() |>
	filter(D != 9, D != 18)

factors <- readxl::read_excel("Data_28.04.2023.xlsx", 
		sheet = "factors", range = "A2:E26") |>
	mutate(D = substr(str_squish(D), 1, 2), 
		  D = as.numeric(D), 
		  log_lux = log10(lux), .after = 2)
factors <- left_join(factors, loggers, by = c("D", "H")) 

df <- readxl::read_excel("Data_28.04.2023.xlsx", 
			    sheet = "animals", skip = 3) |>
	mutate_all(as.character) |>
	pivot_longer(names_to = "dt", values_to = "abu", -1:-2) |> 
	separate(dt, sep = "_", into = c("trap_id", "D", "h")) |> 
	separate(h , sep = "\\:", into = "H", extra = "drop") |> 
	mutate_at(c("trap_id", "D", "H", "abu"), as.numeric)

# Environmental variables graph -------------------------------------------
#library(ggh4x) # dendrograms as well
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
           `2. Temperature, °C` = temp_0cm, 
           `3. Humidity, %` = hum_rel
    )

b <- min(for_pic$`3. Humidity, %`)
for_pic <- for_pic %>% 
    mutate(`3. Humidity, %_2` = `3. Humidity, %`-b) 
a <- max(for_pic$`2. Temperature, °C`) / max(for_pic$`3. Humidity, %_2`)
for_pic <- for_pic %>% 
    mutate(`3. Humidity, %_3` = `3. Humidity, %_2`*a)

(p2_down <- for_pic %>% 
    ggplot(aes(x = H2)) + 
        # 
    geom_vline(xintercept = seq(12, 48, by = 24), color = "red", 
               linetype = "dotted", alpha = 0.5) +
    geom_vline(xintercept = seq(0, 48, by = 24), color = "blue", 
               linetype = "dotted", alpha = 0.5) +
        #
    geom_line(aes(y = `2. Temperature, °C`), color = "#F8766D") + 
    geom_point(aes(y = `2. Temperature, °C`), color = "#F8766D", size = 1) +
    geom_line(aes(y = `3. Humidity, %_3`), color = "#4371C4") +
    geom_point(aes(y = `3. Humidity, %_3`), color = "#4371C4", size = 1) +
    scale_y_continuous(
        sec.axis = sec_axis(trans=~./a+b, name="3. Humidity, %")) +
    facet_wrap(~cell1) + 
        #
    theme_classic() +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 48, by = 4), 
                       labels = c(rep(seq(0, 20, by = 4), 2), 0)) + 
        #
    labs(x = NULL))
ggsave("Fig.2 Environmental variables, down section_raw.pdf", 
       plot = p2_down, width = 7.305, height = 2.5)
(p2_top <- for_pic %>% 
    filter(`1. Light level, klx` > 0) %>% 
    ggplot(aes(H2, `1. Light level, klx`)) + 
    geom_vline(xintercept = seq(12, 48, by = 24), color = "red", 
                   linetype = "dotted", alpha = 0.5) +
    geom_vline(xintercept = seq(0, 48, by = 24), color = "blue", 
                   linetype = "dotted", alpha = 0.5) +
    geom_line(color = "#FFBF00") + 
    geom_point(color = "#FFBF00", size = 1) + 
    facet_wrap(~cell1) + 
    scale_y_log10() +
    theme_classic() +
    theme(panel.grid.minor = element_blank(), 
              panel.grid.major = element_blank()) + 
    scale_x_continuous(breaks = seq(0, 48, by = 4), 
                           labels = c(rep(seq(0, 20, by = 4), 2), 0)) + 
    labs(x = NULL))
ggsave("Fig.2 Environmental variables, top section_raw.pdf", 
       plot = p2_top, width = 7, height = 2.5)

# p2 <- gridExtra::grid.arrange(p2_top, p2_down, ncol = 1)
# ggsave("Fig.2 Environmental variables_raw.pdf", plot = p2, width = 7, height = 4.5)

# H1: is total activity periodic? -------------------------------------------------------
abu <- df %>% 
	filter(taxatype == "general") %>% 
	group_by(D, H) %>% 
	summarise(abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
	pull(abu) %>% 
	ts(frequency = 6) %>% 
	decompose()

abu[c("x", "seasonal", "trend", "random")] %>% 
	as.data.frame() %>% 
	rename(observed = x, periodic = seasonal) %>% 
	mutate(dd = rep(1:(nrow(.)/6), each = 6), 
		  hh = rep(1:6, 4),
		  xx = paste0("d", dd, "_t", hh)) %>% 
	rownames_to_column("id") %>% 
	mutate(id = as.numeric(id), hh = as.factor(hh)) %>% 
	pivot_longer(names_to = "component", values_to = "val", -c("id", "xx", "dd", "hh")) %>% 
	mutate(component = factor(component,
						 levels = c("observed","trend", "periodic", "random"))) %>% 
	ggplot(aes(x = id, y = val)) +
	geom_vline(xintercept = (1:4)*6-3, color = "red", linetype = "dotted") +
	geom_vline(xintercept = (1:4)*6, color = "blue", linetype = "dotted") +
	geom_line() + 
	geom_point(shape = 21, fill = "white") +
	facet_grid(rows = vars(component)) + #, scales = "free") + 
	scale_x_continuous(breaks = 1:24, 
		labels = rep(1:6*4, 4)) + 
	theme(axis.text.x = element_text(angle = 0), 
		 panel.grid.minor = element_blank(), 
		 panel.grid.minor.y = element_blank(), 
		 panel.grid.major.x = element_blank()) +
	labs(x = "Sampling hours", y = "Abundance/activity and its components")

ggsave("Fig.3. Total abundance_raw.pdf", width = 7, height = 4)
# Time series of diurnal arthropods activity 
# and its components (trend, periodicity, residuals). 
# Red line marks noon, blue line marks midnights

per2(abu$x, 1/6)

# Signal to noise ratio is  > 4
# Н0 about non-periodicity is rejected with p<0.001 significance level

for_boxplots <- df %>% 
    filter(taxatype == "general") %>% 
    group_by(trap_id, D, H) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    mutate(taxa = "All_taxa", .before = 1) %>% 
    rbind(select(df, -taxatype)) %>% 
    filter(taxa %in%  c("All_taxa", "Collembola", "Aranei", "Acari", "Diptera", 
        "Masikia indistincta","Brachystomella parvula (Schäffer, 1896)", 
        "Isotomurus chaos Potapov et Babenko, 2011", 
        "Pachyotoma crassicauda (Tullberg, 1871)")) %>% 
    separate(taxa, c("gen", "sp"), " ", extra = "drop") %>%
    mutate(taxa = case_when(is.na(sp) ~ gen, TRUE ~ paste0(gen, " ", sp)),
           taxa = fct_inorder(taxa, ordered = TRUE), 
           H2 = factor(H),
           H,
           .before = 1,
           .keep = "unused")

for_boxplots %>% 
    filter(taxa == "All_taxa") %>% 
    group_by(H) %>% 
    mutate(mean_abu = mean(abu)) %>% 
    ungroup() %>% 
    ggplot() + 
    geom_boxplot(aes(x = H, y = abu, group = H2, fill = H2)) +
    geom_line(aes(x = H, y = mean_abu), color = "blue") + 
    geom_point(aes(x = H, y = mean_abu), shape = 15, color = "blue") +
    theme(legend.position = "none") + 
    labs(x = "Sampling hour", y = "Abundance/Activity", 
         subtitle = "All taxa") +
    scale_x_continuous(breaks = 1:6*4,
                       labels = c("04", "08", "12", "16", "20","24")) +
    scale_fill_manual(values = colorRampPalette(c("darkgrey", "orange"))
                      (6)[c(2, 4, 6, 5, 3, 1)])
ggsave("Fig.4. Part A_raw.pdf", height = 3, width = 3)
# Average activity of terrestrial arthropods at different times of the day 1

for_boxplots %>% 
    filter(taxa %in% c("Collembola", "Aranei", "Acari", "Diptera")) %>% 
    group_by(taxa, H, H2) %>% 
    mutate(mean_abu = mean(abu)) %>% 
    ungroup() %>% 
    ggplot() + 
    geom_boxplot(aes(x = H, y = abu, group = H2, fill = H2)) +
    geom_line(aes(x = H, y = mean_abu), color = "blue") +
    geom_point(aes(x = H, y = mean_abu), shape = 15, color = "blue") +
    theme(legend.position = "none") + 
    labs(x = "Sampling hour", y = "Abundance/Activity", 
         subtitle = "By orders") +
    scale_x_continuous(breaks = 1:6*4,
                      labels = c("04", "08", "12", "16", "20","24")) +
    scale_fill_manual(values = colorRampPalette(c("darkgrey", "orange"))
                      (6)[c(2, 4, 6, 5, 3, 1)]) + 
    facet_wrap(~taxa, scales = "free_y")
ggsave("Fig.4. Part B_raw.pdf", height = 6, width = 6)
# Average activity of terrestrial arthropods at different times of the day 2

for_boxplots %>% 
    filter(taxa %in% c("Masikia indistincta",
                       "Brachystomella parvula", 
                       "Isotomurus chaos",
                       "Pachyotoma crassicauda")) %>% 
    group_by(taxa, H, H2) %>% 
    mutate(mean_abu = mean(abu)) %>% 
    ungroup() %>% 
    ggplot() + 
    geom_boxplot(aes(x = H, y = abu, group = H2, fill = H2)) +
    geom_line(aes(x = H, y = mean_abu), color = "blue") +
    geom_point(aes(x = H, y = mean_abu), shape = 15, color = "blue") +
    theme(legend.position = "none") + 
    labs(x = "Sampling hour", y = "Abundance/Activity", 
         subtitle = "By dominant species") +
    scale_x_continuous(breaks = 1:6*4,
                       labels = c("04", "08", "12", "16", "20","24")) +
    scale_fill_manual(values = colorRampPalette(c("darkgrey", "orange"))
                      (6)[c(2, 4, 6, 5, 3, 1)]) + 
    facet_wrap(~taxa, scales = "free_y") + 
    theme(strip.text = element_text(face = "italic"))
ggsave("Fig.4. Part C_raw.pdf", height = 6, width = 6)
# Average activity of terrestrial arthropods at different times of the day 3

# H#1 by taxa -------------------------------------------------------------
by_taxa <- df %>% 
	filter(taxatype == "general", taxa %in% c("Collembola", "Aranei", "Acari", "Diptera")) %>% 
	split(.$taxa) %>% 
	map(~.x %>% 
		group_by(D, H) %>% 
		summarise(abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
		pull(abu) %>% 
		ts(frequency = 6) %>% 
		decompose()
	)
by_taxa %>% 
	transpose %>% 
	pluck("x") %>% 
	lapply(as.numeric) %>% 
	lapply(function(a){
		per2(a, 1/6)
	})

by_taxa %>% 
    lapply(function(a){
		cbind(
		    x = a$x, trend = a$trend, 
			seasonal = a$seasonal, random = a$random) %>% 
		as_tibble %>% 
		mutate_all(function(b){b/max(b, na.rm = TRUE)}) %>% 
		rownames_to_column("id") %>% 
		mutate(id = as.numeric(id))}) %>%  
	map_df(rbind, .id = "taxa") %>% 
	rename(observed = x, periodic = seasonal) %>% 
	pivot_longer(names_to = "component", values_to = "val", -1:-2) %>% 
	mutate(component = factor(component,
		levels = c("observed","trend", "periodic", "random"))) %>% 
	ggplot(aes(x = id, y = val, color = taxa)) +
	geom_vline(xintercept = (1:4)*6-3, color = "red", linetype = "dotted") +
	geom_vline(xintercept = (1:4)*6, color = "blue", linetype = "dotted") +
	geom_line() + 
	geom_point(shape = 21, fill = "white") +
	theme_bw() +
	facet_grid(rows = vars(component), scales = "free") + 
	scale_x_continuous(breaks = 1:24, labels = rep(1:6*4, 4)) + 
	theme(axis.text.x = element_text(angle = 0), 
		 panel.grid.minor = element_blank(), 
		 panel.grid.minor.y = element_blank(), 
		 panel.grid.major.x = element_blank(), 
		 legend.position = "bottom") +
	labs(x = "Sampling hours", y = "Abundance/activity and its components", 
		subtitle = "all variables are normalized to max = 1")
ggsave("Fig.5 Orders_raw.pdf", width = 7, height = 4.5)
# Time series of diurnal arthropods activity (by orders) 
# and its components (trend, periodicity, residuals). 
# Red line marks middays, blue line marks midnights. 
# Note: all raw variables data were normalized to max=1

# H#1 by dominant species  -------------------------------------------------------
dom_sp <- df %>%
	filter(taxa %in% c("Masikia indistincta",
		"Brachystomella parvula (Schäffer, 1896)", 
		#"Isotomurus stuxbergi (Tullberg, 1876)", # check
		#"Erigone psychrophila", # check
		"Isotomurus chaos Potapov et Babenko, 2011",
		"Pachyotoma crassicauda (Tullberg, 1871)")) %>% 
    separate(taxa, c("gen", "sp"), " ", extra = "drop") %>% 
	unite(taxa, gen, sp, sep = " ") %>% 
	split(.$taxa) %>% 
	map(~.x %>% 
	    	group_by(D, H) %>% 
	    	summarise(abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
	        pull(abu) %>% 
	    	ts(frequency = 6) %>% 
	    	decompose()
	)
dom_sp %>% 
	transpose %>% 
	pluck("x") %>% 
	lapply(as.numeric) %>% 
	lapply(function(a){
		per2(a, 1/6)
	})

dom_sp %>% 
    lapply(function(a){
        cbind(x = a$x, trend = a$trend, 
              seasonal = a$seasonal, random = a$random) %>% 
            as_tibble %>% 
            mutate(random = random/max(x), 
                   seasonal = seasonal/max(x), 
                   trend = trend/max(x), 
                   x = x/max(x)) %>% 
            rownames_to_column("id") %>% 
            mutate(id = as.numeric(id))}) %>% 
    map_df(rbind, .id = "taxa") %>% 
    rename(observed = x, periodic = seasonal) %>% 
    pivot_longer(names_to = "component", values_to = "val", -1:-2) %>% 
    mutate(component = factor(component,
                              levels = c("observed","trend", "periodic", "random"))) %>% 
    ggplot(aes(x = id, y = val, color = taxa)) +
    geom_vline(xintercept = (1:4)*6-3, color = "red", linetype = "dotted") +
    geom_vline(xintercept = (1:4)*6, color = "blue", linetype = "dotted") +
    geom_line() + 
    geom_point(shape = 21, fill = "white") +
    theme_bw() +
    facet_grid(rows = vars(component), scales = "free") + 
    scale_x_continuous(breaks = 1:24, labels = rep(1:6*4, 4)) + 
    theme(axis.text.x = element_text(angle = 0), 
          panel.grid.minor = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_blank(), 
          legend.text = element_text(face = "italic"),
          legend.position = "bottom") +
    labs(x = "Sampling hours", y = "Abundance/activity and its components", 
         subtitle = "Note: all variables are normalized to max=1")
ggsave("Fig.6 Dominant species_raw.pdf", width = 7, height = 4.5)
# Time series of arthropods dominant species 
# and its components (trend, periodicity, residuals). 
# Red line marks middays, blue line marks midnights. 
# Note: all raw variables data awere normalized to max=1

# predictors multicoll ---------------------------------------------------
factors[,-1:-2] %>% 
	cor(method = "spearman") %>% 
	round(2)

cor.val1 <- cor(factors[,-1:-3], method = "spearman")
cor.pval1 <- expand_grid(v1 = 4:10, v2 = 4:10) %>% 
	split(1:nrow(.)) %>% 
	lapply(function(a){
		p = cor.test(as_vector(factors[,a$v1]), as_vector(factors[,a$v2]), 
				   method = "spearman")
		data.frame(p = p$p.value, 
				 v1 = colnames(factors)[a$v1], 
				 v2 = colnames(factors)[a$v2])
	}) %>% 
	map_df(tibble) %>% 
	mutate(p = p.adjust(p, method = "BY")) %>% 
	pivot_wider(names_from = v2, values_from = p) %>% 
	column_to_rownames("v1") %>% 
	as.matrix()

colnames(cor.pval1) <- 
    rownames(cor.pval1) <- 
    colnames(cor.val1) <- 
    rownames(cor.val1) <- c(
        "light level", "wind_20cm", "wind_2m", "temperature_0cm", 
        "temperature_2m", "humidity_relative", "humidity_absolute"
    )

pdf("Fig.7. Factor_multicollin_raw.pdf", width = 4, height = 4)
corrplot::corrplot(
	corr = cor.val1, 
	p.mat = cor.pval1, 
	type="upper", 
	order = "original",
	diag = FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off() 
# Spearman correlation matrix. 
# The Spearman coefficient value is shown on the scale on the right. 
# The significance level (α) = 0.05. 
# Statistically non-significant correlations are crossed out

# correlations ------------------------------------------------------------
df.forcor <- df %>% 
	filter(taxa %in% c("Collembola", "Aranei", "Acari", "Diptera")) %>% 
	group_by(D, H, taxa) %>% 
	summarise(abu = sum(abu), .groups = "drop") %>% 
	pivot_wider(names_from = taxa, values_from = abu) %>% 
    mutate(All = abu$x, .after = 2) 
df.forcor <- df %>%
    filter(taxa %in% c("Masikia indistincta",
                       "Brachystomella parvula (Schäffer, 1896)", 
                       "Isotomurus chaos Potapov et Babenko, 2011",
                       "Pachyotoma crassicauda (Tullberg, 1871)")) %>% 
    separate(taxa, c("gen", "sp"), " ", extra = "drop") %>% 
    mutate(gen = paste0(substr(gen, 1, 1), ".")) %>% 
    unite(taxa, gen, sp, sep = " ") %>% 
    group_by(D, H, taxa) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    pivot_wider(names_from = taxa, values_from = abu) %>% 
    left_join(df.forcor, ., by = c("D", "H")) %>% 
    left_join(factors, by = c("D", "H")) %>% 
	select(-wind_2m, -temp_2m, -hum_abs, -log_lux)

cor.val2 <- cor(df.forcor[,12:15], df.forcor[,3:11], 
	method = "spearman")
cor.pval2 <- expand_grid(v1 = 12:15, v2 = 3:11) %>% 
	split(1:nrow(.)) %>% 
	lapply(function(a){
		p = cor.test(as_vector(df.forcor[,a$v1]), as_vector(df.forcor[,a$v2]), 
				   method = "spearman")
		data.frame(p = p$p.value, 
				 # est = p$estimate,
				 v1 = colnames(df.forcor[a$v1]),
				 v2 = colnames(df.forcor[,a$v2]))
	}) %>% 
	map_df(tibble) %>% 
	mutate(p = p.adjust(p, method = "BH")) %>% 
	pivot_wider(names_from = v2, values_from = p) %>% 
	column_to_rownames("v1") %>% 
	as.matrix()

rownames(cor.val2) <- 
    rownames(cor.pval2) <- 
    c("light level", "wind_20cm", "temperature_0cm", "humidity_relative")

pdf("Fig.8. Orders and factors_raw.pdf", height = 4)
corrplot::corrplot(
	corr = cbind(cbind(All = cor.val2[,1]), 
	             rep(NA, nrow(cor.val2)), 
	             cor.val2[,2:5], 
	             rep(NA, nrow(cor.val2)), 
	             cor.val2[,6:9]),
	p.mat = cbind(cbind(All = cor.pval2[,1]), 
	              rep(NA, nrow(cor.val2)), 
	              cor.pval2[,2:5], 
	              rep(NA, nrow(cor.val2)), 
	              cor.pval2[,6:9]), 
	type="full", 
	order = "original",
	na.label = " ",
	diag = TRUE,
	is.corr=FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off()
# Correlation of taxa’ activity (including dominant species) 
# and measured factors by Spearman method. 
# The Spearman coefficient value is shown on the scale on the right. 
# The significance level (α) = 0.05. 
# Statistically non-significant correlations are crossed out

