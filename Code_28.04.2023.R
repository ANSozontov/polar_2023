library(tidyverse)
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
loggers <- readxl::read_excel("Data_28.04.2023.xlsx", 
		sheet = "logger_data", skip = 1) %>% 
	# select(D, H) %>%
	mutate(D = as.numeric(substr(D, 1, 2)), 
		  H = (floor(H/4)+1)*4) %>%
	group_by(D, H) %>% 
	summarise_all(mean) %>% 
	ungroup() %>% 
	filter(D != 9, D != 18)

factors <- readxl::read_excel("Data_28.04.2023.xlsx", 
		sheet = "factors", range = "A2:E26") %>% 
	mutate(D = substr(str_squish(D), 1, 2), 
		  D = as.numeric(D), 
		  log_lux = log10(lux), .after = 2)
factors <- left_join(factors, loggers, by = c("D", "H")) 

suppressMessages(suppressWarnings(
df <- readxl::read_excel("Data_28.04.2023.xlsx", 
			    sheet = "animals", skip = 3)
))
df <- df %>% 
	mutate_all(as.character) %>% 
	pivot_longer(names_to = "dt", values_to = "abu", -1:-2) %>% 
	separate(dt, sep = "_", into = c("D", "h")) %>% 
	separate(h , sep = "\\:", into = "H", extra = "drop") %>% 
	mutate_at(c("D", "H", "abu"), as.numeric) #%>% filter(abu > 0)

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
	theme_bw() +
	facet_grid(rows = vars(component)) + #, scales = "free") + 
	scale_x_continuous(breaks = 1:24, 
		labels = rep(1:6*4, 4)) + 
	theme(axis.text.x = element_text(angle = 0), 
		 panel.grid.minor = element_blank(), 
		 panel.grid.minor.y = element_blank(), 
		 panel.grid.major.x = element_blank()) +
	labs(x = "Sampling hours", y = "Abundance and its components")
ggsave("total_abundance.svg")

per2(abu$x, 1/6)

# соотношение сигнал/шум > 4
# Н0 о непериодичности активности отклонена на уровне p<0.001

data.frame(x = rep(1:6, 4), y = as.numeric(abu$x)) %>% 
	mutate(x = factor(x)) %>% 
	ggplot(aes(y = y, x = x, group = x, fill = x)) + 
	geom_boxplot() + 
	theme_bw() + 
	theme(legend.position = "none") + 
	labs(x = NULL, y = "Abundance/Activity") + 
	scale_x_discrete(breaks = 1:6, 
		labels = c("04", "08", "12", "16", "20","24"))

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

by_taxa %>% lapply(function(a){
		cbind(x = a$x, trend = a$trend, 
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
	# geom_vline(xintercept = (1:4)*6-3) +
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
	labs(x = "Sampling hours", y = "Abundance and its components", 
		subtitle = "all variables are normalized to max = 1")

# dominant species  -------------------------------------------------------
dom_sp <- df %>%
	filter(taxa %in% c("Erigone psychrophila", "Masikia indistincta",
		"Brachystomella parvula (Schäffer, 1896)", 
		"Isotomurus chaos Potapov et Babenko, 2011",
		"Isotomurus stuxbergi (Tullberg, 1876)", 
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
		 legend.text = element_text(face = "italic"),
		 legend.position = "bottom") +
	labs(x = "Sampling hours", y = "Abundance and its components", 
		subtitle = "all variables are normalized to max = 1")

# I.stuxbergi
I.stuxbergi <- dom_sp %>% 
	pluck("Isotomurus stuxbergi", "x") %>% 
	as.numeric()

per2(I.stuxbergi, 1/2)
per2(I.stuxbergi, 1/3)
per2(I.stuxbergi, 1/4)
per2(I.stuxbergi, 1/5)
per2(I.stuxbergi, 1/6)
per2(I.stuxbergi, 1/7)
per2(I.stuxbergi, 1/8)
per2(I.stuxbergi, 1/9)

# predictors multicoll & periodicity ----------------------------------------------------


df.forgam <- df %>% 
	filter(taxatype == "general") %>% 
	group_by(D, H) %>% 
	summarise(abu = sum(abu), .groups = "drop") %>% 
	left_join(factors, by = c("D", "H")) %>% 
	mutate(decompress(abu), decompress(lux), decompress(log_lux))

factors[,-1:-2] %>% 
	cor(method = "spearman") %>% 
	round(2)

cor.val <- cor(factors[,-1:-3], method = "spearman")
cor.pval <- expand_grid(v1 = 4:10, v2 = 4:10) %>% 
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


svg("multicollin.svg")
corrplot::corrplot(
	corr = cor.val, 
	p.mat = cor.pval, 
	type="upper", 
	order = "original",
	diag = FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off() 

df.forcor <- df %>% 
	filter(taxatype == "general", taxa %in% c("Collembola", "Aranei", "Acari", "Diptera")) %>% 
	group_by(D, H, taxa) %>% 
	summarise(abu = sum(abu), .groups = "drop") %>% 
	pivot_wider(names_from = taxa, values_from = abu) %>% 
	left_join(df.forgam, by = c("D", "H")) %>% 
	select(All = abu, Collembola, Acari, Aranei, Diptera, 
		  lux, wind_20cm, temp_20cm, hum_rel)

cor.val <- cor(df.forcor[,1:5], df.forcor[,6:9], 
	method = "spearman")
cor.pval <- expand_grid(v1 = 1:5, v2 = 6:9) %>% 
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
	mutate(p = p.adjust(p, method = "BY")) %>% 
	pivot_wider(names_from = v2, values_from = p) %>% 
	column_to_rownames("v1") %>% 
	as.matrix()

svg("orders_and_factors.svg")
corrplot::corrplot(
	corr = cor.val, 
	p.mat = cor.pval, 
	type="full", 
	order = "original",
	diag = TRUE,
	is.corr=FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off()


df.mass <- df %>% 
	filter(taxa %in% c("Erigone psychrophila", "Masikia indistincta",
				    "Brachystomella parvula (Schäffer, 1896)", 
				    "Isotomurus chaos Potapov et Babenko, 2011",
				    "Isotomurus stuxbergi (Tullberg, 1876)", 
				    "Pachyotoma crassicauda (Tullberg, 1871)")) %>% 
	separate(taxa, c("gen", "sp"), " ", extra = "drop") %>% 
	unite(taxa, gen, sp, sep = " ") %>% 
	group_by(D, H, taxa) %>% 
	summarise(abu = sum(abu), .groups = "drop") %>% 
	pivot_wider(names_from = taxa, values_from = abu) %>% 
	left_join(df.forgam, by = c("D", "H")) %>% 
	select(3:8, lux, wind_20cm, temp_20cm, hum_rel)

cor.val <- cor(df.mass[,1:6], df.mass[,7:10], 
			method = "spearman")
cor.pval <- expand_grid(v1 = 1:6, v2 = 7:10) %>% 
	split(1:nrow(.)) %>% 
	lapply(function(a){
		p = cor.test(as_vector(df.mass[,a$v1]), as_vector(df.mass[,a$v2]), 
				   method = "spearman")
		data.frame(p = p$p.value, 
				 # est = p$estimate,
				 v1 = colnames(df.mass[a$v1]),
				 v2 = colnames(df.mass[,a$v2]))
	}) %>% 
	map_df(tibble) %>% 
	mutate(p = p.adjust(p, method = "BY")) %>% 
	pivot_wider(names_from = v2, values_from = p) %>% 
	column_to_rownames("v1") %>% 
	as.matrix()

svg("dominants_and_factors.svg")
corrplot::corrplot(
	corr = cor.val, 
	p.mat = cor.pval, 
	type="full", 
	order = "original",
	diag = TRUE,
	is.corr=FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off()

# models ------------------------------------------------------------------
library(mgcv)


# DF <- abu[1:4] %>% 
# 	as.data.frame() %>% 
# 	rename_with(~ paste0("abu_", .x)) %>% 
# 	cbind(DF, .)
	


mgcv::gam(abu ~ log_lux + hum_rel, data = DF) %>% ft
mgcv::gam(abu ~ s(log_lux) + hum_rel, data = DF) %>% ft
mgcv::gam(abu ~ s(log_lux) + s(hum_rel), data = DF) %>% ft
mgcv::gam(abu ~ log_lux + hum_rel + wind_20cm, data = DF) %>% ft
mgcv::gam(abu ~ log_lux + hum_rel + s(D, bs = 're') , data = DF) %>% ft
mgcv::gam(abu ~ s(log_lux) + hum_rel + s(D, bs = 're') , data = DF) %>% ft
mgcv::gam(abu ~ s(lux)     + hum_rel + s(D, bs = 're') , data = DF) %>% ft
mgcv::gam(abu_seasonal ~ lux_seasonal     + hum_rel + s(D, bs = 're') , data = DF) %>% ft


mgcv::gam(abu_random ~ lux + hum_rel , data = DF) %>% ft
mgcv::gam(abu_random ~ hum_rel , data = DF) %>% ft
mgcv::gam(abu_random ~ lux + hum_rel , data = DF) %>% ft


fit <- mgcv::gam(abu ~ s(lux, k = 10)     + hum_rel + s(D,  bs = 're') , data = DF) 
mgcv::gam.check(fit)



fit3 <- lm(abu ~ lux_raw + temp_20cm + hum_abs, data = DF)
summary(fit3)
AIC(fit3)
