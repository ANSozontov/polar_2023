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

suppressWarnings(suppressMessages(
df <- readxl::read_excel("Data_28.04.2023.xlsx", 
			    sheet = "animals", skip = 3)
))
df <- df %>% 
	mutate_all(as.character) %>% 
	pivot_longer(names_to = "dt", values_to = "abu", -1:-2) %>% 
	separate(dt, sep = "_", into = c("trap_id", "D", "h")) %>% 
	separate(h , sep = "\\:", into = "H", extra = "drop") %>% 
	mutate_at(c("trap_id", "D", "H", "abu"), as.numeric)

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
# ggsave("Fig. A. total_abundance.pdf", width = 7, height = 4)

per2(abu$x, 1/6)

# Signal to noise ratio is  > 4
# Н0 about non-periodicity is rejected with p<0.001 significance level

data.frame(x = rep(1:6, 4), y = as.numeric(abu$x)) %>% 
	mutate(x = factor(x)) %>% 
	ggplot(aes(y = y, x = x, group = x, fill = x)) + 
	geom_boxplot() + 
	theme(legend.position = "none") + 
	labs(x = "Sampling hour", y = "Abundance/Activity") + 
	scale_x_discrete(breaks = 1:6, 
		labels = c("04", "08", "12", "16", "20","24")) + 
    scale_fill_manual(values = colorRampPalette(c("darkgrey", "orange"))
                      (6)[c(2, 4, 6, 5, 3, 1)])
# ggsave("Fig. B. Abundance_cycle.pdf", width = 3.5, height = 2)

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
ggsave("activity by hours - all taxa together.svg", height = 3, width = 3)

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
ggsave("activity by hours - by orders.svg", height = 6, width = 6)

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
ggsave("activity by hours - by dominant species.svg", height = 6, width = 6)


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
ggsave("Fig.C. taxa_abundance.pdf", width = 7, height = 4.5)

# dominant species  -------------------------------------------------------
dom_sp <- df %>%
	filter(taxa %in% c("Masikia indistincta",
		"Brachystomella parvula (Schäffer, 1896)", 
		"Isotomurus stuxbergi (Tullberg, 1876)", # check
		"Erigone psychrophila", # check
		"Isotomurus chaos Potapov et Babenko, 2011",
		"Pachyotoma crassicauda (Tullberg, 1871)")) %>% 
    separate(taxa, c("gen", "sp"), " ", extra = "drop") %>% 
	unite(taxa, gen, sp, sep = " ") %>% 
	split(.$taxa) %>% 
	map(~.x %>% 
	    	group_by(D, H) %>% 
	    	summarise(abu = sum(abu, na.rm = TRUE), .groups = "drop") %>% 
	        # mutate(abu = abu/max(abu)) %>% 
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
		# mutate_all(function(b){b/max(b, na.rm = TRUE)}) %>% 
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
						 levels = c("observed","trend", "periodic", "random")), 
            periodicity = case_when(taxa %in% c("Isotomurus stuxbergi", "Erigone psychrophila") ~ "no", 
                             TRUE ~ "yes"),
            taxa = factor(taxa, levels = c("Isotomurus stuxbergi", 
                "Brachystomella parvula", "Isotomurus chaos",
                "Masikia indistincta", "Pachyotoma crassicauda", 
                "Erigone psychrophila"))) %>% 
	ggplot(aes(x = id, y = val, color = taxa, shape = periodicity)) +
	geom_vline(xintercept = (1:4)*6-3, color = "red", linetype = "dotted") +
	geom_vline(xintercept = (1:4)*6, color = "blue", linetype = "dotted") +
	geom_line() + 
	geom_point(fill = "white") +
	theme_bw() +
	facet_grid(rows = vars(component), scales = "free") + 
	scale_x_continuous(breaks = 1:24, labels = rep(1:6*4, 4)) + 
    scale_shape_manual(values = c(20, 21)) +
	theme(axis.text.x = element_text(angle = 0), 
		 panel.grid.minor = element_blank(), 
		 panel.grid.minor.y = element_blank(), 
		 panel.grid.major.x = element_blank(), 
		 legend.text = element_text(face = "italic"),
		 legend.position = "bottom") +
	labs(x = "Sampling hours", y = "Abundance/activity and its components", 
		subtitle = "Note: all variables are normalized to max=1")
ggsave("Fig.7. dom.spec.pdf", width = 297, height = 210, units ="mm")
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

pdf("Fig. 7. Factor_multicollin.pdf", width = 4, height = 4.5)
corrplot::corrplot(
	corr = cor.val1, 
	p.mat = cor.pval1, 
	type="upper", 
	order = "original",
	diag = FALSE,
	col =  corrplot::COL2('RdYlBu', 10)[10:1], 
	sig.level = 0.05)
dev.off() 


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
	mutate(p = p.adjust(p, method = "BY")) %>% 
	pivot_wider(names_from = v2, values_from = p) %>% 
	column_to_rownames("v1") %>% 
	as.matrix()

rownames(cor.val2) <- 
    rownames(cor.pval2) <- 
    c("light level", "wind_20cm", "temperature_0cm", "humidity_relative")

pdf("Fig. 8. Orders_and_factors.pdf", height = 4)
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
