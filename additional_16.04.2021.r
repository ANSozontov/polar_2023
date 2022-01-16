# data.load ---------------------------------------------------------------
library(tidyverse) #
library(iNEXT)
nspec <- function(d, go = 1, finish = ncol(d)){ 
  apply(d[,go:finish], 1, FUN = function(a){length(a[a>0])})
}
extr <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
wide.df <- readxl::read_xlsx("Mukhacheva_data_18022021.xlsx") %>% 
  filter(field == 1) %>%  # there were no samples at 6&30 before 2004
  select(-field)

# spatial beta along the time ------------------------------------------------
B <- wide.df %>% 
  mutate(period = apply(wide.df[,8:21], 1, function(a){length(a[a>0])})) %>% 
  # period is used as n.spec here
  filter(period > 2) %>% 
  mutate(rrow = paste0("v_", year, ".", km)) %>% 
  select(8:ncol(.)) %>% 
  column_to_rownames("rrow") %>% 
  t %>% 
  as.tibble() %>% 
  as.list() %>% 
  map(., .f = function(a){a[a>0]})
j <- 1:110
j <- j[!(j %in% c(6, 16, 29, 30, 36, 49, 64, 65, 79, 82, 88, 95, 97, 99, 100))]
a.rar <- data.frame(v = character(), m = numeric(), qD = numeric())
for(i in j) {
  a.rar <- rbind(a.rar, 
                 iNEXT(B[[i]], q = 0, datatype = "abundance", 
                       size = c(5, 10, 15, 25))$iNextEst %>% 
                   filter(m %in% c(5, 10, 15, 25)) %>% 
                   transmute(v = names(B)[i], m, qD)
  )
  }

n <- a.rar %>% #transmute(v = substr(v, 3, 6)) %>% 
  mutate(v = substr(v, 1, 6)) %>% 
  group_by(v) %>% 
  summarise(nn = n()/4) %>%  #cont., sample size 
  mutate(nn = case_when(nn == 1 ~ 1.5, TRUE ~ nn), 
         n5 = nn*5, n10 = nn*10, n15 = nn*15, n25 = nn*25)

G <- wide.df %>% 
  mutate(year = paste0("v_", year), 
         period = apply(wide.df[,8:21], 1, function(a){length(a[a>0])})) %>% 
  # period is used as n.spec here
  filter(period > 2) %>% 
  select(year, 8:21) %>% 
  group_by(year) %>% 
  summarise_all(mean) %>% 
  as.data.frame()
rownames(G) <- G$year
G <- G[,2:ncol(G)]
G <- G %>% 
  t %>% 
  as.tibble() %>% 
  as.list() %>% 
  map(., .f = function(a){a[a>0]})
j <- 1:27
j <- j[!(j %in% c(5, 10, 19, 25))]
g.rar <- data.frame(v = character(), m = numeric(), qG = numeric())
for(i in j) {
  abu <- as.vector(as.matrix(n[which(n$v == names(G)[15]),3:6]))
  g.rar <- rbind(g.rar, 
                 iNEXT(G[[i]], q = 0, datatype = "abundance", 
                       size = abu)$iNextEst %>% 
                   filter(m %in% abu) %>% 
                   transmute(v = names(G)[i], realm = m, m = c(5, 10, 15, 25), qG = qD)
  )
}

rar <- a.rar %>% #rarefied
  mutate(v = substr(v, 1, 6)) %>% 
  full_join(g.rar, by = c("v", "m")) %>% 
  as_tibble() %>% 
  filter(!is.na(qG)) %>% 
  mutate(v = as.numeric(substr(v, 3, 6))) %>% 
  group_by(v, m) %>% 
  summarise(b = mean(qG) / mean(qD), .groups = "drop") %>% 
  rename(year = v)

bw <- wide.df %>% #Bw as is 
  select(-id, -period, -cycle, -zone, -total) %>% 
  group_by(year, km) %>% 
  summarise_all(mean) %>% 
  ungroup %>% 
  mutate(ss = nspec(., 3, ncol(.))) %>% 
  group_by(year) %>%
  summarise_all(mean) %>% 
  transmute(year, m = "bw.as.is.", b = nspec(., 3, ncol(.) - 1) / ss)

ibc.years <- as.matrix(vegan::vegdist(wide.df[,8:ncol(wide.df)], #ibc as is
    method = "bray", binary = FALSE))
rownames(ibc.years) <- colnames(ibc.years) <- wide.df$id
ibc <- reshape2::melt(ibc.years, varnames = c("id1", "id2")) %>% 
  as_tibble() %>% 
  mutate(id1 = as.character(id1), id2 = as.character(id2)) %>% 
  filter(as.numeric(substr(id1, 4, nchar(id1))) > as.numeric(substr(id2, 4, nchar(id2)))) %>% 
  left_join(select(wide.df, id1 = id, km1 = km, year1 = year, zone1 = zone), by = "id1") %>% 
  left_join(select(wide.df, id2 = id, km2 = km, year2 = year, zone2 = zone), by = "id2") %>% 
  filter(year1 == year2) %>% 
  mutate(value = case_when(is.nan(value) ~ 1, TRUE ~ value)) %>% 
  group_by(year1) %>% 
  summarise(B = mean(value), .groups = "drop") %>% 
  transmute(year = year1, m = "iBC.as.is", b = B)

ibc.years <- as.matrix(vegan::vegdist(wide.df[,8:ncol(wide.df)], #ics
                                      method = "bray", binary = TRUE))
rownames(ibc.years) <- colnames(ibc.years) <- wide.df$id
ics <- reshape2::melt(ibc.years, varnames = c("id1", "id2")) %>% 
  as_tibble() %>% 
  mutate(id1 = as.character(id1), id2 = as.character(id2)) %>% 
  filter(as.numeric(substr(id1, 4, nchar(id1))) > as.numeric(substr(id2, 4, nchar(id2)))) %>% 
  left_join(select(wide.df, id1 = id, km1 = km, year1 = year, zone1 = zone), by = "id1") %>% 
  left_join(select(wide.df, id2 = id, km2 = km, year2 = year, zone2 = zone), by = "id2") %>% 
  filter(year1 == year2) %>% 
  mutate(value = case_when(is.nan(value) ~ 1, TRUE ~ value)) %>% 
  group_by(year1) %>% 
  summarise(B = mean(value), .groups = "drop") %>% 
  transmute(year = year1, m = "iCz.Soe", b = B)

rm(a.rar, g.rar, B, G, n, i, j, abu, ibc.years)

# vis ----
spat <- readxl::read_xlsx("Mukhacheva_data_19042021.xlsx", sheet = "beta_in_space") %>% 
  pivot_longer(names_to = "m", values_to = "b", -year) %>% 
  rbind(bw, rar, ibc, ics) %>% 
  filter(!is.na(b), nchar(m) > 2, 
         !(m %in% c("Cody", "Williams", "Wilson.Shmida", "bw.as.is."))) %>% 
  mutate(m = case_when(m == "Routledge" ~ "a. Routledge", m == "Harrison" ~ "c. Harrison", 
                       m ==  "Harrison_2" ~ "b. Harrison_2", m == "Whittaker" ~ "d. Whittaker (old)", 
                       m == "Mourelle" ~ "e. Mourelle", m == "iBC.as.is" ~ "f. Bray-Curt. (old)", 
                       m == "iCz.Soe" ~ "g. Czek._Soer."
  ))

add0 <- ggplot(spat, aes(x = year, y = b)) + 
  geom_line(color = "black") + 
  geom_point(size = 2.5, color = "black", fill = "darkgrey", shape = 21) + 
  facet_grid(rows = vars(m), scales = "free") +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90), 
        panel.grid.minor.x = element_blank()) + 
  scale_x_continuous(breaks = 1990:2020) + 
  geom_vline(aes(xintercept = d), size = 1, linetype = "dashed", 
             data = data.frame(d = c(1997.5, 2009.5))) +
  labs(x = "", y = "", title = "")
add0
ggsave("add0.png", add0, width = 8, height = 10, dpi = 1200)

d <- spat %>% pivot_wider(names_from = m, values_from = b) %>% 
  select(-year) 
e <- expand.grid(d1 = names(d), d2 = names(d)) %>% 
  mutate(sim = NA) %>% as_tibble()
for(i in 1:nrow(e)){ 
  f <- tibble(
    x = pull(d[,which(colnames(d) == e$d1[i])], 1), 
    y = pull(d[,which(colnames(d) == e$d2[i])], 1)
    ) %>%
    filter(!is.na(x), !is.na(y))
  e$sim[i] <- round(cor(f$x, f$y), 3)
  #e$sim[i] <- round(cor.test(f$x, f$y)$p.value, 4)
}
e %>% pivot_wider(names_from = d2, values_from = sim)
e %>% mutate(sim = round(sim, 2)) %>% 
  pivot_wider(names_from = d2, values_from = sim) %>% 
  writexl::write_xlsx("temp.xlsx")

rm(d, f, e, ibc, ics, bw, rar)
# A tibble: 7 x 8
# d1         bbw.as.is.     b5    b10    b15    b25 biBC.as.is biCz.Soe
# <fct>           <dbl>  <dbl>  <dbl>  <dbl>  <dbl>      <dbl>    <dbl>
# 1 bbw.as.is.      1     -0.463 -0.303 -0.083 -0.006      0.756    0.804
# 2 b5             -0.463  1      0.805  0.572  0.304     -0.352   -0.563
# 3 b10            -0.303  0.805  1      0.849  0.694     -0.195   -0.533
# 4 b15            -0.083  0.572  0.849  1      0.921     -0.226   -0.4  
# 5 b25            -0.006  0.304  0.694  0.921  1         -0.142   -0.286
# 6 biBC.as.is      0.756 -0.352 -0.195 -0.226 -0.142      1        0.739
# 7 biCz.Soe        0.804 -0.563 -0.533 -0.4   -0.286      0.739    1    

# temporal beta along the space -------------------------------------------
bw <- wide.df  %>%  # bw.as.is 
  select(-id, -cycle, -zone, -total) %>% 
  mutate(nsp = nspec(., 3, ncol(.))) %>% 
  group_by(period, km) %>% 
  summarise_all(mean) %>% 
  ungroup() %>% 
  transmute(per = period, km, type = "bw.old", b = nspec(., 3, ncol(.) - 1) / nsp)

ibc.km <- as.matrix( 
  vegan::vegdist(wide.df[,8:ncol(wide.df)], method = "bray", binary = FALSE))
rownames(ibc.km) <- colnames(ibc.km) <- wide.df$id
ibc <- reshape2::melt(ibc.km, varnames = c("id1", "id2")) %>% # bray-curtis
  as_tibble() %>% 
  mutate(id1 = as.character(id1), id2 = as.character(id2)) %>% 
  filter(as.numeric(substr(id1, 4, nchar(id1))) > as.numeric(substr(id2, 4, nchar(id2)))) %>% 
  left_join(select(wide.df, id1 = id, km1 = km, per1 = period), by = "id1") %>% 
  left_join(select(wide.df, id2 = id, km2 = km, per2 = period), by = "id2") %>% 
  filter(km1 == km2, per1 == per2) %>% 
  mutate(value = case_when(is.nan(value) ~ 1, TRUE ~ value)) %>% 
  group_by(per1, km1) %>% 
  summarise(ibc = mean(value), .groups = "drop") %>% 
  transmute(per = per1, km = km1, type = "bray-curtis", b = ibc)

ibc.km <- as.matrix( 
  vegan::vegdist(wide.df[,8:ncol(wide.df)], method = "bray", binary = TRUE))
rownames(ibc.km) <- colnames(ibc.km) <- wide.df$id
ics <- reshape2::melt(ibc.km, varnames = c("id1", "id2")) %>% # Cze-Soe
  as_tibble() %>% 
  mutate(id1 = as.character(id1), id2 = as.character(id2)) %>% 
  filter(as.numeric(substr(id1, 4, nchar(id1))) > as.numeric(substr(id2, 4, nchar(id2)))) %>% 
  left_join(select(wide.df, id1 = id, km1 = km, per1 = period), by = "id1") %>% 
  left_join(select(wide.df, id2 = id, km2 = km, per2 = period), by = "id2") %>% 
  filter(km1 == km2, per1 == per2) %>% 
  mutate(value = case_when(is.nan(value) ~ 1, TRUE ~ value)) %>% 
  group_by(per1, km1) %>% 
  summarise(ibc = mean(value), .groups = "drop") %>% 
  transmute(per = per1, km = km1, type = "soerensen", b = ibc)

## rarefy ------------------------------------------------------------------
A <- wide.df %>% 
  mutate(id = paste0(period, "_", year, km)) %>% 
  select(-year, -period, -cycle, -zone, -km, -total) %>% 
  column_to_rownames("id") %>% 
  t %>% 
  as_tibble() %>% 
  as.list() %>% 
  map(., .f = function(a){a[a>0]}) %>% 
  keep(~ length(.x) > 1)
a <- data.frame(v = character(), m = numeric(), a = numeric())
j <- 1:length(A)
j <- j[!(j %in% c(6, 15, 19, 27, 37:39, 45, 47, 60, 62, 79, 82, 94, 98, 101, 
                  108, 115, 116, 118, 119, 121, 122))]
for(i in j) {
  a <- rbind(a, 
             iNEXT(A[[i]], q = 0, datatype = "abundance", 
                   size = c(25, 35, 45))$iNextEst %>% 
               filter(m %in% c(25, 35, 45)) %>% 
               transmute(v = names(A)[i], m, a = qD)
  )
} # 40, 70, 110

G <- wide.df %>% 
  select(-id, -year, -cycle, -zone, -total) %>% 
  group_by(period, km) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
  mutate(vv = paste0(period, ".", km)) %>% 
  select(-period, -km) %>% 
  column_to_rownames("vv") %>% 
  t %>% 
  as_tibble() %>% 
  as.list() %>% 
  map(., .f = function(a){a[a>0]}) %>% 
  keep(~ length(.x) > 1)
j <- 1:length(G)

g <- data.frame(v = character(), n = numeric(), g = numeric())
for(i in j) {
  g <- rbind(g, 
             iNEXT(G[[i]], q = 0, datatype = "abundance", 
                   size = c(100, 125, 150))$iNextEst %>% 
               filter(m %in% c(100, 125, 150)) %>% 
               transmute(v = names(G)[i], n = m, g = qD)
  )
}

bb <- full_join( #rarefied
  transmute(a, per = substr(v, 1,1), km = substr(v, 7, nchar(v)), m, a), 
  transmute(g, per = substr(v, 1,1), km = substr(v, 3, nchar(v)), n, g), 
  by = c("per", "km")
) %>% as_tibble() %>% 
  group_by(per, km, n, m) %>% 
  summarise(b = mean(g)/mean(a), .groups = "drop") %>% 
  #filter( (n == 40 & m == 5) | (n == 70 & m == 10) | (n == 110 & m == 25)) %>% 
  transmute(per = as.numeric(per), km = as.numeric(km), 
            type = paste0("R. a= ", m, "; b= ", n), b)
#readxl::read_xlsx("Mukhacheva_data_19042021.xlsx", sheet = "beta_in_time")
rm(a, A, g, G, ibc.km, p, rar, i, j, nspec)
tm <- readxl::read_xlsx("Mukhacheva_data_19042021.xlsx", sheet = "beta_in_time") %>% 
  pivot_longer(names_to = "type", values_to = "b", -c("period", "km")) %>% 
  rename(per = period) %>% 
  rbind(., bw, ibc, ics, bb) %>%  
  filter(!is.na(b)) %>% 
  filter(type %in% c("bw.old", "bw.old", "Harrison_2", "Harrison", "Routledge", "soerensen", "bray-curtis")) %>% 
  mutate(type = case_when(type == "bw.old" ~ "a. Whittaker (old)", type == "Routledge" ~ "b. Routledge", 
                          type == "Harrison_2" ~ "c. Harrison_2", type == "bray-curtis" ~ "d. Bray-Curtis (old)", 
                          type == "soerensen" ~ "e. Czek.-Soer.", type ==  "Harrison" ~ "f. Harrison_1"))


# models & vis -------------------------------------------------------------------
res <- expand.grid(per = 1:3, type = unique(tm$type), 
                   int = NA, sl = NA, aic = NA, pval = NA, r2 = NA, r2a = NA)
for(i in 1:nrow(res)) { 
  fit <- lm(b ~ log(km), data = filter(tm, type == res$type[i], per == res$per[i]))
  res[i, 3:4] <- fit$coefficients
  res$aic[i] <- round(AIC(fit), 1)
  res$pval[i] <- round(extr(fit), 4)
  fit <- summary(fit)
  res$r2[i] <- round(fit$r.squared, 2)
  res$r2a[i] <- round(fit$adj.r.squared, 2)
  rm(fit)
  rm(i)
}

add1 <- ggplot() +
  geom_point(aes(x = log(km), y = b, fill = period, shape = period), size = 2,
             data = mutate(tm, period = as.factor(per))) + 
  geom_abline(aes(intercept = int, slope = sl, color = period), 
              data = mutate(res, period = as.factor(per))) +
  facet_wrap(~type, scales = "free") +
  scale_shape_manual(values = c(21, 22, 24)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "")
add1
ggsave("add1.png", add1, width = 7, height = 5, dpi = 1200)







