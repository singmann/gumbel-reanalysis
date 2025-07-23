warmup <- 1000
iter <- 6000
source("fit-6point.R")

#iter <- 4000
source("fit-8point.R")

##### Rhat
source("check-functions.R")
max(vapply(roc6_fits_gumbel, get_max_rhat, 0))
## 1.006571

which.max(vapply(roc6_fits_gumbel, get_max_rhat, 0))
max(vapply(roc8_fits_gumbel, get_max_rhat, 0))
# [1] 1.006891

max(vapply(roc6_fits_uvsdt, get_max_rhat, 0))
## 1.005455
max(vapply(roc8_fits_uvsdt, get_max_rhat, 0))
## 1.008181


#### divergent transitions
#iter <- 2000
xxx <- map(roc6_fits_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset6)) {
  cat(dataset6[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[(warmup+1):dim(.)[1],"divergent__"])))/((iter-warmup)*4), "\n")
}

xxy <- map(roc6_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset6)) {
  cat(dataset6[i], ": ", sum(map_dbl(xxy[[i]], ~sum(.[(warmup+1):dim(.)[1],"divergent__"])))/((iter-warmup)*4), "\n")
}

xxx <- map(roc8_fits_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset8)) {
  cat(dataset8[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[(warmup+1):dim(.)[1],"divergent__"])))/((iter-warmup)*4), "\n")
}

xxy <- map(roc8_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset8)) {
  cat(dataset8[i], ": ", sum(map_dbl(xxy[[i]], ~sum(.[(warmup+1):dim(.)[1],"divergent__"])))/((iter-warmup)*4), "\n")
}

## plots
plot_data6 <- roc6 %>% 
  group_by(exp) %>% 
  summarise(across(c(OLD_3new:NEW_3old), sum)) %>% 
  pivot_longer(-exp, names_to = c("status", "response"), names_sep = "_") %>% 
  group_by(exp, status) %>% 
  mutate(observed = value / sum(value)) %>% 
  mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("3new", "2new", "1new", 
                                           "1old", "2old", "3old")))

pred_gumbel <- lapply(roc6_fits_gumbel, posterior_epred)
pred_uvsdt <- lapply(roc6_fits_uvsdt, posterior_epred)

plot_data6$gumbel <- unlist(map(pred_gumbel, ~apply(., c(3), mean)))
plot_data6$uvsd <- unlist(map(pred_uvsdt, ~apply(., c(3), mean)))

plot_data8 <- roc8 %>% 
  #filter(exp != dataset8[1]) %>% 
  group_by(exp) %>% 
  summarise(across(c(OLD_4new:NEW_4old), sum)) %>% 
  pivot_longer(-exp, names_to = c("status", "response"), names_sep = "_") %>% 
  group_by(exp, status) %>% 
  mutate(observed = value / sum(value)) %>% 
  mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("4new", "3new", "2new", "1new", 
                                           "1old", "2old", "3old", "4old")))

pred_gumbel <- lapply(roc8_fits_gumbel, posterior_epred)
pred_uvsdt <- lapply(roc8_fits_uvsdt, posterior_epred)

plot_data8$gumbel <- unlist(map(pred_gumbel, ~apply(., c(3), mean)))
plot_data8$uvsd <- unlist(map(pred_uvsdt, ~apply(., c(3), mean)))
# 
# plot_dat <- bind_rows(plot_data6, plot_data8)
# plot_dat <- plot_dat %>% 
#   mutate(newexp = factor(
#     exp, 
#     levels = 
#       c("Dube_2012-P", "Dube_2012-W", "Heathcote_2006_e1", 
#         "Heathcote_2006_e2", "Jaeger_2012", "Jang_2009", "Koen_2010_pure", 
#         "Koen_2011", "Koen-2013_full", "Koen-2013_immediate", "Pratte_2010", 
#         "Smith_2004", "Benjamin_2013", "Onyper_2010-Pics", "Onyper_2010-Words"), 
#     labels = 
#       c("Dube (2012, E1b-P)", "Dube (2012 E1b-W)", "Heathcote (2006, E1)", 
#         "Heathcote (2006, E2)", "Jaeger (2012, E1)", "Jang (2009)", "Koen (2010)", 
#         "Koen (2011)", "Koen (2013, E2)", "Koen (2013, E4)", "Pratte (2010)", 
#         "Smith (2004)", "Benjamin (2013)", "Onyper (2010, E1-P)", "Onyper (2010, E1-W"))) %>% 
#   ungroup()


library(showtext)
font_paths("fonts")

# Add font
font_add("Palatino Linotype", 
         #regular="pala.ttf", 
         regular = "asana-math.otf",
         #regular = "palatinolinotype_roman.ttf",
         italic = "palatinolinotype_italic.ttf", 
         bold = "palatinolinotype_bold.ttf", 
         bolditalic = "palatinolinotype_bolditalic.ttf")
showtext_auto()

theme_set(theme_bw(base_size = 12, base_family = "Palatino Linotype") + 
            theme(legend.position="bottom", 
                  panel.grid = element_blank()))



plot_data8_2 <- plot_data8 %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "4new") 

plot_data6_2 <-plot_data6 %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "3new")

plot_data_2 <- bind_rows(plot_data6_2, plot_data8_2) %>% 
  mutate(newexp = factor(
    exp, 
    levels = 
      c("Dube_2012-P", "Dube_2012-W", "Heathcote_2006_e1", 
        "Heathcote_2006_e2", "Jaeger_2012", "Jang_2009", "Koen_2010_pure", 
        "Koen_2011", "Koen-2013_full", "Koen-2013_immediate", "Pratte_2010", 
        "Smith_2004", "Benjamin_2013", "Onyper_2010-Pics", "Onyper_2010-Words"), 
    labels = 
      c("Dube (2012, E1b-P)", "Dube (2012 E1b-W)", "Heathcote (2006, E1)", 
        "Heathcote (2006, E2)", "Jaeger (2012, E1)", "Jang (2009)", "Koen (2010)", 
        "Koen (2011)", "Koen (2013, E2)", "Koen (2013, E4)", "Pratte (2010)", 
        "Smith (2004)", "Benjamin (2013)", "Onyper (2010, E1-P)", "Onyper (2010, E1-W)"))) %>% 
  ungroup()

bin_n <- bind_rows(roc6, roc8) %>% 
  group_by(exp) %>% 
  summarise(n = n_distinct(id)) %>% 
  mutate(n_text = paste0("italic(N) == ", n))
bin_n <- left_join(bin_n, unique(select(plot_data_2, exp, newexp)))

psize <- 3.5
lsize <- 1.5
ssize <- 1.0
plot_data_2 %>%
  ggplot(aes(x =  observed_NEW, y = observed_OLD)) +
  geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "white") +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x = gumbel_NEW, y = gumbel_OLD, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize, stroke = ssize) +
  geom_point(aes(x = uvsd_NEW, y = uvsd_OLD, 
                 shape = "UVSD", colour = "UVSD"), size = psize, stroke = ssize) + 
  geom_label(mapping = aes(x = 0.75, y = 0.15, label = n_text), 
             data = bin_n, hjust = "center", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(newexp), nrow = 3) + 
  scale_color_manual(
    name = '',
    breaks = c('Data', 'UVSD', 'Gumbel'),
    values = c('Data' = 'black', 'UVSD' = "#0072B2", 'Gumbel' = "#E69F00"),
    labels = c("Data", "Gaussian", expression(Gumbel[min]))
  ) +
  scale_shape_manual(
    name = '',
    breaks = c('Data', 'UVSD', 'Gumbel'),
    values = c('Data' = 19, 'Gumbel' = 3, 'UVSD' = 5),
    labels = c("Data", "Gaussian", expression(Gumbel[min]))
  ) + 
  theme(legend.title = NULL) +
  labs(x = expression(italic(p)[FA]), y = expression(italic(p)[H]))
ggsave("confroc-plot1.pdf", width = 22, height = 17.25, units = "cm")


