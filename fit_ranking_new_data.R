library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
library("patchwork")
source("check-functions.R")

source("data_from_david.R")

source("gumbelrank-stan.R")
source("uvsdtranknew-stan.R")

gumbel_priors <- prior(student_t(3, 1, 2), class = "b", coef = "Intercept")
#gumbel_priors <- prior(student_t(3, 1, 2), class = Intercept)
uvsdt_priors <- prior(student_t(3, 0.5, 1), 
                      class = "Intercept", dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = "b", coef = "Intercept")

# xxx <- 1/exp(extraDistr::rlst(1e4, 20, 0.5, 2))
# mean(xxx > 20)
# xxx <- xxx[xxx < 20]
# plot(density(xxx))


##----------------------------------------------------------------
##                              Data                             -
##----------------------------------------------------------------

kks12l <- kks12 %>% 
  mutate(exp = "Kellen (2012)") %>% 
  pivot_longer(cols = -c(exp, id), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 4L) %>% 
  mutate(strength = "r") %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

kk14_e1_use <- kk14_e1 %>% 
  mutate(exp = "Kellen (2014, E1)") %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s"))) %>% 
  pivot_longer(cols = -c(exp, id, strength), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 4L) %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

mhe_e1_use <- mhe_e1 %>% 
  mutate(exp = "Malejka (2022, E1)") %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s"))) %>% 
  pivot_longer(cols = -c(exp, id, strength), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 4L) %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

kk14_e2_use <- kk14_e2 %>% 
  mutate(exp = "Kellen (2014, E2)") %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s"))) %>% 
  pivot_longer(cols = -c(exp, id, strength), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 3L) %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

mg16_e1
mg16_e1_use <- mg16_e1 %>% 
  mutate(exp = "McAdoo (2016, E1)") %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))  %>% 
  pivot_longer(cols = -c(exp, id, strength), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 3L) %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

mg16_e2
mg16_e2_use <- mg16_e2 %>% 
  mutate(exp = "McAdoo (2016, E2)") %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s"))) %>% 
  pivot_longer(cols = -c(exp, id, strength), 
               names_to = "rank", values_to = "observed") %>% 
  mutate(rank = as.integer(as.numeric(str_extract(rank, "\\d"))),
         maxrank = 3L) %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

dmgj24 <- read_csv("data_mj2024.csv")
mgj24 <- dmgj24 %>% 
  rename(id = ID) %>% 
  mutate(exp = "Meyer-Grant (2024)") %>% 
  group_by(exp, id, n_images, rank_target) %>% 
  summarise(n=n(), .groups="drop") %>% 
  rename(maxrank = n_images, rank = rank_target, observed = n) %>% 
  mutate(strength = "r") %>% 
  select(exp, id, strength, rank, observed, maxrank) %>% 
  mutate(id = as.character(id))

##---------------------------------------------------------------
##                            4 Ranks                           -
##---------------------------------------------------------------

##------------
##  KKS 2012  
##------------

gumbel_formula_kks <- brmsformula(
  observed | vint(rank, maxrank) ~ 1 + (1|p|id), 
  family = gumbelrank_family, center = FALSE
)

# stancode(gumbel_formula_kks, data = kks12, 
#          stanvars = sv_gumbelrank, prior = gumbel_priors)

fit_kks_gumbel <- brm(
    gumbel_formula_kks, data = kks12l, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kks <- brmsformula(
  observed | vint(rank, maxrank) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, center = FALSE
)

# stancode(uvsdt_formula_kks, data = kks12, 
#          stanvars = sv_uvsdtrank, prior = uvsdt_priors)

fit_kks_uvsdt <- brm(
    uvsdt_formula_kks, data = kks12l, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.1
  )

kks12


##-----------
##  KK14 E1  
##-----------



gumbel_formula_kke1 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  family = gumbelrank_family, center = FALSE
)

fit_kke1_gumbel <- brm(
    gumbel_formula_kke1, data = kk14_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kke1 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, center = FALSE
)

fit_kke1_uvsdt <- brm(
    uvsdt_formula_kke1, data = kk14_e1_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##------------
##  MHE22 E1  
##------------


gumbel_formula_mhe1 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  family = gumbelrank_family, center = FALSE
)

fit_mhe1_gumbel <- brm(
    gumbel_formula_mhe1, data = mhe_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


uvsdt_formula_mhe1 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, center = FALSE
)

fit_mhe1_uvsdt <- brm(
    uvsdt_formula_mhe1, data = mhe_e1_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

save(fit_kke1_gumbel, fit_kke1_uvsdt, fit_kks_gumbel, fit_kks_uvsdt, 
     fit_mhe1_gumbel, fit_mhe1_uvsdt, file = "fit-4rank.rda", compress = "xz")
load("fit-4rank.rda")

##---------------------------------------------------------------
##                            3 Ranks                           -
##---------------------------------------------------------------

##-----------
##  KK14 E2  
##-----------



gumbel_formula_kke2 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  family = gumbelrank_family, center = FALSE
)

fit_kke2_gumbel <- brm(
    gumbel_formula_kke2, data = kk14_e2_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kke2 <- brmsformula(
  observed | vint(rank, maxrank) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, center = FALSE
)

fit_kke2_uvsdt <- brm(
    uvsdt_formula_kke2, data = kk14_e2_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##-----------
##  MG16 E1  
##-----------
#McAdoo and Gronlund (2016)

fit_mge1_gumbel <- brm(
    gumbel_formula_kke2, data = mg16_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )
fit_mge1_uvsdt <- brm(
    uvsdt_formula_kke2, data = mg16_e1_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##-----------
##  MG16 E2  
##-----------



fit_mge2_gumbel <- brm(
    gumbel_formula_kke2, data = mg16_e2_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )
fit_mge2_uvsdt <- brm(
    uvsdt_formula_kke2, data = mg16_e2_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

save(fit_kke2_gumbel, fit_kke2_uvsdt, 
     fit_mge1_gumbel, fit_mge1_uvsdt,
     fit_mge2_gumbel, fit_mge2_uvsdt, file = "fit-3rank.rda", compress = "xz")
load("fit-3rank.rda")


##----------------------------------------------------------------
##                        Multiple Ranks                         -
##----------------------------------------------------------------
### Meyer-Grant and Jakob (2024)


gumbel_formula_mgj <- brmsformula(
  observed | vint(rank, maxrank) ~ 1 + (1|p|id), 
  family = gumbelrank_family, center = FALSE
)

fit_mgj_gumbel <- brm(
  gumbel_formula_mgj, data = mgj24, 
  stanvars = sv_gumbelrank, 
  prior = gumbel_priors,
  init_r = 0.25
  #control = list(adapt_delta = 0.99)
)

uvsdt_formula_mgj <- brmsformula(
  observed | vint(rank, maxrank) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, center = FALSE
)

fit_mgj_uvsdt <- brm(
  uvsdt_formula_mgj, data = mgj24, 
  stanvars = sv_uvsdtrank, 
  prior = uvsdt_priors,
  init_r = 0.25
  #control = list(adapt_delta = 0.99)
)

save(fit_mgj_gumbel, fit_mgj_uvsdt, 
     file = "fit-multirank.rda", compress = "xz")
load("fit-multirank.rda")

##----------------------------------------------------------------
##                              Plot                             -
##----------------------------------------------------------------

all_rank_data <- bind_rows(
  kks12l,
  kk14_e1_use, kk14_e2_use,
  mg16_e1_use, mg16_e2_use,
  mhe_e1_use,
  mgj24
) %>% 
  mutate(
    exp = factor(exp, levels = c("Kellen (2012)", "Kellen (2014, E1)", 
                                 "Malejka (2022, E1)",
                                 "Kellen (2014, E2)", 
                                 "McAdoo (2016, E1)", "McAdoo (2016, E2)",
                                 "Meyer-Grant (2024)"))
  )

rankfit_gumbel <- list(
  fit_kks_gumbel,
  fit_kke1_gumbel, fit_kke2_gumbel,
  fit_mge1_gumbel, fit_mge2_gumbel,
  fit_mhe1_gumbel,
  fit_mgj_gumbel
)
names(rankfit_gumbel) <- unique(all_rank_data$exp)

gumbel_samp_params <- map(rankfit_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(gumbel_samp_params)) {
  if (i == 1) cat("Number divergent transistions Gumbel:\n")
  cat(names(rankfit_gumbel)[i], ": ", 
      sum(map_dbl(gumbel_samp_params[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
}
# Number divergent transistions Gumbel:
# Kellen (2012) :  0 
# Kellen (2014, E1) :  0 
# Kellen (2014, E2) :  2 
# McAdoo (2016, E1) :  0 
# McAdoo (2016, E2) :  0 
# Malejka (2022, E1) :  1 
# Meyer-Grant (2024) :  0 

max(vapply(rankfit_gumbel, get_max_rhat, 0)) 
# 1.009762

rankfit_uvsdt <- list(
  fit_kks_uvsdt,
  fit_kke1_uvsdt, fit_kke2_uvsdt,
  fit_mge1_uvsdt, fit_mge2_uvsdt,
  fit_mhe1_uvsdt,
  fit_mgj_uvsdt
)
names(rankfit_uvsdt) <- unique(all_rank_data$exp)

max(vapply(rankfit_uvsdt, get_max_rhat, 0)) 
# 1.007562

uvsdt_samp_params <- map(rankfit_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(uvsdt_samp_params)) {
  if (i == 1) cat("Number divergent transistions UVSDT:\n")
  cat(names(rankfit_uvsdt)[i], ": ", 
      sum(map_dbl(uvsdt_samp_params[[i]], ~sum(.[1001:2000,"divergent__"]))),
      "\n")
}



pred_uvsdt <- lapply(rankfit_uvsdt, posterior_epred)
pred_gumbel <- lapply(rankfit_gumbel, posterior_epred)

all_rank_data$uvsd <- unlist(map(pred_uvsdt, ~apply(., c(2), mean)))
all_rank_data$gumbel <- unlist(map(pred_gumbel, ~apply(., c(2), mean)))

all_rank_data <- all_rank_data %>% 
  group_by(exp, id, strength, maxrank) %>% 
  mutate(prob = observed/sum(observed)) %>% 
  ungroup()

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

plot_rank_data <- all_rank_data %>% 
  group_by(exp, strength, rank, maxrank) %>% 
  summarise(across(c(prob, gumbel, uvsd), mean)) %>% 
  ungroup() %>% 
  mutate(rank = factor(rank), maxrank = factor(maxrank)) %>% 
  mutate(st_mxrnk = paste0(strength, "-", maxrank))

bin_n <- all_rank_data %>% 
  group_by(exp) %>% 
  summarise(n = n_distinct(id)) %>% 
  mutate(n_text = paste0("italic(N) == ", n))

psize <- 3.5
lsize <- 1.5

plrd1 <- plot_rank_data %>%
  filter(exp %in% c("Kellen (2014, E2)", 
                    "McAdoo (2016, E1)", "McAdoo (2016, E2)"))
plrd2 <- plot_rank_data %>%
  filter(exp %in% c("Kellen (2012)", "Kellen (2014, E1)", 
                    "Malejka (2022, E1)"))
plrd3 <- plot_rank_data %>%
  filter(exp %in% c("Meyer-Grant (2024)")) %>% 
  mutate(newexp = paste0(substr(exp, 1, 17), ", K", maxrank, ")"))
newn <- plrd3 %>% 
  select(exp, newexp) %>% 
  unique() %>% 
  left_join(bin_n) %>% 
  slice(1)

width_box <- 0.5
ylim <- c(0.07, 0.83)
ybreaks <- c(0.25, 0.5, 0.75)
########### v2 #########

p1b <- plrd1 %>%
  ggplot(aes(x = rank, y = prob)) +
  annotate(geom = "rect", xmin = 1 - width_box/2, xmax = 1 + width_box/2, 
           ymin = 0, ymax = 1/3,
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  # geom_hline(yintercept = 1/3, colour = rgb(0.7, 0.7, 0.7, alpha = 0.4), 
  #            linetype = 2) +
  geom_line(aes(group = strength, linetype = strength), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = Inf, y = Inf, label = n_text),
             data = filter(bin_n, exp %in% unique(plrd1$exp)),
             hjust = 1.1,
             vjust = 1.2,
             parse = TRUE, family = "Palatino Linotype") +
  #facet_wrap(vars(exp), ncol = 1, dir = "v", scales = "free_y") +
  facet_wrap(vars(exp), nrow = 1, dir = "v") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = expression(italic(R)[italic(i)]^{scriptstyle(paste("\u27E8", italic(K) == 3, "\u27E9"))})) +
  scale_linetype_manual(breaks = c("r", "s", "w"), 
                        values = c(1,1,2), guide = NULL) +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = ybreaks)

p2b <- plrd2 %>%
  ggplot(aes(x = rank, y = prob)) +
  annotate(geom = "rect", xmin = 1 - width_box/2, xmax = 1 + width_box/2, 
           ymin = 0, ymax = 1/4,
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  # geom_hline(yintercept = 1/4, colour = rgb(0.7, 0.7, 0.7, alpha = 0.4),
  #            linetype = 2) +
  geom_line(aes(group = strength, linetype = strength), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = Inf, y = Inf, label = n_text),
             data = filter(bin_n, exp %in% unique(plrd2$exp)),
             hjust = 1.1,
             vjust = 1.2,
             parse = TRUE, family = "Palatino Linotype") +
  #facet_wrap(vars(exp), ncol = 1, dir = "v", scales = "free_y") + 
  facet_wrap(vars(exp), nrow = 1, dir = "v") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = expression(italic(R)[italic(i)]^{scriptstyle(paste("\u27E8", italic(K) == 4, "\u27E9"))})) +
  #theme(axis.title.y = element_blank()) +
  scale_linetype_manual(breaks = c("r", "s", "w"), 
                        values = c(1,1,2), guide = NULL) +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = ybreaks)

(p1b / p2b ) +
  plot_layout(guides = 'collect', axes = "collect") 
#plot_annotation(tag_levels = list(c("A", "", "B")))
ggsave("rank-plot-4a.pdf", width = 19, height = 12, units = "cm")

plrd3 <- plrd3 %>% 
  mutate(newexp2 = factor(maxrank, 
                          levels = as.character(3:5), 
                          labels = c("italic(K) == 3", 
                                     "italic(K) == 4",
                                     "italic(K) == 5")
                          ))

p3b_1 <- plrd3 %>%
  filter(maxrank == 3) %>% 
  ggplot(aes(x = rank, y = prob)) +
  annotate(geom = "rect", xmin = 1 - width_box/2, xmax = 1 + width_box/2, 
           ymin = 0, ymax = 1/3,
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_line(aes(group = maxrank, linetype = maxrank), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = Inf, y = Inf, label = n_text),
             data = newn,
             hjust = 1.1,
             vjust = 1.2,
             parse = TRUE, family = "Palatino Linotype") +
  #facet_wrap(vars(newexp), ncol = 1, dir = "v", scales = "free_y") +
  facet_wrap(vars(newexp2), nrow = 1, dir = "v", scales = "free_x", 
             labeller = "label_parsed") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = expression(italic(R)[italic(i)]^{scriptstyle(paste("\u27E8", italic(K) == 3, "\u27E9"))})) +
  #theme(axis.title.y = element_blank()) +
  scale_linetype_manual(breaks = c("3", "4", "5"), 
                        values = c(1,1,1), guide = NULL) +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = ybreaks)

p3b_2 <- plrd3 %>%
  filter(maxrank == 4) %>% 
  ggplot(aes(x = rank, y = prob)) +
  annotate(geom = "rect", xmin = 1 - width_box/2, xmax = 1 + width_box/2, 
           ymin = 0, ymax = 1/4,
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_line(aes(group = maxrank, linetype = maxrank), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  #facet_wrap(vars(newexp), ncol = 1, dir = "v", scales = "free_y") +
  facet_wrap(vars(newexp2), nrow = 1, dir = "v", scales = "free_x", 
             labeller = "label_parsed") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = expression(italic(R)[italic(i)]^{scriptstyle(paste("\u27E8", italic(K) == 4, "\u27E9"))})) +
  #theme(axis.title.y = element_blank()) +
  scale_linetype_manual(breaks = c("3", "4", "5"), 
                        values = c(1,1,1), guide = NULL)+
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(breaks = ybreaks)
  #coord_cartesian(ylim = ylim, xlim = c(1, 5))

p3b_3 <- plrd3 %>%
  filter(maxrank == 5) %>% 
  ggplot(aes(x = rank, y = prob)) +
  annotate(geom = "rect", xmin = 1 - width_box/2, xmax = 1 + width_box/2, 
           ymin = 0, ymax = 1/5,
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_line(aes(group = maxrank, linetype = maxrank), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  #facet_wrap(vars(newexp), ncol = 1, dir = "v", scales = "free_y") +
  #facet_wrap(vars(newexp), nrow = 1, dir = "v", scales = "free_x") + 
  facet_wrap(vars(newexp2), nrow = 1, dir = "v", scales = "free_x", 
             labeller = "label_parsed") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = expression(italic(R)[italic(i)]^{scriptstyle(paste("\u27E8", italic(K) == 5, "\u27E9"))})) +
  #theme(axis.title.y = element_blank()) +
  scale_linetype_manual(breaks = c("3", "4", "5"), 
                        values = c(1,1,1), guide = NULL)+
  coord_cartesian(ylim = ylim, xlim = c(1, 5)) +
  scale_y_continuous(breaks = ybreaks)

p3b_1/p3b_2/p3b_3 + 
  plot_layout(guides = 'collect', axes = "collect") 
ggsave("rank-plot-4b.pdf", width = 8.5, height = 16.5, units = "cm")

p3b <- plrd3 %>%
  ggplot(aes(x = rank, y = prob)) +
  geom_hline(aes(yintercept = guess), 
             colour = rgb(0.7, 0.7, 0.7, alpha = 0.4), data = plrd3mr,
             linetype = 2) +
  geom_line(aes(group = maxrank, linetype = maxrank), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(y = gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(y = uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = Inf, y = Inf, label = n_text),
             data = newn,
             hjust = 1.1,
             vjust = 1.2,
             parse = TRUE, family = "Palatino Linotype") +
  #facet_wrap(vars(newexp), ncol = 1, dir = "v", scales = "free_y") +
  facet_wrap(vars(newexp), nrow = 1, dir = "v", scales = "free_x") + 
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
  theme(legend.title = NULL)  +
  labs(x = expression("Old-Item Rank" ~ group("(", italic(i), ")")), 
       y = "Pr(Old item rank)") +
  #theme(axis.title.y = element_blank()) +
  scale_linetype_manual(breaks = c("3", "4", "5"), 
                        values = c(1,1,1), guide = NULL)+
  coord_cartesian(ylim = ylim)

# (p1b / p2b / p3b) +
#   plot_layout(guides = 'collect', axes = "collect") 
# #plot_annotation(tag_levels = list(c("A", "", "B")))
# ggsave("rank-plot-3.pdf", width = 19, height = 18, units = "cm")



p3
ggsave("rank-plot-4b.pdf", width = 9, height = 18, units = "cm")
