
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 12) + 
            theme(legend.position="bottom"))

load("malejka-broeder.rda")
source("bin-roc-data.R")
source("gumbelbin-stan.R")
source("uvsdtbin-stan.R")

sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "pid")

#### new data sets

head(dbin5point) ## miss, hit, cr, false alarm 

colnames(dbin5point) <- outer(c("miss", "hit", "cr", "fa"), c("cr1", "cr2", "cr3", "cr4", "cr5"), 
      FUN = "paste", sep = "_") %>% 
  as.vector()

bdin_all <- dbin5point %>% 
  as_tibble() 
bdin_all$id <- rownames(dbin5point)
bdin_all <- bdin_all %>% 
  separate(id, into = c("exp", "pid"), sep = ", ") %>% 
  mutate(pid = str_trim(pid))

bdin_long <- bdin_all %>%
  pivot_longer(cols = -c(exp, pid), 
               names_to = c("type", "baserate"), names_sep = "_") %>% 
  pivot_wider(names_from = "type") %>% 
  mutate(Nold = hit + miss, 
         Nnew = fa + cr)

mb_extra <- bind_rows(mbe1, mbe3) %>% 
  select(experiment, Subject, BaseRate, everything()) %>% 
  mutate(Subject = as.character(Subject))
colnames(mb_extra) <- colnames(bdin_long)
bdin_long <- bind_rows(bdin_long, mb_extra)
bdin_long <- bdin_long %>% 
  mutate(exp = factor(
    exp, 
    levels = c("Broder E3", 
               "Dube E1 pictures", "Dube E1 word", "Dube E2", 
                "e1", "Malejka_e2", "e3",
               "Van Zandt E1 slow", "Van Zandt E1 fast", "Van Zandt E2"), 
    labels = c("Broder (2009, E3)", 
               "Dube (2012, E1a-P)", "Dube (2012, E1a-W)", "Dube (2012, E2)", 
                "Malejka (2019, E1)", "Malejka (2019, E2)", "Malejka (2019, E3)",
               "Van Zandt (2000, E1-F)", "Van Zandt (2000, E1-S)", 
               "Van Zandt (2000, E2)")))

dataset_all <- levels(bdin_long$exp)
dataset_all <- dataset_all[-which(dataset_all %in% c("Malejka (2019, E1)", 
                                                     "Malejka (2019, E3)"))] 

gumbel_formula_2 <- brmsformula(
  hit | vint(Nold, fa, Nnew) ~ 1 + (1|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = gumbelbin_family
)

get_prior(gumbel_formula_2, 
          data = filter(bdin_long, exp == dataset_all[[1]]))

gumbel_priors <- prior(student_t(3, 1, 2), class = "Intercept") +
  sd_priors
## prior(normal(0,0.5), class = b, dpar = "cr") + 

uvsdt_formula_2 <- brmsformula(
  hit | vint(Nold, fa, Nnew) ~ 1 + (1|p|pid), 
  discsignal ~ 1 + (1|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = uvsdtbin_family
)

uvsdt_priors <- prior(student_t(3, 1, 2), class = "Intercept") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  sd_priors
# prior(normal(0,0.5), class = b, dpar = "cr") + 


rocbin_data <- vector("list", length(dataset_all))

rocbin_fits_gumbel <- vector("list", length(dataset_all))
rocbin_fits_uvsdt <- vector("list", length(dataset_all))

iter <- 3000
warmup <- 1000

#set.seed(4567123)
for (i in seq_along(dataset_all)) {
#for (i in c(1,2)) {
  print(i)
  rocbin_data[[i]] <- bdin_long %>% 
    filter(exp == dataset_all[i])
  rocbin_fits_gumbel[[i]] <- brm(
    gumbel_formula_2, data = rocbin_data[[i]], 
    stanvars = sv_gumbelbin ,
    prior = gumbel_priors,
    init_r = 0.5, 
    iter = iter, warmup = warmup,
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  )
  rocbin_fits_uvsdt[[i]] <- brm(
    uvsdt_formula_2, data = rocbin_data[[i]], 
    stanvars = sv_uvsdtbin, 
    prior = uvsdt_priors,
    init_r = 0.5, 
    iter = iter, warmup = warmup,
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  )
}

### convergence stats
source("check-functions.R")
max(vapply(rocbin_fits_gumbel, get_max_rhat, 0))
# 1.00445
max(vapply(rocbin_fits_uvsdt, get_max_rhat, 0))
# 1.009024

xxx <- map(rocbin_fits_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset_all)) {
  cat(dataset_all[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[(warmup+1):iter,"divergent__"])))/((iter-warmup)*4), "\n")
}
# Broder (2009, E3) :  0 
# Dube (2012, E1a-P) :  0 
# Dube (2012, E1a-W) :  0 
# Dube (2012, E2) :  0 
# Malejka (2019, E2) :  0 
# Van Zandt (2000, E1-F) :  0 
# Van Zandt (2000, E1-S) :  0 
# Van Zandt (2000, E2) :  0.00025 

xxx <- map(rocbin_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset_all)) {
  cat(dataset_all[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[(warmup+1):iter,"divergent__"])))/((iter-warmup)*4), "\n")
}
# Broder (2009, E3) :  0 
# Dube (2012, E1a-P) :  0 
# Dube (2012, E1a-W) :  0 
# Dube (2012, E2) :  0 
# Malejka (2019, E2) :  0 
# Van Zandt (2000, E1-F) :  0.00075 
# Van Zandt (2000, E1-S) :  0.000125 
# Van Zandt (2000, E2) :  0.000125 


### make plots

pred_gumbel <- lapply(rocbin_fits_gumbel, posterior_epred)
names(pred_gumbel) <- dataset_all
pred_uvsd <- lapply(rocbin_fits_uvsdt, posterior_epred)
names(pred_uvsd) <- dataset_all

str(pred_gumbel, 1)

plot_dat <- bdin_long %>% 
  filter(exp %in% dataset_all) %>% 
  group_by(exp, pid, baserate) %>% 
  summarise(
    hit = hit/Nold,
    fa = fa/Nnew,
  ) 

plot_dat %>% 
  ungroup() %>% 
  count(exp) 

plot_dat <- bind_cols(
  plot_dat, 
  map_dfr(pred_gumbel, ~as.data.frame(apply(., c(2, 3), mean)))
) %>% 
  rename(
    hit_gumbel = old, fa_gumbel = new
  )
plot_dat <- bind_cols(
  plot_dat, 
  map_dfr(pred_uvsd, ~as.data.frame(apply(., c(2, 3), mean)))
) %>% 
  rename(
    hit_uvsd = old, fa_uvsd = new
  )

plot_dat_2 <- plot_dat %>% 
  group_by(exp, baserate) %>% 
  summarise(across(c(hit,fa, hit_gumbel, fa_gumbel, hit_uvsd, fa_uvsd), mean))

plot_data_bin <- plot_dat_2


bin_n <- plot_dat %>% 
  group_by(exp) %>% 
  summarise(n = n_distinct(pid)) %>% 
  mutate(n_text = paste0("italic(N) == ", n))

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

psize <- 3.5
lsize <- 1.5
plot_dat_2 %>%
  ggplot(aes(x =  fa, y = hit)) +
  geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "white") +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x = fa_gumbel, y = hit_gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(x = fa_uvsd, y = hit_uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = 0.75, y = 0.15, label = n_text), 
             data = bin_n, hjust = "center", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(exp), nrow = 2) + 
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
ggsave("binroc-plot1.pdf", width = 22, height = 14.5, units = "cm")
