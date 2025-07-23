
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel6agg-strength-stan.R")
source("uvsdt6agg-strength-stan.R")
afex::set_sum_contrasts()

sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "id")

dratcliff <- read_delim("data_ratcliff.txt", 
                        col_names = paste0(rep(c("weak_r", "strong_r", "new_r"), each = 6), 1:6))
dratcliff <- dratcliff %>% 
  mutate(id = as.character(rep(1:11, 2))) %>% 
  mutate(frequency = rep(c("high", "low"), each = 11)) %>% 
  select(id, frequency, everything()) %>% 
  pivot_longer(cols = -c(id, frequency), 
               names_to = c("strength", "response"), names_sep = "_") %>% 
  pivot_wider(names_from = response, values_from = value) %>% 
  mutate(type = if_else(strength == "new", 0, 1)) %>% 
  mutate(newstrength = if_else(strength == "weak", "weak", "strong-new"))

gumbel_formula <- brmsformula(
  r1 | vint(r2, r3, r4, r5, r6, type) ~ newstrength + (newstrength|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id),
  family = gumbel6aggreg_family, cmc = FALSE
)

get_prior(gumbel_formula, data = dratcliff)

gumbel_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

uvsdt_formula <- brmsformula(
  r1 | vint(r2, r3, r4, r5, r6, type) ~ newstrength + (newstrength|p|id), 
  discsignal ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id),
  family = uvsdt6aggreg_family, cmc = FALSE
)

uvsdt_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

dset_strength <- c("high", "low")
rocstrength_data <- vector("list", length(dset_strength))

rocstrength_fits_gumbel <- vector("list", length(dset_strength))
rocstrength_fits_uvsdt <- vector("list", length(dset_strength))

control1 <- list(adapt_delta = 0.99, max_treedepth = 20)
iter <- 2000
warmup <- 1000

for (i in seq_along(dset_strength)) {
  print(i)
  rocstrength_data[[i]] <- dratcliff %>% 
    filter(frequency == dset_strength[i])
  rocstrength_fits_gumbel[[i]] <- brm(
    gumbel_formula, data = rocstrength_data[[i]], 
    stanvars = sv_gumbel6aggreg,
    prior = gumbel_priors,
    init_r = 0.5, 
    iter = iter, warmup = warmup,
    control = control1
  )

  rocstrength_fits_uvsdt[[i]] <- brm(
    uvsdt_formula, data = rocstrength_data[[i]], 
    stanvars = sv_uvsdt6aggreg, 
    prior = uvsdt_priors,
    init_r = 0.5, 
    iter = iter, warmup = warmup,
    control = control1
  )

}

### convergence stats
source("check-functions.R")
max(vapply(rocstrength_fits_gumbel, get_max_rhat, 0))
# 1.005573
max(vapply(rocstrength_fits_uvsdt, get_max_rhat, 0))
# 1.006043

xxx <- map(rocstrength_fits_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dset_strength)) {
  cat(dset_strength[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
}

xxy <- map(rocstrength_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dset_strength)) {
  cat(dset_strength[i], ": ", sum(map_dbl(xxy[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
}


########### PLOT ########

pred_gumbel <- lapply(rocstrength_fits_gumbel, posterior_epred)
names(pred_gumbel) <- dset_strength
pred_uvsd <- lapply(rocstrength_fits_uvsdt, posterior_epred)
names(pred_uvsd) <- dset_strength

str(pred_gumbel, 1)

plot_dat <- dratcliff %>% 
  mutate(sum = r1 + r2 + r3 + r4 + r5 + r6) %>% 
  mutate(across(r1:r6, ~./sum)) %>% 
  select(-sum)

pred_gumbel2 <- do.call("rbind", map(pred_gumbel, ~apply(., c(2, 3), mean))) %>% 
  as.data.frame()
colnames(pred_gumbel2) <- str_replace(colnames(pred_gumbel2), "r", "g")
pred_uvsd2 <- do.call("rbind", map(pred_uvsd, ~apply(., c(2, 3), mean))) %>% 
  as.data.frame()
colnames(pred_uvsd2) <- str_replace(colnames(pred_uvsd2), "r", "u")

plot_dat <- plot_dat %>% 
  bind_cols(pred_gumbel2) %>% 
  bind_cols(pred_uvsd2)

plot_dat_2 <- plot_dat %>% 
  pivot_longer(cols = c(r1:r6, g1:g6, u1:u6)) %>% 
  mutate(type = substr(name, 1, 1)) %>% 
  mutate(response = substr(name, 2, 2)) %>% 
  group_by(frequency, strength, type, response) %>% 
  summarise(prob = mean(value)) %>% 
  arrange(frequency, strength, type, rev(response)) %>% 
  group_by(frequency, strength, type) %>% 
  mutate(prob = cumsum(prob)) %>% 
  filter( response != "1") %>% 
  pivot_wider(names_from = type, values_from = prob) %>% 
  ungroup()


plot_dat_2_new <- plot_dat_2 %>% 
  filter( strength == "new") %>% 
  select(-strength) %>% 
  rename(
    fa_prob = r, fa_gumbel = g, fa_uvsd = u
  )
plot_dat_2_notnew <- plot_dat_2 %>% 
  filter( strength != "new")
plot_dat_2 <- left_join(plot_dat_2_notnew, plot_dat_2_new)

bin_n <- plot_dat %>% 
  group_by(frequency) %>% 
  summarise(n = n_distinct(id)) %>% 
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

#hcl.colors(5, "Plasma")

psize <- 3.5
lsize <- 1.5
stsize <- 1.0
plot_dat_2 %>%
  ggplot(aes(x =  fa_prob, y = r)) +
  geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "white") +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(aes(group = strength), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x = fa_uvsd, y = u,
                 shape = "UVSD", colour = "UVSD"), 
             size = psize, stroke = stsize) +
  geom_point(aes(x = fa_gumbel, y = g,
                 shape = "Gumbel", colour = "Gumbel"), 
             size = psize, stroke = stsize) +
  geom_label(mapping = aes(x = 0.75, y = 0.15, label = n_text), 
             data = bin_n, hjust = "center", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(frequency), ncol = 1, 
             labeller = as_labeller(c("high" = "High Frequency", 
                                      "low" = "Low Frequency"))) + 
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
ggsave("rocstrength-plot1.pdf", 
       width = 8.1, height = 13.5, units = "cm")

