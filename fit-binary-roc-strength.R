
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 12) + 
            theme(legend.position="bottom"))

source("dube-bin-strength-data.R")
## Dubet et al. (2012, JML)
source("gumbelbinsep-stan.R")
source("uvsdtbinsep-stan.R")

#### new data sets

head(dbin_dube)

dbin_dube <- dbin_dube %>% 
  mutate(
    total = O+N,
    newstrength = if_else(strength %in% c("S","N"), "SN", "W"),
    oldnew = if_else(strength %in% c("S","W"), 0, 1)
  )

datasets_all <- unique(dbin_dube$exp)

gumbel_formula <- brmsformula(
  O | vint(total, oldnew) ~ 0 + newstrength + (0 + newstrength|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = gumbelbinsep_family, cmc = TRUE
)

get_prior(gumbel_formula, data = dbin_dube)

gumbel_priors <- prior(student_t(3, 1, 2), class = b)
### prior(normal(0,0.5), class = b, dpar = "cr") + 

# make_stancode(gumbel_formula, family = gumbelbinsep_family, 
#               priors = gumbel_priors, data = dbin_dube)
# 
# xxx <- make_standata(gumbel_formula, family = gumbelbinsep_family, 
#                      priors = gumbel_priors, data = dbin_dube)
# str(xxx, 1)
# head(xxx$X, 20)
# head(xxx$X_cr, 20)

uvsdt_formula <- brmsformula(
  O | vint(total, oldnew) ~ 0 + newstrength + (0 + newstrength|p|pid),  
  discsignal ~ 1 + (1|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = uvsdtbinsep_family, cmc = TRUE
)

uvsdt_priors <- prior(student_t(3, 1, 2), class = "b") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal")
# prior(normal(0,0.5), class = b, dpar = "cr") + 


dube_bin_data <- vector("list", length(datasets_all))

dube_bin_fits_gumbel <- vector("list", length(datasets_all))
dube_bin_fits_uvsdt <- vector("list", length(datasets_all))


#set.seed(4567123)
for (i in seq_along(datasets_all)) {
#for (i in c(1,2)) {
  print(i)
  dube_bin_data[[i]] <- dbin_dube %>% 
    filter(exp == datasets_all[i])
  dube_bin_fits_gumbel[[i]] <- brm(
    gumbel_formula, data = dube_bin_data[[i]], 
    stanvars = sv_gumbelbinsep ,
    prior = gumbel_priors,
    init_r = 0.5, 
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  )
  dube_bin_fits_uvsdt[[i]] <- brm(
    uvsdt_formula, data = dube_bin_data[[i]], 
    stanvars = sv_uvsdtbinsep, 
    prior = uvsdt_priors,
    init_r = 0.5, 
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  )
}



### make plots

pred_gumbel <- lapply(dube_bin_fits_gumbel, posterior_epred)
names(pred_gumbel) <- datasets_all
pred_uvsd <- lapply(dube_bin_fits_uvsdt, posterior_epred)
names(pred_uvsd) <- datasets_all

str(pred_gumbel, 1)

plot_dat <- dbin_dube %>% 
  filter(exp %in% datasets_all) %>% 
  mutate(
    prob = O/total
  ) 

plot_dat %>% 
  group_by(exp) %>% 
  summarise(n = n()) 

plot_dat$gumbel <- unname(unlist(map(pred_gumbel, ~apply(., c(2), mean))))
plot_dat$uvsd <- unname(unlist(map(pred_uvsd, ~apply(., c(2), mean))))

plot_dat_2 <- plot_dat %>% 
  group_by(exp, strength, newstrength, baserate, oldnew) %>%
  summarise(across(c(prob, gumbel, uvsd), mean)) %>% 
  ungroup()
plot_dat_2_new <- plot_dat_2 %>% 
  filter( strength == "N") %>% 
  select(-strength, -newstrength, -oldnew) %>% 
  rename(
    fa_prob = prob, fa_gumbel = gumbel, fa_uvsd = uvsd
  )
plot_dat_2_notnew <- plot_dat_2 %>% 
  filter( strength != "N")
plot_dat_2 <- left_join(plot_dat_2_notnew, plot_dat_2_new)

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

#hcl.colors(5, "Plasma")

psize <- 3.5
lsize <- 1.5
stsize <- 0.5
plot_dat_2 %>%
  ggplot(aes(x =  fa_prob, y = prob)) +
  geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "white") +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(aes(group = strength), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x = fa_uvsd, y = uvsd,
                 shape = "UVSD", colour = "UVSD"), 
             size = psize, stroke = stsize) +
  geom_point(aes(x = fa_gumbel, y = gumbel,
                 shape = "Gumbel", colour = "Gumbel"), 
             size = psize, stroke = stsize) +
  geom_label(mapping = aes(x = 0.75, y = 0.15, label = n_text), 
             data = bin_n, hjust = "center", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(exp), nrow = 1) + 
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
ggsave("binstrengthroc-plot1.pdf", 
       width = 14, height = 9.75, units = "cm")
