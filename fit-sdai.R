library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

source("gumbelsdai-stan.R")
source("uvsdtsdai-stan.R")

rawrev <- read_csv("data_norev.csv")

d_sdai <- rawrev %>% 
  mutate(ident_resp = as.integer(ident_resp)) %>% 
  group_by(subject_id, ident_resp, conf_resp) %>%
  summarise(count = n(), .groups = "drop") %>% 
  mutate(
    conf_resp = paste0("c", conf_resp),
    ident_resp = paste0("i", ident_resp+2)
  ) %>% 
  pivot_wider(names_from = c(ident_resp, conf_resp), 
              values_from = count, values_fill = 0)
### i1: ident_resp = -1 
### i2: ident_resp = 0
### i3: ident_resp = 1

##### UVSD ####

uvsdt_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crl") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crh") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept)


uvsdt_formula <- brmsformula(
  i1_c1 | vint(i1_c2, i1_c3, i1_c4, i2_c1, i2_c2, i2_c3, i2_c4, i3_c1, i3_c2, i3_c3, i3_c4) ~ 1 + (1|p|subject_id), 
  discsignal ~ 1 + (1|p|subject_id), 
  crc ~ (1|p|subject_id), 
  crl ~ (1|p|subject_id), 
  crh ~ (1|p|subject_id), 
  family = uvsdtsdai_family, cmc = FALSE
)

fit_uvsd <- brm(
  uvsdt_formula, data = d_sdai, 
  stanvars = sv_uvsdtsdai ,
  prior = uvsdt_priors,
  init_r = 0.5
)

#prep_u <- prepare_predictions(fit_uvsd)
# xxx <- calc_posterior_predictions_uvsdtsdai(2, prep_u)
# str(xxx)

gumbel_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crl") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crh") +
  prior(student_t(3, 1, 2), class = Intercept)


gumbel_formula <- brmsformula(
  i1_c1 | vint(i1_c2, i1_c3, i1_c4, i2_c1, i2_c2, i2_c3, i2_c4, i3_c1, i3_c2, i3_c3, i3_c4) ~ 1 + (1|p|subject_id), 
  crc ~ (1|p|subject_id), 
  crl ~ (1|p|subject_id), 
  crh ~ (1|p|subject_id), 
  family = gumbelsdai_family
)

fit_gumbel <- brm(
  gumbel_formula, data = d_sdai, 
  stanvars = sv_gumbelsdai,
  prior = gumbel_priors,
  init_r = 0.1,
  control = list(adapt_delta = 0.9999999, max_treedepth = 20)
)

save(fit_uvsd, fit_gumbel, file = "sdai-fits.rda", compress = "xz")
load("sdai-fits.rda")

##### plot

res_agg <- d_sdai %>% 
  summarise(across(c(i1_c1, i1_c2, i1_c3, i1_c4, i2_c1, i2_c2, i2_c3, i2_c4, i3_c1, i3_c2, i3_c3, i3_c4), sum)) %>% 
  as.data.frame()
res_agg[,1:4] <- res_agg[,1:4]/sum(res_agg[,1:4])
res_agg[,5:12] <- res_agg[,5:12]/sum(res_agg[,5:12])
dplot2 <- res_agg %>% 
  pivot_longer(everything(), names_to = c("type", "crit"), names_sep = "_") %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  arrange(desc(crit)) %>% 
  mutate(fa = cumsum(i1),
         hi = cumsum(i3),
         h = cumsum(i2 + i3)) %>% 
  select(-i1, -i2, -i3) %>% 
  pivot_longer(c(hi, h))

pred_uvsdt <- posterior_epred(fit_uvsd)
res_uvsdt <- as.data.frame(t(colMeans(apply(pred_uvsdt, c(2,3), mean))))
colnames(res_uvsdt) <- colnames(res_agg)
dplot2_uvsdt <- res_uvsdt %>% 
  pivot_longer(everything(), names_to = c("type", "crit"), names_sep = "_") %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  arrange(desc(crit)) %>% 
  mutate(fa = cumsum(i1),
         hi = cumsum(i3),
         h = cumsum(i2 + i3)) %>% 
  select(-i1, -i2, -i3) %>% 
  pivot_longer(c(hi, h))

pred_gumbel <- posterior_epred(fit_gumbel)
res_gumbel <- as.data.frame(t(colMeans(apply(pred_gumbel, c(2,3), mean))))
colnames(res_gumbel) <- colnames(res_agg)
dplot2_gumbel <- res_gumbel %>% 
  pivot_longer(everything(), names_to = c("type", "crit"), names_sep = "_") %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  arrange(desc(crit)) %>% 
  mutate(fa = cumsum(i1),
         hi = cumsum(i3),
         h = cumsum(i2 + i3)) %>% 
  select(-i1, -i2, -i3) %>% 
  pivot_longer(c(hi, h))


# dplot <- res_agg %>% 
#   pivot_longer(everything(), names_to = c("type", "crit"), names_sep = "_") %>% 
#   group_by(type) %>% 
#   mutate(roc = cumsum(rev(value))) %>% 
#   ungroup()
# dplot2 <- dplot %>% 
#   filter(type == "i1") %>% 
#   select(-type, - value) %>% 
#   rename(fa = roc) %>% 
#   right_join(filter(dplot, type != "i1"))
# dplot2 %>% 
#   arrange(type)

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
stsize <- 0.5
dplot2 %>%
  ggplot(aes(x =  fa, y = value)) +
  #geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "white") +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), 
           fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, -Inf, 1/2 + 0.025), 
           fill = rgb(0.3, 0.3, 0.3, alpha = 0.4)) +
  geom_abline(slope = 1/2, intercept = 0, linetype = 2) +
  geom_line(aes(group = name), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x =  fa, y = value,
                 shape = "UVSD", colour = "UVSD"),
             size = psize, stroke = stsize, data = dplot2_uvsdt) +
  geom_point(aes(x =  fa, y = value,
                 shape = "Gumbel", colour = "Gumbel"),
             size = psize, stroke = stsize, data = dplot2_gumbel) +
  annotate("label", x = 0.75, y = 0.15, label = "italic(N) == 48",
             hjust = "center", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1.05), ylim = c(0, 1.05), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
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
  labs(x = expression(italic(p)["D-FA"]), 
       y = expression(italic(p)["D-HI"]~"&"~italic(p)["D-H"]))
ggsave("fit-sdai.pdf", width = 8.1, height = 7.5, units = "cm")
