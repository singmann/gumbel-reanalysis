library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

load("dat-kantner-full.rda")
dat_kantner_prep <- dat_kantner_prep %>% 
  mutate(afc = case_when(
    trial_type == "TwoAFCProc" ~ 2L,
    trial_type == "FourAFCProc" ~ 4L
  ))

dat_kantner_prep3 <- dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  filter(z_response_dur < 15)

###### replicate paper (Table 1):

dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  select(z_onset,z_response_dur,z_mean_pitch,z_mean_intensity) %>% 
  summarise(across(everything(), c(min = min, max = max))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(measure = if_else(str_detect(name, "min"), "min", "max"),
         name = str_remove(name, "_min|_max")) %>% 
  pivot_wider(names_from = measure) %>% 
  mutate(range = max - min)

dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  select(z_onset,z_response_dur,z_mean_pitch,z_mean_intensity) %>% 
  summarise(across(everything(), c(mean)))

dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  select(z_onset,z_response_dur,z_mean_pitch,z_mean_intensity) %>% 
  GGally::ggpairs()


m1a <-  lme4::glmer(acc~z_onset+z_response_dur+z_mean_pitch+z_mean_intensity+trial_type+
             #(1|subject)+
             (0+z_onset|subject)+
             (0+z_response_dur|subject)+
             (0+z_mean_pitch|subject)+
             (0+z_mean_intensity|subject)
           ,family=binomial,
           data=dat_kantner_prep,
           lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m1a))

round(exp(lme4::fixef(m1a)), 2)


m2a <-  lme4::glmer(acc ~ z_onset + z_response_dur + poly(z_mean_pitch, 2) + poly(z_mean_intensity, 2) + 
                      trial_type+
                      #(1|subject)+
                      (0+z_onset|subject)+
                      (0+z_response_dur|subject)+
                      (0+z_mean_pitch|subject)+
                      (0+z_mean_intensity|subject)
                    ,family=binomial,
                    data=dat_kantner_prep3,
                    lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
print(summary(m2a))

round(exp(lme4::fixef(m1a)), 2)
#####

## all files for fitting 2-AFC and 4-AFC Gumbel model 
source("gumbelafctrial-stan.R")
sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "subject")

### goal: include trial-level predictors and crossed random effects for Gumbel_min
### result: provides clearer interpretation compared to unequal-variance Gaussian model,
### as trial-level predictors have to go on both parameters in Gaussian model 
### (i.e., mean and variance)


gumbel_priors <- prior(student_t(3, 1, 2), class = b, coef = Intercept) +
  sd_priors

gumbel_formula <- brmsformula(
  acc | vint(afc) ~ z_onset+z_response_dur+z_mean_pitch+z_mean_intensity +
    (z_onset+z_response_dur+z_mean_pitch+z_mean_intensity|p|subject),
  family = gumbelafctrial_family, center = FALSE
)

# gumbel_formula <- brmsformula(
#   acc | vint(afc) ~ 1 + (1|p|subject),
#   family = gumbel24afctrial_family
# )

# get_prior(gumbel_formula, data = dat_kantner_prep,
#               stanvars = sv_gumbelafctrial,
#               #prior = gumbel_priors
# )

# make_stancode(gumbel_formula, data = dat_kantner_prep,
#               stanvars = sv_gumbelafctrial,
#               prior = gumbel_priors
#               )
fit_acoustic_gumbel <- brm(
  gumbel_formula, data = dat_kantner_prep, 
  stanvars = sv_gumbelafctrial, 
  prior = gumbel_priors,
  init_r = 0.25,
  #control = list(adapt_delta = 0.99, max_treedepth = 20)
)
summary(fit_acoustic_gumbel)

### comparison of estimates across models
cbind(kantner = exp(lme4::fixef(m1a))[-c(1,6)],
      gumbel = exp(lme4::fixef(fit_acoustic_gumbel)[-1,1])) %>% 
  round(2)
#                  kantner gumbel
# z_onset             0.68   0.69
# z_response_dur      0.78   0.82
# z_mean_pitch        1.05   1.05
# z_mean_intensity    1.02   1.06

pred_gumbel <- posterior_epred(fit_acoustic_gumbel)
str(pred_gumbel)
pred_gumbel_mean <- apply(pred_gumbel, MARGIN = 2, mean)

### create data that is actually fitted without missings:
dat_kantner_prep2 <- dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  mutate(gumbel = pred_gumbel_mean)

nrow(dat_kantner_prep2)
# [1] 7333
# same as: 
# summary(fit_acoustic_gumbel)
# [...]
# Data: dat_kantner_prep (Number of observations: 7333) 

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

### ,
### compare observed vs. predicted
dat_plot <- dat_kantner_prep2 %>% 
  group_by(subject, afc) %>% 
  summarise(observed = mean(acc),
            gumbel = mean(gumbel)) %>% 
  mutate(mafc = paste0(afc, "-AFC"))


dat_plot2 <- dat_plot %>% 
  group_by(mafc) %>% 
  summarise(afc = afc[1])

dat_plot %>% 
  ggplot(aes(x = gumbel, y = observed)) +
  geom_rect(aes(xmin = -Inf, xmax = 1/afc, ymin = -Inf, ymax = Inf), 
              data = dat_plot2, inherit.aes = FALSE,
            fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_rect(aes(xmin = 1/afc, xmax = Inf, ymin = -Inf, ymax = 1/afc), 
            data = dat_plot2, inherit.aes = FALSE,
            fill = rgb(0.7, 0.7, 0.7, alpha = 0.4)) +
  geom_segment(aes(x = 1/afc, xend = Inf, y = 1/afc),
               linetype = 2,
               data = dat_plot2, inherit.aes = FALSE) +
  geom_segment(aes(y = 1/afc, yend = Inf, x = 1/afc),
               linetype = 2,
               data = dat_plot2, inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(shape = 21, colour = "black", fill = "grey") +
  facet_wrap(vars(mafc), nrow = 2) + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 1), labels = c("0", ".25", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 1), labels = c("0", ".25", ".5", "1")) +
  labs(x = expression(
    paste("Predicted ", italic(R)[1]^{scriptstyle(paste("\u27E8", italic(K), "\u27E9"))})), 
       y = expression(
         paste("Observed ", italic(R)[1]^{scriptstyle(paste("\u27E8", italic(K), "\u27E9"))})))
ggsave("kantner-fit.pdf", 
       width = 8.1, height = 13.5, units = "cm")

dat_kantner_prep2 %>% 
  group_by(subject, afc) %>% 
  summarise(observed = mean(acc),
            gumbel = mean(gumbel)) %>% 
  ungroup() %>% 
  summarise(cor(observed, gumbel))


#### more complicated model with quadratic effects of all covariates
#### shows no quadratic effect of any covariate
gumbel_formula3 <- brmsformula(
  acc | vint(afc) ~ 
    z_onset + I(z_onset^2) +
    z_response_dur + I(z_response_dur^2) +
    z_mean_pitch + I(z_mean_pitch^2) +
    z_mean_intensity + I(z_mean_intensity^2) +
    (z_onset+z_response_dur+z_mean_pitch+z_mean_intensity|p|subject),
  family = gumbel24afctrial_family
)
fit_acoustic_gumbel3 <- brm(
  gumbel_formula3, data = dat_kantner_prep3, 
  stanvars = sv_gumbel24afctrial, 
  prior = gumbel_priors,
  init_r = 0.25
  #control = list(adapt_delta = 0.99)
)
summary(fit_acoustic_gumbel3)
pred_gumbel3 <- posterior_epred(fit_acoustic_gumbel3)
save(fit_acoustic_gumbel3, 
     file = "fit-acoustic-indiv3.rda", compress = "xz")