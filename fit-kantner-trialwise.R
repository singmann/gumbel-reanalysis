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

source("gumbel24afctrial-stan.R")

### goal: include trial-level predictors and crossed random effects and show it gives clearer results as UVSD where trial-level predictors have to go on both parameters
### 


gumbel_priors <- prior(student_t(3, 1, 2), class = Intercept)

gumbel_formula <- brmsformula(
  acc | vint(afc) ~ z_onset+z_response_dur+z_mean_pitch+z_mean_intensity +
    (z_onset+z_response_dur+z_mean_pitch+z_mean_intensity|p|subject),
  family = gumbel24afctrial_family
)

# gumbel_formula <- brmsformula(
#   acc | vint(afc) ~ 1 + (1|p|subject),
#   family = gumbel24afctrial_family
# )

# make_stancode(gumbel_formula, data = dat_kantner_prep, 
#               stanvars = sv_gumbel24afctrial, 
#               prior = gumbel_priors)
fit_acoustic_gumbel <- brm(
  gumbel_formula, data = dat_kantner_prep, 
  stanvars = sv_gumbel24afctrial, 
  prior = gumbel_priors,
  init_r = 0.25
  #control = list(adapt_delta = 0.99)
)
summary(fit_acoustic_gumbel)

cbind(kantner = exp(lme4::fixef(m1a))[-6],
      gumbel = exp(lme4::fixef(fit_acoustic_gumbel)[,1])) %>% 
  round(2)

pred_gumbel <- posterior_epred(fit_acoustic_gumbel, cores = 10)
save(fit_acoustic_gumbel, pred_gumbel_mean, 
     file = "fit-acoustic-indiv.rda", compress = "xz")

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

load("fit-acoustic-indiv.rda")
# xxx <- prepare_predictions(fit_acoustic_gumbel)
# calc_posterior_predictions_gumbel24afctrial(1, xxx)
str(pred_gumbel)
pred_gumbel_mean <- apply(pred_gumbel, MARGIN = 2, mean)

dat_kantner_prep2 <- dat_kantner_prep %>% 
  filter(!is.na(z_onset), !is.na(z_response_dur), !is.na(z_mean_pitch), !is.na(z_mean_intensity)) %>% 
  mutate(gumbel = pred_gumbel_mean)

dat_kantner_prep2 %>% 
  group_by(subject, afc) %>% 
  summarise(observed = mean(acc),
            gumbel = mean(gumbel)) %>% 
  ggplot(aes(x = observed, y = gumbel)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(vars(afc)) + 
  coord_fixed(xlim = c(0.2, 1), ylim = c(0.2, 1)) +
  labs(x = "data", y = "prediction")

dat_kantner_prep2 %>% 
  group_by(subject, afc) %>% 
  summarise(observed = mean(acc),
            gumbel = mean(gumbel)) %>% 
  ungroup() %>% 
  summarise(cor(observed, gumbel))
