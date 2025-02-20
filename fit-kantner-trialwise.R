library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

load("dat-kantner-full.rda")

###### replicate paper (Table 1):

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

#####

source("gumbel24afctrial-stan.R")

### goal: include trial-level predictors and crossed random effects and show it gives clearer results as UVSD where trial-level predictors have to go on both parameters
### 

dat_kantner_prep <- dat_kantner_prep %>% 
  mutate(afc = case_when(
    trial_type == "TwoAFCProc" ~ 2L,
    trial_type == "FourAFCProc" ~ 4L
  ))

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

pred_gumbel <- posterior_epred(fit_acoustic_gumbel)
str(pred_gumbel)

pred_gumbel_mean <- apply(pred_gumbel, MARGIN = c(2,3), mean)

dacoustic <- dacoustic %>% 
  mutate(acc_2 = corr_2/n_2,
         acc_4 = corr_4/n_4)

dacoustic_pred <- dacoustic %>% 
  bind_cols(as_tibble(pred_gumbel_mean))

dacoustic_pred %>% 
  select(-c(corr_2:n_4)) %>% 
  pivot_longer(cols = c(acc_2:AFC4)) %>% 
  mutate(afc = str_extract(name, "\\d"),
         type = substr(name, 1, 3)) %>% 
  select(-name) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  ggplot(aes(x = acc, y = AFC)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(vars(afc)) + 
  coord_fixed(xlim = c(0.2, 1), ylim = c(0.2, 1)) +
  labs(x = "data", y = "prediction")

