
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

load("malejka-broeder.rda")
source("gumbelbin-stan.R")
source("uvsdtbin-stan.R")

gumbel_formula <- brmsformula(
  OLD_old | vint(Nold, NEW_old, Nnew) ~ 1 + (1|p|Subject), 
  cr ~ 0 + BaseRate + (0 + BaseRate|p|Subject),
  family = gumbelbin_family, cmc = FALSE
)

get_prior(gumbel_formula, data = mbe1)

gumbel_priors <- prior(normal(0,0.5), class = b, dpar = "cr") + 
  prior(student_t(3, 1, 2), class = Intercept)

stancode(gumbel_formula, data = mbe1,
         stanvars = sv_gumbelbin,
         prior = gumbel_priors)

uvsdt_formula <- brmsformula(
  OLD_old | vint(Nold, NEW_old, Nnew) ~ 1 + (1|p|Subject), 
  discsignal ~ 1 + (1|p|Subject), 
  cr ~ 0 + BaseRate + (0 + BaseRate|p|Subject),
  family = uvsdtbin_family, cmc = FALSE
)

uvsdt_priors <- prior(normal(0,0.5), class = b, dpar = "cr") + 
  prior(student_t(3, 1, 2), class = Intercept) +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal")

stancode(uvsdt_formula, data = mbe1,
         stanvars = sv_uvsdtbin,
         prior = uvsdt_priors)

### fitting

mbfit_e1_gumbel <- brm(
  gumbel_formula, data = mbe1,
  stanvars = sv_gumbelbin,
  prior = gumbel_priors,
  init_r = 0.5
)

mbfit_e1_uvsd <- brm(
  uvsdt_formula, data = mbe1,
  stanvars = sv_uvsdtbin,
  prior = uvsdt_priors,
  init_r = 0.5
)

mbfit_e2_gumbel <- brm(
  gumbel_formula, data = mbe2,
  stanvars = sv_gumbelbin,
  prior = gumbel_priors,
  init_r = 0.5
)

mbfit_e2_uvsd <- brm(
  uvsdt_formula, data = mbe2,
  stanvars = sv_uvsdtbin,
  prior = uvsdt_priors,
  init_r = 0.5
)

mbfit_e3_gumbel <- brm(
  gumbel_formula, data = mbe3,
  stanvars = sv_gumbelbin,
  prior = gumbel_priors,
  init_r = 0.5
)

mbfit_e3_uvsd <- brm(
  uvsdt_formula, data = mbe3,
  stanvars = sv_uvsdtbin,
  prior = uvsdt_priors,
  init_r = 0.5
)


## loo fitting
all_dat %>% 
  group_by(experiment) %>% 
  summarise(n = n_distinct(Subject))


### k-fold CV in brms per default runs log_lik with:
### allow_new_levels = TRUE, sample_new_levels = "gaussian"
### gaussian means: If "gaussian", sample new levels from the (multivariate)
### normal distribution implied by the group-level standard deviations and
### correlations. see ?prepare_predictions
## it then aggregates the returned matrix xxx across samples using: 
## apply(xxx, 2, brms:::log_mean_exp)

library(future)
plan(multisession, workers = 12)
exloo_mbfit_e1_gumbel <- kfold(
  x = mbfit_e1_gumbel, group = "Subject", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_gumbelbin", "calc_posterior_predictions_gumbelbin", 
                                        "posterior_epred_gumbelbin", "posterior_predict_gumbelbin")))
# Based on 20-fold cross-validation.
# 
#            Estimate   SE
# elpd_kfold  -1016.3 26.2
# p_kfold       256.5 27.3
# kfoldic      2032.7 52.4
exloo_mbfit_e1_uvsd <- kfold(
  x = mbfit_e1_uvsd, group = "Subject", sample_new_levels = "uncertainty",
  future_args = list(future.globals = c("log_lik_uvsdtbin", "calc_posterior_predictions_uvsdtbin", 
                                        "posterior_epred_uvsdtbin", "posterior_predict_uvsdtbin")))
# Based on 20-fold cross-validation.
# 
#            Estimate   SE
# elpd_kfold  -1009.3 25.6
# p_kfold       299.3 26.0
# kfoldic      2018.6 51.2
exloo_bm_e1 <- loo_compare(exloo_mbfit_e1_gumbel, exloo_mbfit_e1_uvsd)
#                  elpd_diff se_diff
# mbfit_e1_uvsd    0.0       0.0   
# mbfit_e1_gumbel -7.0       8.8  

exloo_mbfit_e2_gumbel <- kfold(
  mbfit_e2_gumbel, group = "Subject", sample_new_levels = "uncertainty",
  future_args = list(future.globals = c("log_lik_gumbelbin", "calc_posterior_predictions_gumbelbin", 
                     "posterior_epred_gumbelbin", "posterior_predict_gumbelbin")))
#            Estimate   SE
# elpd_kfold   -835.9 16.9
# p_kfold       240.7 18.2
# kfoldic      1671.8 33.7
exloo_mbfit_e3_gumbel <- kfold(
  mbfit_e3_gumbel, group = "Subject", sample_new_levels = "uncertainty",
  future_args = list(future.globals = c("log_lik_gumbelbin", "calc_posterior_predictions_gumbelbin", 
                     "posterior_epred_gumbelbin", "posterior_predict_gumbelbin")))
#            Estimate   SE
# elpd_kfold   -635.9 16.0
# p_kfold       166.4 16.7
# kfoldic      1271.8 32.0
exloo_mbfit_e2_uvsd <- kfold(
  mbfit_e2_uvsd, group = "Subject", sample_new_levels = "uncertainty",
  future_args = list(future.globals = c("log_lik_uvsdtbin", "calc_posterior_predictions_uvsdtbin", 
                     "posterior_epred_uvsdtbin", "posterior_predict_uvsdtbin")))
#            Estimate   SE
# elpd_kfold   -831.0 17.9
# p_kfold       279.5 17.2
# kfoldic      1662.0 35.7
exloo_mbfit_e3_uvsd <- kfold(
  mbfit_e3_uvsd, group = "Subject", sample_new_levels = "uncertainty",
  future_args = list(future.globals = c("log_lik_uvsdtbin", "calc_posterior_predictions_uvsdtbin", 
                     "posterior_epred_uvsdtbin", "posterior_predict_uvsdtbin")))
#            Estimate   SE
# elpd_kfold   -640.2 16.7
# p_kfold       182.8 16.7
# kfoldic      1280.5 33.4

exloo_bm_e2 <- loo_compare(exloo_mbfit_e2_gumbel, exloo_mbfit_e2_uvsd)
#                 elpd_diff se_diff
# mbfit_e2_uvsd    0.0       0.0   
# mbfit_e2_gumbel -4.9       8.6   

exloo_bm_e3 <- loo_compare(exloo_mbfit_e3_gumbel, exloo_mbfit_e3_uvsd)
#                 elpd_diff se_diff
# mbfit_e3_gumbel  0.0       0.0   
# mbfit_e3_uvsd   -4.3       2.6   

# pptmp <- prepare_predictions(mbfit_e1_gumbel)
# log_lik_gumbelbin(2, pptmp)
# pptmp <- prepare_predictions(mbfit_e1_uvsd)

all_gumbel_fit <- list(
  mbfit_e1_gumbel, mbfit_e2_gumbel, mbfit_e3_gumbel
)
all_uvsd_fit <- list(
  mbfit_e1_uvsd, mbfit_e2_uvsd, mbfit_e3_uvsd
)

mb_loo_gumbel <- lapply(all_gumbel_fit, loo)
mb_waic_gumbel <- lapply(all_gumbel_fit, waic)

mb_loo_uvsdt <- lapply(all_uvsd_fit, loo)
mb_waic_uvsdt <- lapply(all_uvsd_fit, waic)

for (i in seq_along(mb_loo_gumbel)) {
  attr(mb_loo_gumbel[[i]], "model_name") <- "gumbel"
  attr(mb_waic_gumbel[[i]], "model_name") <- "gumbel"
  
  attr(mb_loo_uvsdt[[i]], "model_name") <- "uvsdt"
  attr(mb_waic_uvsdt[[i]], "model_name") <- "uvsdt"
  
}

waic_comp <- mapply(loo_compare, mb_waic_gumbel, mb_waic_uvsdt, SIMPLIFY = FALSE)
# [[1]]
#        elpd_diff se_diff
# uvsdt    0.0       0.0  
# gumbel -53.8      14.2  
# 
# [[2]]
#        elpd_diff se_diff
# uvsdt    0.0       0.0  
# gumbel -58.4      13.4  
# 
# [[3]]
#        elpd_diff se_diff
# uvsdt   0.0       0.0   
# gumbel -9.9       5.8 

loo_comp <- mapply(loo_compare, mb_loo_gumbel, mb_loo_uvsdt, SIMPLIFY = FALSE)
# [[1]]
#        elpd_diff se_diff
# uvsdt    0.0       0.0  
# gumbel -50.9      14.1  
# 
# [[2]]
#        elpd_diff se_diff
# uvsdt    0.0       0.0  
# gumbel -57.1      13.6  
# 
# [[3]]
#        elpd_diff se_diff
# uvsdt   0.0       0.0   
# gumbel -7.8       5.7   


# loo_gumbel_e1 <- loo(mbfit_e1_gumbel)
# loo_uvsd_e1 <- loo(mbfit_e1_uvsd)
# loo_compare(loo_gumbel_e1, loo_uvsd_e1)
# 
# loo_gumbel_e2 <- loo(mbfit_e2_gumbel)
# loo_uvsd_e2 <- loo(mbfit_e2_uvsd)
# loo_compare(loo_gumbel_e2, loo_uvsd_e2)
# 
# loo_gumbel_e3 <- loo(mbfit_e3_gumbel)
# loo_uvsd_e3 <- loo(mbfit_e3_uvsd)
# loo_compare(loo_gumbel_e3, loo_uvsd_e3)

### make plots

pred_gumbel <- lapply(all_gumbel_fit, posterior_epred)
pred_uvsd <- lapply(all_uvsd_fit, posterior_epred)

plot_dat <- all_dat %>% 
  group_by(experiment, Subject, BaseRate) %>% 
  summarise(
    hit = OLD_old/Nold,
    fa = NEW_old/Nnew,
  ) 

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

plot_dat <- plot_dat %>% 
  group_by(experiment, BaseRate) %>% 
  summarise(across(c(hit,fa, hit_gumbel, fa_gumbel, hit_uvsd, fa_uvsd), mean))

plot_data_mb <- plot_dat

save(plot_data_mb, 
     exloo_mbfit_e1_uvsd, exloo_mbfit_e1_gumbel, 
     exloo_mbfit_e2_uvsd, exloo_mbfit_e2_gumbel, 
     exloo_mbfit_e3_uvsd, exloo_mbfit_e3_gumbel, 
     exloo_bm_e1, exloo_bm_e2, exloo_bm_e3,
     file = "mb_exloo_res.rda")

psize <- 3.5
lsize <- 1.5
plot_dat %>%
  ggplot(aes(x =  fa, y = hit)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize) +
  geom_point(aes(x = fa_gumbel, y = hit_gumbel), shape = 3, colour = "blue", size = psize) +
  geom_point(aes(x = fa_uvsd, y = hit_uvsd), shape = 2, colour = "red", size = psize) + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(experiment)) +
  labs(x = "False alarms", y = "Hits")
ggsave("malejka-broder-plot1.pdf", width = 22, height = 10, units = "cm")
