library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

source("data_from_david.R")

source("gumbelrank2-stan.R")
source("uvsdtrank-stan.R")

gumbel_priors <- prior(student_t(3, 1, 2), class = Intercept)
uvsdt_priors <- prior(student_t(3, 0.5, 1), 
                      class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept)

##---------------------------------------------------------------
##                            4 Ranks                           -
##---------------------------------------------------------------

##------------
##  KKS 2012  
##------------

gumbel_formula_kks <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

# stancode(gumbel_formula_kks, data = kks12, 
#          stanvars = sv_gumbelrank, prior = gumbel_priors)

fit_kks_gumbel <- brm(
    gumbel_formula_kks, data = kks12, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kks <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
)

# stancode(uvsdt_formula_kks, data = kks12, 
#          stanvars = sv_uvsdtrank, prior = uvsdt_priors)

fit_kks_uvsdt <- brm(
    uvsdt_formula_kks, data = kks12, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.5
  )

kks12


##-----------
##  KK14 E1  
##-----------

kk14_e1_use <- kk14_e1 %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

gumbel_formula_kke1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

fit_kke1_gumbel <- brm(
    gumbel_formula_kke1, data = kk14_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kke1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
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


mhe_e1_use <- mhe_e1 %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

gumbel_formula_mhe1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

fit_mhe1_gumbel <- brm(
    gumbel_formula_mhe1, data = mhe_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_mhe1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
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


##-------------
##  Exact LOO  
##-------------

library(future)
plan(multisession, workers = 12)

exloo_kks_gumbel <- kfold(
  x = fit_kks_gumbel, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_gumbelrank", "calc_posterior_predictions_gumbelrank", 
                                        "posterior_epred_gumbelrank", "posterior_predict_gumbelrank")))

exloo_kke1_gumbel <- kfold(
  x = fit_kke1_gumbel, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_gumbelrank", "calc_posterior_predictions_gumbelrank", 
                                        "posterior_epred_gumbelrank", "posterior_predict_gumbelrank")))

exloo_mhe1_gumbel <- kfold(
  x = fit_mhe1_gumbel, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_gumbelrank", "calc_posterior_predictions_gumbelrank", 
                                        "posterior_epred_gumbelrank", "posterior_predict_gumbelrank")))

exloo_kks_uvsdt <- kfold(
  x = fit_kks_uvsdt, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_uvsdtrank", "calc_posterior_predictions_uvsdtrank", 
                                        "posterior_epred_uvsdtrank", "posterior_predict_uvsdtrank")))

exloo_kke1_uvsdt <- kfold(
  x = fit_kke1_uvsdt, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_uvsdtrank", "calc_posterior_predictions_uvsdtrank", 
                                        "posterior_epred_uvsdtrank", "posterior_predict_uvsdtrank")))

exloo_mhe1_uvsdt <- kfold(
  x = fit_mhe1_uvsdt, group = "id", sample_new_levels = "uncertainty", 
  future_args = list(future.globals = c("log_lik_uvsdtrank", "calc_posterior_predictions_uvsdtrank", 
                                        "posterior_epred_uvsdtrank", "posterior_predict_uvsdtrank")))


exloo_kks <- loo_compare(exloo_kks_uvsdt, exloo_kks_gumbel)
exloo_kke1 <- loo_compare(exloo_kke1_uvsdt, exloo_kke1_gumbel)
exloo_mhe1 <- loo_compare(exloo_mhe1_uvsdt, exloo_mhe1_gumbel)

save(exloo_kke1, exloo_kke1_gumbel, exloo_kke1_uvsdt, 
     exloo_kks, exloo_kks_gumbel, exloo_kks_uvsdt, 
     exloo_mhe1, exloo_mhe1_gumbel, exloo_mhe1_uvsdt, file = "exloo_4rank.rda")


##---------------------------------------------------------------
##                            3 Ranks                           -
##---------------------------------------------------------------

kk14_e2
mg16_e1
mg16_e2