
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
load("malejka-broeder.rda")
source("gumbelbin-stan.R")
source("uvsdtbin-stan.R")

mbe1 <- mbe1 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old)

mbe2 <- mbe2 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old)

mbe3 <- mbe3 %>% 
  mutate(Nold = OLD_new + OLD_old,
         Nnew = NEW_new + NEW_old)

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


# pptmp <- prepare_predictions(mbfit_e1_gumbel)
# log_lik_gumbelbin(2, pptmp)
# pptmp <- prepare_predictions(mbfit_e1_uvsd)

loo_gumbel_e1 <- loo(mbfit_e1_gumbel)
loo_uvsd_e1 <- loo(mbfit_e1_uvsd)
loo_compare(loo_gumbel_e1, loo_uvsd_e1)

loo_gumbel_e2 <- loo(mbfit_e2_gumbel)
loo_uvsd_e2 <- loo(mbfit_e2_uvsd)
loo_compare(loo_gumbel_e2, loo_uvsd_e2)

loo_gumbel_e3 <- loo(mbfit_e3_gumbel)
loo_uvsd_e3 <- loo(mbfit_e3_uvsd)
loo_compare(loo_gumbel_e3, loo_uvsd_e3)
