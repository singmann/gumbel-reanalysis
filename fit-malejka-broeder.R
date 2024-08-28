
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

plot_dat %>%
  ggplot(aes(x =  fa, y = hit)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_point(aes(x = fa_gumbel, y = hit_gumbel), shape = 3, colour = "blue") +
  geom_point(aes(x = fa_uvsd, y = hit_uvsd), shape = 2, colour = "red") + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(vars(experiment))
