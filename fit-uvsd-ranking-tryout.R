
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
source("uvsdtrank-stan.R")

### data
?MPTinR::fit.mptinr
ranking.data <- structure(c(39, 80, 75, 35, 61, 54, 73, 52, 44, 63, 40, 48, 80,
49, 43, 80, 68, 53, 81, 60, 60, 65, 49, 58, 69, 75, 71, 47, 44,
85, 23, 9, 11, 21, 12, 21, 14, 20, 19, 15, 29, 13, 14, 15, 22,
11, 12, 16, 13, 20, 20, 9, 26, 19, 13, 9, 14, 15, 24, 9, 19,
7, 9, 26, 16, 14, 6, 17, 21, 14, 20, 18, 5, 19, 17, 5, 11, 21,
4, 9, 15, 17, 7, 17, 11, 11, 9, 19, 20, 3, 19, 4, 5, 18, 11,
11, 7, 11, 16, 8, 11, 21, 1, 17, 18, 4, 9, 10, 2, 11, 5, 9, 18,
6, 7, 5, 6, 19, 12, 3), .Dim = c(30L, 4L)) %>% 
  as.data.frame() %>% 
  mutate(id = 1:n())

uvsdt_formula <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
)

get_prior(uvsdt_formula, data = ranking.data)

uvsdt_priors <- prior(student_t(3, 0.5, 1), 
                      class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept)

stancode(uvsdt_formula, data = ranking.data, 
         stanvars = sv_uvsdtrank, prior = uvsdt_priors)


fit_uvsdt_rank <- brm(
    uvsdt_formula, data = ranking.data, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.5
  )
# Multilevel Hyperparameters:
# ~id (Number of levels: 30) 
#                                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)                           0.66      0.11     0.47     0.91 1.00      836     1265
# sd(discsignal_Intercept)                0.09      0.05     0.00     0.21 1.00      996     1446
# cor(Intercept,discsignal_Intercept)    -0.48      0.40    -0.97     0.56 1.00     2536     2224
# 
# Regression Coefficients:
#                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                1.29      0.13     1.04     1.54 1.01      579     1064
# discsignal_Intercept    -0.36      0.04    -0.44    -0.28 1.00     2371     2725


gumbel_formula <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

get_prior(gumbel_formula, data = ranking.data)

gumbel_priors <- prior(student_t(3, 1, 2), class = Intercept)

stancode(gumbel_formula, data = ranking.data, 
         stanvars = sv_gumbelrank, prior = gumbel_priors)

fit_gumbel_rank <- brm(
    gumbel_formula, data = ranking.data, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25, 
    control = list(adapt_delta = 0.99)
  )
