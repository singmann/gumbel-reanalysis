
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
#                                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# sd(Intercept)                           0.66      0.11     0.48     0.90 1.01      920
# sd(discsignal_Intercept)                0.09      0.05     0.01     0.20 1.00     1034
# cor(Intercept,discsignal_Intercept)     0.48      0.41    -0.64     0.97 1.00     2360
#                                     Tail_ESS
# sd(Intercept)                           1654
# sd(discsignal_Intercept)                1629
# cor(Intercept,discsignal_Intercept)     2232
# 
# Regression Coefficients:
#                      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                1.29      0.13     1.05     1.56 1.00      697     1136
# discsignal_Intercept     0.36      0.04     0.28     0.45 1.00     2298     2515

uvsdt_pred <- posterior_epred(fit_uvsdt_rank)
str(gumbel_pred)

## UVSDT predictions:
up2 <- apply(uvsdt_pred, c(1,3), mean)
apply(up2, 2, mean) %>% 
  print(digits = 2)
#   R1   R2   R3   R4 
# 0.60 0.18 0.12 0.11 
apply(up2, 2, quantile, probs = c(0.027, 0.975)) %>% 
  print(digits = 2)
#         R1   R2   R3    R4
# 2.7%  0.58 0.17 0.11 0.097
# 97.5% 0.61 0.19 0.12 0.117

source("gumbelrank2-stan.R")
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
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )
# Multilevel Hyperparameters:
# ~id (Number of levels: 30) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.57      0.09     0.42     0.78 1.00      867     1766
# 
# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     1.21      0.11     1.00     1.44 1.00      594      903

gumbel_pred <- posterior_epred(fit_gumbel_rank)
str(gumbel_pred)

## Gumbel predictions:
gp2 <- apply(gumbel_pred, c(1,3), mean)
apply(gp2, 2, mean) %>% 
  print(digits = 2)
#   R1   R2   R3   R4 
# 0.60 0.18 0.12 0.10 
apply(gp2, 2, quantile, probs = c(0.027, 0.975)) %>% 
  print(digits = 2)
#         R1   R2   R3    R4
# 2.7%  0.58 0.17 0.12 0.094
# 97.5% 0.61 0.19 0.13 0.106

## mean of data
(colSums(ranking.data[,-5]) / sum(ranking.data[,-5])) %>% 
  print(digits = 2)
## data:
#   V1   V2   V3   V4 
# 0.60 0.16 0.14 0.10 


### loo
lu <- loo(fit_uvsdt_rank)
lg <- loo(fit_gumbel_rank)
loo_compare(lu, lg)
#                 elpd_diff se_diff
# fit_gumbel_rank  0.0       0.0   
# fit_uvsdt_rank  -3.1       3.6   

