
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel6agg-stan.R")

data("roc6", package = "MPTinR")
head(roc6)

dataset6 <- levels(roc6$exp)

d6_1 <- roc6 %>% 
  filter(exp == dataset6[1])

d6_1_oldmat <- d6_1 %>% 
  select(OLD_3new:OLD_3old) %>% 
  as.matrix()

d6_1_newmat <- d6_1 %>% 
  select(NEW_3new:NEW_3old) %>% 
  as.matrix()

gumbel_formula <- brmsformula(
  OLD_3new ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id),
  family = gumbel6agg_family, cmc = FALSE
)

get_prior(gumbel_formula, data = d6_1)

gumbel_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(student_t(3, 0, 2), class = Intercept)
sv_gumbel6agg_1 <- sv_gumbel6agg +
  stanvar(d6_1_oldmat, name = "oldmat", block = "data") +
  stanvar(scode = "array[N, 6] int oldmat2;", block = "tdata") +
  stanvar(scode = "oldmat2 = to_int(to_array_2d(oldmat));", block = "tdata") +
  stanvar(d6_1_newmat, name = "newmat", block = "data") +
  stanvar(scode = "array[N, 6] int newmat2;", block = "tdata") +
  stanvar(scode = "newmat2 = to_int(to_array_2d(newmat));", block = "tdata")


stancode(gumbel_formula, data = d6_1, 
         stanvars = sv_gumbel6agg_1, 
         prior = gumbel_priors)

f6_1_gumbel_1 <- brm(
  gumbel_formula, data = d6_1, 
  stanvars = sv_gumbel6agg_1, 
  prior = gumbel_priors,
  file = "f6e1_gumbel_1", init_r = 0.5, 
  save_pars = save_pars(all = TRUE)
)
f6_1_gumbel_1

# pptmp <- prepare_predictions(f6_1_gumbel_1)
# log_lik_gumbel6agg(2, pptmp)
#posterior_epred_gumbel6agg(pptmp)
posterior_predict_gumbel6agg(20, pptmp)
#loo_e1_gumbel <- loo(f6_1_gumbel_1)

preds <- posterior_epred(f6_1_gumbel_1)
str(preds)

preds[1:5, 1, ]
rowSums(preds[1:5, 1, ])

ppreds <- posterior_predict(f6_1_gumbel_1)
str(ppreds)

do.call("rbind", ppreds[1:10, 27])
