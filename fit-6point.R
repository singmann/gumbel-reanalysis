
library("tidyverse")
library("brms")
load("dat-prep.rda")

source("gumbel6agg-stan.R")

dataset6 <- levels(d6$exp)

d6_1 <- d6 %>% 
  filter(exp == dataset6[1])


d6_1_dmat <- d6_1 %>% 
  select(r3new:r3old) %>% 
  as.matrix()

gumbel_formula <- brmsformula(
  r3new | vint(statusnum)  ~ 1 + (1|p|id), 
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
sv_gumbel6agg <- sv_gumbel6agg +
  stanvar(d6_1_dmat, name = "dmat", block = "data")

stancode(gumbel_formula, data = d6_1, 
         stanvars = sv_gumbel6agg, 
         prior = gumbel_priors)
