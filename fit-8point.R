
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel8agg-stan.R")
source("uvsdt8agg-stan.R")
source("gumbel8agglog-stan.R")
source("uvsdt8agglog-stan.R")

sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "id")

data("roc8", package = "MPTinR")
head(roc8)

dataset8 <- levels(roc8$exp)

all_perf <- roc8 %>% 
  mutate(
    hit = rowSums(select(., OLD_1old:OLD_4old)) /
           rowSums(select(., OLD_4new:OLD_4old)),
    fa = rowSums(select(., NEW_1old:NEW_4old)) /
           rowSums(select(., NEW_4new:NEW_4old))
  ) %>% 
  mutate(
    acc = (hit + (1-fa)) / 2
  ) %>% 
  mutate(
    empty = rowSums(select(., OLD_4new:NEW_4old) == 0)
  )

all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  arrange(desc(empty))

# all_perf %>% 
#   filter(exp == dataset8[1]) %>% 
#   arrange(desc(fa))
# 
# all_perf %>% 
#   filter(exp == dataset8[1]) %>% 
#   arrange(hit)
# 
# all_perf %>% 
#   filter(exp == dataset8[1]) %>% 
#   arrange(desc(acc))
# 
# all_perf %>% 
#   filter(exp == dataset8[1]) %>% 
#   arrange(acc)

low_perf <- all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  filter(acc < .59) ## select only significant above chance
  #filter(acc < .65 | empty > 6)
  #filter(empty > 6)

roc8_use <- roc8 %>% 
  filter(!(id %in% low_perf$id))

gumbel_formula_8 <- brmsformula(
  OLD_4new | vint(OLD_3new, OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, OLD_4old, NEW_4new, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old, NEW_4old) ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = gumbel8agglog_family
)

gumbel_formula_8_nolog <- brmsformula(
  OLD_4new | vint(OLD_3new, OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, OLD_4old, NEW_4new, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old, NEW_4old) ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = gumbel8agg_family
)


gumbel_priors_8 <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlx") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhx") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

uvsdt_formula_8 <- brmsformula(
  OLD_4new | vint(OLD_3new, OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, OLD_4old, NEW_4new, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old, NEW_4old) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = uvsdt8agglog_family
)

uvsdt_formula_8_nolog <- brmsformula(
  OLD_4new | vint(OLD_3new, OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, OLD_4old, NEW_4new, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old, NEW_4old) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = uvsdt8agg_family
)

uvsdt_priors_8 <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlx") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhx") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

roc8_data <- vector("list", length(dataset8))

roc8_fits_gumbel <- vector("list", length(dataset8))
roc8_fits_uvsdt <- vector("list", length(dataset8))

#i <- 1
# iter <- 4000
# warmup <- 1000

for (i in seq_along(dataset8)) {
  print(i)
  roc8_data[[i]] <- roc8_use %>% 
    filter(exp == dataset8[i])
  
  roc8_fits_gumbel[[i]] <- brm(
    gumbel_formula_8, data = roc8_data[[i]], 
    stanvars = sv_gumbel8agglog, 
    prior = gumbel_priors_8,
    iter = iter, warmup = warmup,
    init_r = 0.25, 
    #control = list(adapt_delta = 0.999999, max_treedepth = 20)
  )
  roc8_fits_uvsdt[[i]] <- brm(
    uvsdt_formula_8, data = roc8_data[[i]], 
    stanvars = sv_uvsdt8agglog, 
    prior = uvsdt_priors_8,
    iter = iter, warmup = warmup,
    init_r = 0.5, 
    #control = list(adapt_delta = 0.999999, max_treedepth = 20)
  )

  # roc8_fits_uvsdt[[i]] <- brm(
  #   uvsdt_formula_8_nolog, data = roc8_data[[i]], 
  #   stanvars = sv_uvsdt8agg, 
  #   prior = uvsdt_priors_8,
  #   iter = iter, warmup = warmup,
  #   init_r = 0.5, 
  #   #control = list(adapt_delta = 0.999999, max_treedepth = 20)
  # )
}
# xxx <- map(roc8_fits_gumbel, ~rstan::get_sampler_params(.$fit))
# for (i in seq_along(dataset8)) {
#   cat(dataset8[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
# }
# 
# xxy <- map(roc8_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
# for (i in seq_along(dataset8)) {
#   cat(dataset8[i], ": ", sum(map_dbl(xxy[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
# }



## plots
