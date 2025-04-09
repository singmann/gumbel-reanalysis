
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel6agg-stan.R")
source("uvsdt6agg-stan.R")

sd_priors <- set_prior("student_t(5, 0, 2.5)", class = "sd", group = "id")

data("roc6", package = "MPTinR")
head(roc6)

roc6_use <- roc6

str(roc6_use)

dataset6 <- levels(roc6_use$exp)

gumbel_formula <- brmsformula(
  OLD_3new | vint(OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old) ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id),
  family = gumbel6agg_family, cmc = FALSE
)

gumbel_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

uvsdt_formula <- brmsformula(
  OLD_3new | vint(OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id),
  family = uvsdt6agg_family, cmc = FALSE
)

uvsdt_priors <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept) +
  sd_priors

roc6_data <- vector("list", length(dataset6))

roc6_fits_gumbel <- vector("list", length(dataset6))
roc6_fits_uvsdt <- vector("list", length(dataset6))

control1 <- list(adapt_delta = 0.99, max_treedepth = 20)
control2 <- list(adapt_delta = 0.9999999, max_treedepth = 20)

for (i in seq_along(dataset6)) {
  print(i)
  roc6_data[[i]] <- roc6_use %>% 
    filter(exp == dataset6[i])
  roc6_fits_gumbel[[i]] <- brm(
    gumbel_formula, data = roc6_data[[i]], 
    stanvars = sv_gumbel6agg ,
    prior = gumbel_priors,
    init_r = 0.5, 
    control = if (i %in% c(5, 6, 7)) control2 else control1
  )

  roc6_fits_uvsdt[[i]] <- brm(
    uvsdt_formula, data = roc6_data[[i]], 
    stanvars = sv_uvsdt6agg, 
    prior = uvsdt_priors,
    init_r = 0.5, 
    control = control1
  )

}

xxx <- map(roc6_fits_gumbel, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset6)) {
  cat(dataset6[i], ": ", sum(map_dbl(xxx[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
}

xxy <- map(roc6_fits_uvsdt, ~rstan::get_sampler_params(.$fit))
for (i in seq_along(dataset6)) {
  cat(dataset6[i], ": ", sum(map_dbl(xxy[[i]], ~sum(.[1001:2000,"divergent__"]))), "\n")
}


###########
