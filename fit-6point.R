
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel6agg-stan.R")
source("uvsdt6agg-stan.R")

data("roc6", package = "MPTinR")
head(roc6)

dataset6 <- levels(roc6$exp)

gumbel_formula <- brmsformula(
  OLD_3new ~ 1 + (1|p|id), 
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
  prior(student_t(3, -1, 2), class = Intercept)

uvsdt_formula <- brmsformula(
  OLD_3new ~ 1 + (1|p|id), 
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
  prior(student_t(3, 1, 2), class = Intercept)

roc6_data <- vector("list", length(dataset6))
roc6_oldmat <- vector("list", length(dataset6))
roc6_newmat <- vector("list", length(dataset6))
roc6_sv <- vector("list", length(dataset6))

roc6_fits_gumbel <- vector("list", length(dataset6))
roc6_fits_uvsdt <- vector("list", length(dataset6))


for (i in seq_along(dataset6)) {
  print(i)
  roc6_data[[i]] <- d6_1 <- roc6 %>% 
    filter(exp == dataset6[i])
  roc6_oldmat[[i]] <- roc6_data[[i]] %>% 
                                  select(OLD_3new:OLD_3old) %>% 
                                  as.matrix()
  roc6_newmat[[i]] <- roc6_data[[i]] %>% 
                                  select(NEW_3new:NEW_3old) %>% 
                                  as.matrix()
  roc6_sv[[i]] <- stanvar(roc6_oldmat[[i]], name = "oldmat", block = "data") +
    stanvar(scode = "array[N, 6] int oldmat2;", block = "tdata") +
    stanvar(scode = "oldmat2 = to_int(to_array_2d(oldmat));", block = "tdata") +
    stanvar(roc6_newmat[[i]], name = "newmat", block = "data") +
    stanvar(scode = "array[N, 6] int newmat2;", block = "tdata") +
    stanvar(scode = "newmat2 = to_int(to_array_2d(newmat));", block = "tdata")
  
  roc6_fits_gumbel[[i]] <- brm(
    gumbel_formula, data = roc6_data[[i]], 
    stanvars = sv_gumbel6agg + roc6_sv[[i]], 
    prior = gumbel_priors,
    init_r = 0.5
  )
  
  roc6_fits_uvsdt[[i]] <- brm(
    uvsdt_formula, data = roc6_data[[i]], 
    stanvars = sv_uvsdt6agg + roc6_sv[[i]], 
    prior = uvsdt_priors,
    init_r = 0.5
  )
}

# stancode(gumbel_formula, data = roc6_data[[i]], 
#          stanvars = roc6_sv[[i]], 
#          prior = gumbel_priors)


roc6_loo_gumbel <- lapply(roc6_fits, loo)
roc6_waic_gumbel <- lapply(roc6_fits, waic)

