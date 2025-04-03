
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
  prior(student_t(3, 1, 2), class = Intercept)

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
  prior(student_t(3, 1, 2), class = Intercept)

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


# data("roc8", package = "MPTinR")
# head(roc8)
# 
# benjamin_6p <- roc8 %>% 
#   filter(exp == "Benjamin_2013") %>% 
#   mutate(
#     OLD_2new_new = OLD_3new + OLD_2new,
#     OLD_2old_new = OLD_3old + OLD_2old,
#     NEW_2new_new = NEW_3new + NEW_2new,
#     NEW_2old_new = NEW_3old + NEW_2old
#   ) %>% 
#   mutate(
#     OLD_3new = OLD_4new, 
#     OLD_3old = OLD_4old,
#     NEW_3new = NEW_4new, 
#     NEW_3old = NEW_4old
#   ) %>% 
#   mutate(
#     OLD_2new = OLD_2new_new, 
#     OLD_2old = OLD_2old_new,
#     NEW_2new = NEW_2new_new, 
#     NEW_2old = NEW_2old_new
#   ) %>% 
#   select(OLD_3new:OLD_3old, NEW_3new:NEW_3old, exp, id)
# 
# roc8 %>%
#   filter(exp == "Benjamin_2013") %>%
#   select(-exp, -id) %>%
#   rowSums()
# 
# benjamin_6p %>% 
#   select(-exp, -id) %>% 
#   rowSums()
# 
# benjamin_6p %>% 
#   mutate(
#     hit = rowSums(cbind(OLD_1old, OLD_2old, OLD_3old)) / 
#       rowSums(cbind(OLD_3new, OLD_2new, OLD_1new, 
#                     OLD_1old, OLD_2old, OLD_3old)),
#     fa = rowSums(cbind(NEW_1old, NEW_2old, NEW_3old)) / 
#       rowSums(cbind(NEW_3new, NEW_2new, NEW_1new, 
#                     NEW_1old, NEW_2old, NEW_3old))
#   ) %>% 
#   mutate(acc = (hit + (1-fa)) / 2) %>% 
#   arrange(acc) %>% 
#   filter(acc < .6) %>% 
#   select(id) %>% 
#   unlist() %>% 
#   unname() %>% 
#   as.character() %>% 
#   dput()
# 
# benjamin_6p <- benjamin_6p %>% 
#   filter(!(id %in% c("74:Benjamin", "91:Benjamin", "32:Benjamin", "35:Benjamin", 
#                      "88:Benjamin", "109:Benjamin", "67:Benjamin", "84:Benjamin", 
#                      "108:Benjamin", "61:Benjamin", "64:Benjamin", "110:Benjamin")))
# 
# roc6_use <- bind_rows(
#   benjamin_6p, roc6
# )
