
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


roc6_loo_gumbel <- lapply(roc6_fits_gumbel, loo)
roc6_waic_gumbel <- lapply(roc6_fits_gumbel, waic)

roc6_loo_uvsdt <- lapply(roc6_fits_uvsdt, loo)
roc6_waic_uvsdt <- lapply(roc6_fits_uvsdt, waic)

for (i in seq_along(roc6_loo_gumbel)) {
  attr(roc6_loo_gumbel[[i]], "model_name") <- "gumbel"
  attr(roc6_waic_gumbel[[i]], "model_name") <- "gumbel"
  
  attr(roc6_loo_uvsdt[[i]], "model_name") <- "uvsdt"
  attr(roc6_waic_uvsdt[[i]], "model_name") <- "uvsdt"
  
}

waic_comp <- mapply(loo_compare, roc6_waic_gumbel, roc6_waic_uvsdt, SIMPLIFY = FALSE)
loo_comp <- mapply(loo_compare, roc6_loo_gumbel, roc6_loo_uvsdt, SIMPLIFY = FALSE)

str(roc6_waic_gumbel[[1]])

roc6_waic_gumbel[[1]]$estimates["waic","Estimate"]

#mod_comp <- 
tibble(
  dataset = dataset6,
  waic_g = map_dbl(roc6_waic_gumbel, ~.$estimates["waic","Estimate"]),
  waic_uv = map_dbl(roc6_waic_uvsdt, ~.$estimates["waic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(waic_g, waic_uv)) %>% 
  mutate(across(c(waic_g, waic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(waic_diff_sig = map_lgl(waic_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))

tibble(
  dataset = dataset6,
  looic_g = map_dbl(roc6_loo_gumbel, ~.$estimates["looic","Estimate"]),
  looic_uv = map_dbl(roc6_loo_uvsdt, ~.$estimates["looic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(looic_g, looic_uv)) %>% 
  mutate(across(c(looic_g, looic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(waic_diff_sig = map_lgl(loo_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))

## plots
plot_data <- roc6 %>% 
  group_by(exp) %>% 
  summarise(across(c(OLD_3new:NEW_3old), sum)) %>% 
  pivot_longer(-exp, names_to = c("status", "response"), names_sep = "_") %>% 
  group_by(exp, status) %>% 
  mutate(observed = value / sum(value)) %>% 
   mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("3new", "2new", "1new", 
                                           "1old", "2old", "3old")))
  
pred_gumbel <- lapply(roc6_fits_gumbel, posterior_epred)
pred_uvsdt <- lapply(roc6_fits_uvsdt, posterior_epred)

plot_data$gumbel <- unlist(map(pred_gumbel, ~apply(., c(3), mean)))
plot_data$uvsd <- unlist(map(pred_uvsdt, ~apply(., c(3), mean)))

plot_data2 <- plot_data
plot_data2$gumbel_low <- unlist(map(pred_gumbel, 
                                   ~apply(apply(., c(1, 3), mean), 
                                          2, quantile, probs = 0.025)))
plot_data2$gumbel_high <- unlist(map(pred_gumbel, 
                                    ~apply(apply(., c(1, 3), mean), 
                                           2, quantile, probs = 0.975)))

plot_data2$uvsd_low <- unlist(map(pred_uvsdt, 
                                   ~apply(apply(., c(1, 3), mean), 
                                          2, quantile, probs = 0.025)))
plot_data2$uvsd_high <- unlist(map(pred_uvsdt, 
                                    ~apply(apply(., c(1, 3), mean), 
                                           2, quantile, probs = 0.975)))

# str(apply(pred_gumbel[[1]], c(1, 3), mean))
# str(pred_gumbel)

plot_data %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "3new") %>% 
  ggplot(aes(x =  observed_NEW, y = observed_OLD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_point(aes(x = gumbel_NEW, y = gumbel_OLD), shape = 3, colour = "blue") +
  geom_point(aes(x = uvsd_NEW, y = uvsd_OLD), shape = 2, colour = "red") + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(vars(exp))


plot_data2 %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd, 
                                                   gumbel_low, gumbel_high, 
                                                   uvsd_low, uvsd_high)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "3new") %>% 
  ggplot(aes(x =  observed_NEW, y = observed_OLD)) +
  geom_abline(slope = 1, intercept = 0) +
   
    geom_linerange(aes(x = uvsd_NEW, y = uvsd_OLD, ymin = uvsd_low_OLD, ymax = uvsd_high_OLD), colour = "red") +
   geom_linerange(aes(x = uvsd_NEW, y = uvsd_OLD, xmin = uvsd_low_NEW, xmax = uvsd_high_NEW), colour = "red") +
  geom_linerange(aes(x = gumbel_NEW, y = gumbel_OLD, ymin = gumbel_low_OLD, ymax = gumbel_high_OLD), colour = "blue") +
   geom_linerange(aes(x = gumbel_NEW, y = gumbel_OLD, xmin = gumbel_low_NEW, xmax = gumbel_high_NEW), colour = "blue") +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_point(aes(x = uvsd_NEW, y = uvsd_OLD), shape = 2, colour = "red") +
  geom_point(aes(x = gumbel_NEW, y = gumbel_OLD), shape = 4, colour = "blue") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(vars(exp))

