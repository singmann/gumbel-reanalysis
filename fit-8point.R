
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
#load("dat-prep.rda")
source("gumbel8agg-stan.R")
source("uvsdt8agg-stan.R")

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

all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  arrange(desc(fa))

all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  arrange(hit)

all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  arrange(desc(acc))

all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  arrange(acc)

low_perf <- all_perf %>% 
  filter(exp == dataset8[1]) %>% 
  filter(acc < .59)
  #filter(acc < .65 | empty > 6)
  #filter(empty > 6)

roc8_use <- roc8 %>% 
  filter(!(id %in% low_perf$id))

gumbel_formula_8 <- brmsformula(
  OLD_3new ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = gumbel8agg_family, cmc = FALSE
)

gumbel_priors_8 <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlx") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhx") +
  prior(student_t(3, -1, 2), class = Intercept)

uvsdt_formula_8 <- brmsformula(
  OLD_3new ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id), 
  crc ~ (1|p|id), 
  crlm ~ (1|p|id), crll ~ (1|p|id), crlx ~ (1|p|id), 
  crhm ~ (1|p|id), crhh ~ (1|p|id), crhx ~ (1|p|id),
  family = uvsdt8agg_family, cmc = FALSE
)

uvsdt_priors_8 <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlx") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhx") +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept)

roc8_data <- vector("list", length(dataset8))
roc8_oldmat <- vector("list", length(dataset8))
roc8_newmat <- vector("list", length(dataset8))
roc8_sv <- vector("list", length(dataset8))

roc8_fits_gumbel <- vector("list", length(dataset8))
roc8_fits_uvsdt <- vector("list", length(dataset8))

#i <- 1
for (i in seq_along(dataset8)) {
  print(i)
  roc8_data[[i]] <- roc8_use %>% 
    filter(exp == dataset8[i])
  roc8_oldmat[[i]] <- roc8_data[[i]] %>% 
                                  select(OLD_4new:OLD_4old) %>% 
                                  as.matrix()
  roc8_newmat[[i]] <- roc8_data[[i]] %>% 
                                  select(NEW_4new:NEW_4old) %>% 
                                  as.matrix()
  roc8_sv[[i]] <- stanvar(roc8_oldmat[[i]], name = "oldmat", block = "data") +
    stanvar(scode = "array[N, 8] int oldmat2;", block = "tdata") +
    stanvar(scode = "oldmat2 = to_int(to_array_2d(oldmat));", block = "tdata") +
    stanvar(roc8_newmat[[i]], name = "newmat", block = "data") +
    stanvar(scode = "array[N, 8] int newmat2;", block = "tdata") +
    stanvar(scode = "newmat2 = to_int(to_array_2d(newmat));", block = "tdata")
  
  roc8_fits_gumbel[[i]] <- brm(
    gumbel_formula_8, data = roc8_data[[i]], 
    stanvars = sv_gumbel8agg + roc8_sv[[i]], 
    prior = gumbel_priors_8,
    init_r = 0.25, control = list(adapt_delta = 0.9999)
  )
  
  roc8_fits_uvsdt[[i]] <- brm(
    uvsdt_formula_8, data = roc8_data[[i]], 
    stanvars = sv_uvsdt8agg + roc8_sv[[i]], 
    prior = uvsdt_priors_8,
    init_r = 0.5, control = list(adapt_delta = 0.9999)
  )
}

# stancode(gumbel_formula, data = roc6_data[[i]], 
#          stanvars = roc6_sv[[i]], 
#          prior = gumbel_priors)


## exclude data set 1 which sows problems

roc8_loo_gumbel <- lapply(roc8_fits_gumbel, loo)
roc8_waic_gumbel <- lapply(roc8_fits_gumbel, waic)

roc8_loo_uvsdt <- lapply(roc8_fits_uvsdt, loo)
roc8_waic_uvsdt <- lapply(roc8_fits_uvsdt, waic)

for (i in seq_along(roc8_loo_gumbel)) {
  attr(roc8_loo_gumbel[[i]], "model_name") <- "gumbel"
  attr(roc8_waic_gumbel[[i]], "model_name") <- "gumbel"
  
  attr(roc8_loo_uvsdt[[i]], "model_name") <- "uvsdt"
  attr(roc8_waic_uvsdt[[i]], "model_name") <- "uvsdt"
  
}

waic_comp <- mapply(loo_compare, roc8_waic_gumbel, roc8_waic_uvsdt, SIMPLIFY = FALSE)
loo_comp <- mapply(loo_compare, roc8_loo_gumbel, roc8_loo_uvsdt, SIMPLIFY = FALSE)

str(roc8_waic_gumbel[[1]])

roc8_waic_gumbel[[1]]$estimates["waic","Estimate"]

#mod_comp <- 
tibble(
  dataset = dataset8[-1],
  waic_g = map_dbl(roc8_waic_gumbel, ~.$estimates["waic","Estimate"]),
  waic_uv = map_dbl(roc8_waic_uvsdt, ~.$estimates["waic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(waic_g, waic_uv)) %>% 
  mutate(across(c(waic_g, waic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(waic_diff_sig = map_lgl(waic_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))

tibble(
  dataset = dataset8[-1],
  looic_g = map_dbl(roc8_loo_gumbel, ~.$estimates["looic","Estimate"]),
  looic_uv = map_dbl(roc8_loo_uvsdt, ~.$estimates["looic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(looic_g, looic_uv)) %>% 
  mutate(across(c(looic_g, looic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(waic_diff_sig = map_lgl(loo_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))

## plots
plot_data <- roc8 %>% 
  filter(exp != dataset8[1]) %>% 
  group_by(exp) %>% 
  summarise(across(c(OLD_4new:NEW_4old), sum)) %>% 
  pivot_longer(-exp, names_to = c("status", "response"), names_sep = "_") %>% 
  group_by(exp, status) %>% 
  mutate(observed = value / sum(value)) %>% 
   mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("4new", "3new", "2new", "1new", 
                                           "1old", "2old", "3old", "4old")))
  
pred_gumbel <- lapply(roc8_fits_gumbel[-1], posterior_epred)
pred_uvsdt <- lapply(roc8_fits_uvsdt[-1], posterior_epred)

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
  filter(response != "4new") %>% 
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
  filter(response != "4new") %>% 
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

