
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
  family = gumbel8agg_family, cmc = FALSE
)

gumbel_priors_8 <- prior(normal(0,0.5), class = Intercept, dpar = "crc") + 
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crll") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crlx") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhm") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhh") +
  prior(normal(-0.5,0.5), class = Intercept, dpar = "crhx") +
  prior(student_t(3, 1, 2), class = Intercept)

uvsdt_formula_8 <- brmsformula(
  OLD_4new | vint(OLD_3new, OLD_2new, OLD_1new, OLD_1old, OLD_2old, OLD_3old, OLD_4old, NEW_4new, NEW_3new, NEW_2new, NEW_1new, NEW_1old, NEW_2old, NEW_3old, NEW_4old) ~ 1 + (1|p|id), 
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

roc8_fits_gumbel <- vector("list", length(dataset8))
roc8_fits_uvsdt <- vector("list", length(dataset8))

roc8_exloo_gumbel <- vector("list", length(dataset8))
roc8_exloo_uvsdt <- vector("list", length(dataset8))

library(future)
plan(multisession, workers = 16)

start_time <- Sys.time()
#i <- 1
for (i in seq_along(dataset8)) {
  print(i)
  roc8_data[[i]] <- roc8_use %>% 
    filter(exp == dataset8[i])
  
  roc8_fits_gumbel[[i]] <- brm(
    gumbel_formula_8, data = roc8_data[[i]], 
    stanvars = sv_gumbel8agg, 
    prior = gumbel_priors_8,
    init_r = 0.25, control = list(adapt_delta = 0.9999)
  )
  roc8_exloo_gumbel[[i]] <- kfold(
    x = roc8_fits_gumbel[[i]], group = "id", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_gumbel8agg", "calc_posterior_predictions_gumbel8agg", 
                                          "posterior_epred_gumbel8agg", "posterior_predict_gumbel8agg")))

  roc8_fits_uvsdt[[i]] <- brm(
    uvsdt_formula_8, data = roc8_data[[i]], 
    stanvars = sv_uvsdt8agg, 
    prior = uvsdt_priors_8,
    init_r = 0.5, control = list(adapt_delta = 0.9999)
  )
  roc8_exloo_uvsdt[[i]] <- kfold(
    x = roc8_fits_uvsdt[[i]], group = "id", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_uvsdt8agg", "calc_posterior_predictions_uvsdt8agg", 
                                          "posterior_epred_uvsdt8agg", "posterior_predict_uvsdt8agg")))

}
end_time <- Sys.time()
# Time difference
time_elapsed <- end_time - start_time
print(time_elapsed) ## 5 days

save(roc8_exloo_gumbel, roc8_exloo_uvsdt, file = "roc8-exloo.rda")

exloo_8roc <- mapply(loo_compare, roc8_exloo_gumbel, roc8_exloo_uvsdt, SIMPLIFY = FALSE)
tibble(
  dataset = dataset8,
  elpd_g = map_dbl(roc8_exloo_gumbel, ~.$estimates["elpd_kfold","Estimate"]),
  elpd_uv = map_dbl(roc8_exloo_uvsdt, ~.$estimates["elpd_kfold","Estimate"]),
) %>% 
  mutate(max_elpd = pmax(elpd_g, elpd_uv)) %>% 
  mutate(across(c(elpd_g, elpd_uv), ~ sprintf(.-max_elpd, fmt = '%#.1f'))) %>% 
  mutate(diff_SE = map_dbl(exloo_8roc, ~ .[2, "se_diff"])) %>% 
  mutate(diff_sig = map_lgl(exloo_8roc, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 3 × 6
#   dataset           elpd_g elpd_uv max_elpd diff_SE diff_sig
#   <chr>             <chr>  <chr>      <dbl>   <dbl> <lgl>   
# 1 Benjamin_2013     0.0    -3.0      -3778.    11.4 FALSE   
# 2 Onyper_2010-Pics  0.0    -77.0    -10221.    61.4 FALSE   
# 3 Onyper_2010-Words 0.0    -61.4    -10351.    56.0 FALSE  


# stancode(gumbel_formula, data = roc6_data[[i]], 
#          stanvars = roc6_sv[[i]], 
#          prior = gumbel_priors)



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
  dataset = dataset8,
  waic_g = map_dbl(roc8_waic_gumbel, ~.$estimates["waic","Estimate"]),
  waic_uv = map_dbl(roc8_waic_uvsdt, ~.$estimates["waic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(waic_g, waic_uv)) %>% 
  mutate(across(c(waic_g, waic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(waic_SE = map_dbl(waic_comp, ~ .[2, "se_diff"])) %>% 
  mutate(waic_diff_sig = map_lgl(waic_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 3 × 6
#   dataset           waic_g waic_uv min_waic waic_SE waic_diff_sig
#   <chr>             <chr>  <chr>      <dbl>   <dbl> <lgl>        
# 1 Benjamin_2013     "  0"  " 2.3"     6251.    8.05 FALSE        
# 2 Onyper_2010-Pics  "  0"  "72.9"    12418.   16.9  TRUE         
# 3 Onyper_2010-Words "107"  " 0.0"    11885.   31.3  FALSE    

tibble(
  dataset = dataset8,
  looic_g = map_dbl(roc8_loo_gumbel, ~.$estimates["looic","Estimate"]),
  looic_uv = map_dbl(roc8_loo_uvsdt, ~.$estimates["looic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(looic_g, looic_uv)) %>% 
  mutate(across(c(looic_g, looic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(looic_SE = map_dbl(loo_comp, ~ .[2, "se_diff"])) %>% 
  mutate(looic_diff_sig = map_lgl(loo_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 3 × 6
#   dataset           looic_g looic_uv min_waic looic_SE looic_diff_sig
#   <chr>             <chr>   <chr>       <dbl>    <dbl> <lgl>         
# 1 Benjamin_2013     " 5.6"  " 0"        6441.     10.4 FALSE         
# 2 Onyper_2010-Pics  " 0.0"  "53"       12724.     19.2 FALSE         
# 3 Onyper_2010-Words "87.5"  " 0"       12174.     32.2 FALSE


## plots
plot_data <- roc8 %>% 
  #filter(exp != dataset8[1]) %>% 
  group_by(exp) %>% 
  summarise(across(c(OLD_4new:NEW_4old), sum)) %>% 
  pivot_longer(-exp, names_to = c("status", "response"), names_sep = "_") %>% 
  group_by(exp, status) %>% 
  mutate(observed = value / sum(value)) %>% 
   mutate(
    status = factor(status, levels = c("OLD", "NEW")), 
    response = factor(response, levels = c("4new", "3new", "2new", "1new", 
                                           "1old", "2old", "3old", "4old")))
  
pred_gumbel <- lapply(roc8_fits_gumbel, posterior_epred)
pred_uvsdt <- lapply(roc8_fits_uvsdt, posterior_epred)

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

plot_data_roc8 <- plot_data
plot_data2_roc8 <- plot_data2

save(plot_data_roc8, plot_data2_roc8, 
     roc8_exloo_gumbel, roc8_exloo_uvsdt, exloo_8roc, 
     file = "roc8_exloo_res.rda")

# str(apply(pred_gumbel[[1]], c(1, 3), mean))
# str(pred_gumbel)

psize <- 3.5
lsize <- 1.5

plot_data %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "4new") %>% 
  ggplot(aes(x =  observed_NEW, y = observed_OLD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize) +
  geom_point(aes(x = gumbel_NEW, y = gumbel_OLD), shape = 3, colour = "blue", size = psize) +
  geom_point(aes(x = uvsd_NEW, y = uvsd_OLD), shape = 2, colour = "red", size = psize) + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(exp)) +
  labs(x = "False alarms", y = "Hits") +
  facet_wrap(vars(exp))
ggsave("roc8-plot1.pdf", width = 16, height = 8, units = "cm")


psize <- 2.5
lsize <- 1.25
ebsize <- 1.15
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
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize) +
    geom_linerange(aes(x = uvsd_NEW, y = uvsd_OLD, ymin = uvsd_low_OLD, ymax = uvsd_high_OLD), 
                   colour = "red", linewidth = ebsize) +
   geom_linerange(aes(x = uvsd_NEW, y = uvsd_OLD, xmin = uvsd_low_NEW, xmax = uvsd_high_NEW), 
                  colour = "red", linewidth = ebsize) +
    geom_linerange(aes(x = gumbel_NEW, y = gumbel_OLD, ymin = gumbel_low_OLD, ymax = gumbel_high_OLD), 
                 colour = "blue", linewidth = ebsize) +
   geom_linerange(aes(x = gumbel_NEW, y = gumbel_OLD, xmin = gumbel_low_NEW, xmax = gumbel_high_NEW), 
                  colour = "blue", linewidth = ebsize) +
  geom_point(aes(x = uvsd_NEW, y = uvsd_OLD), shape = 2, colour = "red", size = psize) +
  geom_point(aes(x = gumbel_NEW, y = gumbel_OLD), shape = 4, colour = "blue", size = psize) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(vars(exp)) +
  labs(x = "False alarms", y = "Hits")
ggsave("roc8-plot2.pdf", width = 16, height = 8, units = "cm")
