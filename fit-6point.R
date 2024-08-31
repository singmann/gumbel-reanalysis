
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

roc6_exloo_gumbel <- vector("list", length(dataset6))
roc6_exloo_uvsdt <- vector("list", length(dataset6))

library(future)
plan(multisession, workers = 16)

start_time <- Sys.time()

for (i in seq_along(dataset6)) {
  print(i)
  roc6_data[[i]] <- roc6_use %>% 
    filter(exp == dataset6[i])
  roc6_fits_gumbel[[i]] <- brm(
    gumbel_formula, data = roc6_data[[i]], 
    stanvars = sv_gumbel6agg ,
    prior = gumbel_priors,
    init_r = 0.5
  )
  roc6_exloo_gumbel[[i]] <- kfold(
    x = roc6_fits_gumbel[[i]], group = "id", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_gumbel6agg", "calc_posterior_predictions_gumbel6agg", 
                                          "posterior_epred_gumbel6agg", "posterior_predict_gumbel6agg")))
  
  roc6_fits_uvsdt[[i]] <- brm(
    uvsdt_formula, data = roc6_data[[i]], 
    stanvars = sv_uvsdt6agg, 
    prior = uvsdt_priors,
    init_r = 0.5
  )
  roc6_exloo_uvsdt[[i]] <- kfold(
    x = roc6_fits_uvsdt[[i]], group = "id", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_uvsdt6agg", "calc_posterior_predictions_uvsdt6agg", 
                                          "posterior_epred_uvsdt6agg", "posterior_predict_uvsdt6agg")))
  
}
end_time <- Sys.time()
# Time difference
time_elapsed <- end_time - start_time
print(time_elapsed)

exloo_6roc <- mapply(loo_compare, roc6_exloo_gumbel, roc6_exloo_uvsdt, SIMPLIFY = FALSE)

str(roc6_exloo_gumbel[[1]])

roc6_exloo_gumbel[[1]]$estimates
dimnames(exloo_6roc[[1]])

as.data.frame(exloo_6roc[[1]])

tibble(
  dataset = dataset6,
  elpd_g = map_dbl(roc6_exloo_gumbel, ~.$estimates["elpd_kfold","Estimate"]),
  elpd_uv = map_dbl(roc6_exloo_uvsdt, ~.$estimates["elpd_kfold","Estimate"]),
) %>% 
  mutate(max_elpd = pmax(elpd_g, elpd_uv)) %>% 
  mutate(across(c(elpd_g, elpd_uv), ~ sprintf(.-max_elpd, fmt = '%#.1f'))) %>% 
  mutate(diff_SE = map_dbl(exloo_6roc, ~ .[2, "se_diff"])) %>% 
  mutate(diff_sig = map_lgl(exloo_6roc, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 12 × 6
#    dataset             elpd_g elpd_uv max_elpd diff_SE diff_sig
#    <chr>               <chr>  <chr>      <dbl>   <dbl> <lgl>   
#  1 Dube_2012-P         -7.8   0.0       -1032.   13.8  FALSE   
#  2 Dube_2012-W         0.0    -1.4      -1155.   19.4  FALSE   
#  3 Heathcote_2006_e1   0.0    -11.2      -769.   13.7  FALSE   
#  4 Heathcote_2006_e2   -32.5  0.0       -1018.   21.1  FALSE   
#  5 Jaeger_2012         -35.4  0.0       -1639.   15.7  TRUE    
#  6 Jang_2009           0.0    -3.9      -1030.    7.12 FALSE   
#  7 Koen_2010_pure      0.0    -0.9      -1528.   17.6  FALSE   
#  8 Koen_2011           -3.6   0.0       -1270.   20.2  FALSE   
#  9 Koen-2013_full      -16.7  0.0       -1507.   10.7  FALSE   
# 10 Koen-2013_immediate 0.0    -3.8      -1769.   16.1  FALSE   
# 11 Pratte_2010         0.0    -9.7      -4403.   30.8  FALSE   
# 12 Smith_2004          0.0    -1.3       -848.    6.96 FALSE   


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
  mutate(waic_SE = map_dbl(waic_comp, ~ .[2, "se_diff"])) %>% 
  mutate(waic_diff_sig = map_lgl(waic_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 12 × 6
#    dataset             waic_g waic_uv min_waic waic_SE waic_diff_sig
#    <chr>               <chr>  <chr>      <dbl>   <dbl> <lgl>        
#  1 Dube_2012-P         " 4.9" " 0.0"     1237.    6.53 FALSE        
#  2 Dube_2012-W         " 0.0" "12.0"     1263.    5.01 FALSE        
#  3 Heathcote_2006_e1   " 0.0" " 1.8"     1051.    7.15 FALSE        
#  4 Heathcote_2006_e2   "86.9" " 0.0"     1438.   16.0  TRUE         
#  5 Jaeger_2012         "85.5" " 0.0"     2645.   11.5  TRUE         
#  6 Jang_2009           " 0.0" "11.2"     1600.    4.28 FALSE        
#  7 Koen_2010_pure      " 0.0" " 3.5"     1845.    8.07 FALSE        
#  8 Koen_2011           "53.6" " 0.0"     1270.   11.9  TRUE         
#  9 Koen-2013_full      " 9.5" " 0.0"     2465.    7.05 FALSE        
# 10 Koen-2013_immediate "16.2" " 0.0"     2690.    7.51 FALSE        
# 11 Pratte_2010         "58.9" " 0.0"     6066.   16.2  FALSE        
# 12 Smith_2004          "14.9" " 0.0"     1446.    5.22 FALSE  

tibble(
  dataset = dataset6,
  looic_g = map_dbl(roc6_loo_gumbel, ~.$estimates["looic","Estimate"]),
  looic_uv = map_dbl(roc6_loo_uvsdt, ~.$estimates["looic","Estimate"]),
) %>% 
  mutate(min_waic = pmin(looic_g, looic_uv)) %>% 
  mutate(across(c(looic_g, looic_uv), ~ format(.-min_waic, digits = 2, format = "g"))) %>% 
  mutate(looic_SE = map_dbl(loo_comp, ~ .[2, "se_diff"])) %>% 
  mutate(looic_diff_sig = map_lgl(loo_comp, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
# # A tibble: 12 × 6
#    dataset             looic_g looic_uv min_waic looic_SE looic_diff_sig
#    <chr>               <chr>   <chr>       <dbl>    <dbl> <lgl>         
#  1 Dube_2012-P         " 7.2"  " 0.0"      1288.     6.92 FALSE         
#  2 Dube_2012-W         " 0.0"  " 6.9"      1310.     4.90 FALSE         
#  3 Heathcote_2006_e1   " 2.6"  " 0.0"      1082.     8.06 FALSE         
#  4 Heathcote_2006_e2   "80.8"  " 0.0"      1485.    16.0  TRUE          
#  5 Jaeger_2012         "90.5"  " 0.0"      2752.    12.3  TRUE          
#  6 Jang_2009           " 0.0"  "18.7"      1654.     6.10 FALSE         
#  7 Koen_2010_pure      " 0.0"  " 1.5"      1913.     9.01 FALSE         
#  8 Koen_2011           "53.4"  " 0.0"      1315.    12.1  TRUE          
#  9 Koen-2013_full      "13.7"  " 0.0"      2545.     7.44 FALSE         
# 10 Koen-2013_immediate "15.8"  " 0.0"      2773.     7.47 FALSE         
# 11 Pratte_2010         "74.7"  " 0.0"      6251.    16.6  TRUE          
# 12 Smith_2004          "22.1"  " 0.0"      1489.     5.48 TRUE    

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

plot_data_roc6 <- plot_data
plot_data2_roc6 <- plot_data2

save(plot_data_roc6, plot_data2_roc6, 
     roc6_exloo_gumbel, roc6_exloo_uvsdt, exloo_6roc, 
     file = "roc6_exloo_res.rda")

# str(apply(pred_gplot_data# str(apply(pred_gumbel[[1]], c(1, 3), mean))
# str(pred_gumbel)

psize <- 3.5
lsize <- 1.5

plot_data %>% 
  select(-value) %>% 
  pivot_wider(names_from = status, values_from = c(observed, gumbel, uvsd)) %>% 
  arrange(exp, desc(response)) %>% 
  mutate(across(-c(response), cumsum)) %>% 
  filter(response != "3new") %>% 
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
  labs(x = "False alarms", y = "Hits")
ggsave("roc6-plot1.pdf", width = 18, height = 17, units = "cm")

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
  filter(response != "3new") %>% 
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
  facet_wrap(vars(exp))
ggsave("roc6-plot2.pdf", width = 18, height = 17, units = "cm")

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
