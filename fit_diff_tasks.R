library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

source("non-ranking-data-from-david.R")
source("gumbel6agg-stan.R")
source("uvsdt6agg-stan.R")
source("6point_other_predictions.R") ## this override theprediction functions with new ones for the 2AFC task
library("tidybayes")

(colSums(dobbins23_fc1[,-3]) / sum(dobbins23_fc1[,-3])) %>% 
  print(digits = 3)
  # correct incorrect 
  #   0.742     0.258 

(colSums(dobbins23_fc2[,-3]) / sum(dobbins23_fc2[,-3])) %>% 
  print(digits = 3)
  # correct incorrect 
  #   0.742     0.258 
  

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


##----------------------------------------------------------------
##                            Exp. 1                             -
##----------------------------------------------------------------



fit_dob23_1_gumbel <- brm(
    gumbel_formula, data = dobbins23_yn1, 
    stanvars = sv_gumbel6agg ,
    prior = gumbel_priors,
    init_r = 0.5
  )

pred_dob1_gumbel <- posterior_epred(fit_dob23_1_gumbel)
gp2 <- apply(pred_dob1_gumbel, c(1,3), mean)
apply(gp2, 2, mean) %>% 
  print(digits = 2)
  # correct incorrect 
  #    0.73      0.27 
apply(gp2, 2, quantile, probs = c(0.027, 0.975)) %>% 
  print(digits = 2)
#       correct incorrect
# 2.7%     0.72      0.27
# 97.5%    0.73      0.28

str(pred_dob1_gumbel)

gp2b <- apply(pred_dob1_gumbel, c(2,3), mean) %>% 
  as.data.frame()
colnames(gp2b) <- c("gumbel", "incorrect_gumbel")

fit_dob23_1_uvsd <-  brm(
    uvsdt_formula, data = dobbins23_yn1, 
    stanvars = sv_uvsdt6agg, 
    prior = uvsdt_priors,
    init_r = 0.5
  )

pred_dob1_uvsd <- posterior_epred(fit_dob23_1_uvsd)
up2 <- apply(pred_dob1_uvsd, c(1,3), mean)
apply(up2, 2, mean) %>% 
  print(digits = 2)
  # correct incorrect 
  #    0.74      0.26 
apply(up2, 2, quantile, probs = c(0.027, 0.975)) %>% 
  print(digits = 2)
#       correct incorrect
# 2.7%     0.73      0.25
# 97.5%    0.75      0.27

up2b <- apply(pred_dob1_uvsd, c(2,3), mean) %>% 
  as.data.frame()
colnames(up2b) <- c("uvsd", "incorrect_uvsd")

head(dobbins23_fc1)

dobb1_pred <- dobbins23_fc1 %>% 
  mutate(data = correct / (correct + incorrect)) %>% 
  bind_cols(gp2b, up2b) %>% 
  as_tibble()

cor.test(~ data + gumbel, dobb1_pred, method = "spearman")
cor.test(~ data + uvsd, dobb1_pred, method = "spearman")

p1 <- ggplot(dobb1_pred, aes(x = data, y = gumbel)) +
  geom_point() +
  geom_smooth()

p2 <- ggplot(dobb1_pred, aes(x = data, y = uvsd)) +
  geom_point() +
  geom_smooth()

cowplot::plot_grid(p1, p2)

##----------------------------------------------------------------
##                            Exp. 2                             -
##----------------------------------------------------------------

fit_dob23_2_gumbel <- brm(
    gumbel_formula, data = dobbins23_yn2, 
    stanvars = sv_gumbel6agg ,
    prior = gumbel_priors,
    init_r = 0.5, 
    control = list(adapt_delta = 0.99)
  )

fit_dob23_2_uvsd <-  brm(
    uvsdt_formula, data = dobbins23_yn2, 
    stanvars = sv_uvsdt6agg, 
    prior = uvsdt_priors,
    init_r = 0.5
  )

dob2pred <- left_join(dobbins23_yn2, dobbins23_fc2)

dobbins23_fc2_b <- dob2pred %>% 
  add_epred_draws(fit_dob23_2_uvsd)

dobbins23_fc2_b <- dobbins23_fc2_b %>% 
  mutate(data = correct / (correct + incorrect)) %>% 
  ungroup() %>% 
  select(-c(OLD_3new:NEW_3old))
  
dobbins23_fc2_b %>% 
  filter(.category == "correct") %>% 
  ggplot(aes(x = data, y = .epred)) +
  geom_point(alpha = 0.01)

dobb2_p_uvsd <- dobbins23_fc2_b %>% 
  group_by(id, .category) %>% 
  summarise(
    acc = correct[1] / (correct[1] + incorrect[1]),
    uvsd_m = mean(.epred),
    uvsd_l = quantile(.epred, prob = 0.025),
    uvsd_h = quantile(.epred, prob = 0.975)
  ) %>% 
  filter(.category == "correct") 

dobb2_p_uvsd %>% 
  #mutate(acc = if_else(acc < .5,  1- acc, acc)) %>% 
  arrange(uvsd_l) %>% 
  ggplot(aes(y = factor(id, levels = unique(id)))) +
  geom_linerange(aes(xmin = uvsd_l, xmax = uvsd_h)) +
  geom_point(aes(x = acc))
