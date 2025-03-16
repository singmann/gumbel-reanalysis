
library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 12) + 
            theme(legend.position="bottom"))

load("malejka-broeder.rda")
source("bin-roc-data.R")
source("gumbelbin-stan.R")
source("uvsdtbin-stan.R")

#### new data sets

head(dbin5point) ## miss, hit, cr, false alarm 

colnames(dbin5point) <- outer(c("miss", "hit", "cr", "fa"), c("cr1", "cr2", "cr3", "cr4", "cr5"), 
      FUN = "paste", sep = "_") %>% 
  as.vector()

bdin_all <- dbin5point %>% 
  as_tibble() 
bdin_all$id <- rownames(dbin5point)
bdin_all <- bdin_all %>% 
  separate(id, into = c("exp", "pid"), sep = ", ") %>% 
  mutate(pid = str_trim(pid))

bdin_long <- bdin_all %>%
  pivot_longer(cols = -c(exp, pid), 
               names_to = c("type", "baserate"), names_sep = "_") %>% 
  pivot_wider(names_from = "type") %>% 
  mutate(Nold = hit + miss, 
         Nnew = fa + cr)

mb_extra <- bind_rows(mbe1, mbe3) %>% 
  select(experiment, Subject, BaseRate, everything()) %>% 
  mutate(Subject = as.character(Subject))
colnames(mb_extra) <- colnames(bdin_long)
bdin_long <- bind_rows(bdin_long, mb_extra)
bdin_long <- bdin_long %>% 
  mutate(exp = factor(
    exp, 
    levels = c("Broder E3", 
               "Dube E1 pictures", "Dube E1 word", "Dube E2", 
                "e1", "Malejka_e2", "e3",
               "Van Zandt E1 slow", "Van Zandt E1 fast", "Van Zandt E2"), 
    labels = c("Broder E3", 
               "Dube E1-P", "Dube E1-W", "Dube E2", 
                "Malejka E1", "Malejka E2", "Malejka E3",
               "Van Zandt E1 slow", "Van Zandt E1 fast", "Van Zandt E2")))

dataset_all <- unique(bdin_long$exp)

gumbel_formula_2 <- brmsformula(
  hit | vint(Nold, fa, Nnew) ~ 1 + (1|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = gumbelbin_family, cmc = FALSE
)

get_prior(gumbel_formula_2, data = bdin_long)

gumbel_priors <- prior(normal(0,0.5), class = b, dpar = "cr") + 
  prior(student_t(3, 1, 2), class = Intercept)

uvsdt_formula_2 <- brmsformula(
  hit | vint(Nold, fa, Nnew) ~ 1 + (1|p|pid), 
  discsignal ~ 1 + (1|p|pid), 
  cr ~ 0 + baserate + (0 + baserate|p|pid),
  family = uvsdtbin_family, cmc = FALSE
)

uvsdt_priors <- prior(normal(0,0.5), class = b, dpar = "cr") + 
  prior(student_t(3, 1, 2), class = Intercept) +
  prior(student_t(3, 0.5, 1), class = Intercept, dpar = "discsignal")


rocbin_data <- vector("list", length(dataset_all))

rocbin_fits_gumbel <- vector("list", length(dataset_all))
rocbin_fits_uvsdt <- vector("list", length(dataset_all))

rocbin_exloo_gumbel <- vector("list", length(dataset_all))
rocbin_exloo_uvsdt <- vector("list", length(dataset_all))

library(future)
plan(multisession, workers = 16)

start_time <- Sys.time()

for (i in seq_along(dataset_all)) {
  print(i)
  rocbin_data[[i]] <- bdin_long %>% 
    filter(exp == dataset_all[i])
  rocbin_fits_gumbel[[i]] <- brm(
    gumbel_formula_2, data = rocbin_data[[i]], 
    stanvars = sv_gumbelbin ,
    prior = gumbel_priors,
    init_r = 0.5
  )
  rocbin_exloo_gumbel[[i]] <- kfold(
    x = rocbin_fits_gumbel[[i]], group = "pid", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_gumbelbin",
                                          "calc_posterior_predictions_gumbelbin", 
                                          "posterior_epred_gumbelbin", "posterior_predict_gumbelbin")))
  
  rocbin_fits_uvsdt[[i]] <- brm(
    uvsdt_formula_2, data = rocbin_data[[i]], 
    stanvars = sv_uvsdtbin, 
    prior = uvsdt_priors,
    init_r = 0.5
  )
  rocbin_exloo_uvsdt[[i]] <- kfold(
    x = rocbin_fits_uvsdt[[i]], group = "pid", sample_new_levels = "uncertainty",
    future_args = list(future.globals = c("log_lik_uvsdtbin", "calc_posterior_predictions_uvsdtbin", 
                                          "posterior_epred_uvsdtbin", "posterior_predict_uvsdtbin")))
  
}
end_time <- Sys.time()
# Time difference
time_elapsed <- end_time - start_time
print(time_elapsed)


rocbin_fits_gumbel[[1]]

### make plots

pred_gumbel <- lapply(rocbin_fits_gumbel, posterior_epred)
names(pred_gumbel) <- dataset_all
pred_uvsd <- lapply(rocbin_fits_uvsdt, posterior_epred)
names(pred_uvsd) <- dataset_all

str(pred_gumbel, 1)

plot_dat <- bdin_long %>% 
  group_by(exp, pid, baserate) %>% 
  summarise(
    hit = hit/Nold,
    fa = fa/Nnew,
  ) 

plot_dat %>% 
  ungroup() %>% 
  count(exp) 

plot_dat <- bind_cols(
  plot_dat, 
  map_dfr(pred_gumbel[levels(bdin_long$exp)], ~as.data.frame(apply(., c(2, 3), mean)))
) %>% 
  rename(
    hit_gumbel = old, fa_gumbel = new
  )
plot_dat <- bind_cols(
  plot_dat, 
  map_dfr(pred_uvsd[levels(bdin_long$exp)], ~as.data.frame(apply(., c(2, 3), mean)))
) %>% 
  rename(
    hit_uvsd = old, fa_uvsd = new
  )

plot_dat_2 <- plot_dat %>% 
  group_by(exp, baserate) %>% 
  summarise(across(c(hit,fa, hit_gumbel, fa_gumbel, hit_uvsd, fa_uvsd), mean))

plot_data_bin <- plot_dat_2

save(plot_data_bin, exloo_bin, 
     rocbin_exloo_gumbel, rocbin_exloo_uvsdt, 
     file = "rocbin_exloo_res.rda")

exloo_bin <- mapply(loo_compare, rocbin_exloo_gumbel, rocbin_exloo_uvsdt, SIMPLIFY = FALSE)

bin_comp <- tibble(
  exp = dataset_all,
  elpd_g = map_dbl(rocbin_exloo_gumbel, ~.$estimates["elpd_kfold","Estimate"]),
  elpd_uv = map_dbl(rocbin_exloo_uvsdt, ~.$estimates["elpd_kfold","Estimate"]),
) %>% 
  mutate(max_elpd = pmax(elpd_g, elpd_uv)) %>% 
  mutate(across(c(elpd_g, elpd_uv), ~ sprintf(.-max_elpd, fmt = '%#.1f'))) %>% 
  mutate(diff_SE = map_dbl(exloo_bin, ~ .[2, "se_diff"])) %>% 
  mutate(diff_sig = map_lgl(exloo_bin, ~ abs(.[2, "elpd_diff"]) > (2*.[2, "se_diff"])))
bin_comp
#    dataset           elpd_g elpd_uv max_elpd diff_SE diff_sig
#    <fct>             <chr>  <chr>      <dbl>   <dbl> <lgl>   
#  1 Broder E3         -0.8   0.0        -992.    1.77 FALSE   
#  2 Dube E1-P         -3.0   0.0       -1114.    4.38 FALSE   
#  3 Dube E1-W         -1.1   0.0       -1155.    3.43 FALSE   
#  4 Dube E2           0.0    -0.3       -846.    1.92 FALSE   
#  5 Van Zandt E1 slow 0.0    -11.1      -336.    6.98 FALSE   
#  6 Van Zandt E1 fast -3.7   0.0        -307.    3.26 FALSE   
#  7 Van Zandt E2      -3.8   0.0        -404.    8.12 FALSE   
#  8 Malejka E2        -2.5   0.0        -833.    8.88 FALSE   
#  9 Malejka E1        -6.1   0.0       -1010.    7.55 FALSE   
# 10 Malejka E3        0.0    -2.4       -637.    2.72 FALSE   

bin_n <- plot_dat %>% 
  group_by(exp) %>% 
  summarise(n = n_distinct(pid))

bin_comp <- left_join(bin_comp, bin_n)

bin_comp <- bin_comp %>% 
  mutate(
    n_text = paste0("italic(N) == ", n),
    elpd_diff = paste0(
      "paste(Δ[ELPD] == '", 
      if_else(elpd_g == "0.0", str_replace(elpd_uv, "-", "+") , elpd_g),
      "', ' (±", sprintf(diff_SE, fmt = '%#.1f'), ")')"
    )
  )

library(showtext)
font_paths("~/.fonts/")

# Add font
font_add("Palatino Linotype", 
         #regular="pala.ttf", 
         regular = "asana-math.otf",
         #regular = "palatinolinotype_roman.ttf",
         italic = "palatinolinotype_italic.ttf", 
         bold = "palatinolinotype_bold.ttf", 
         bolditalic = "palatinolinotype_bolditalic.ttf")
showtext_auto()

theme_set(theme_bw(base_size = 12, base_family = "Palatino Linotype") + 
            theme(legend.position="bottom", 
                  panel.grid = element_blank()))

psize <- 3.5
lsize <- 1.5
plot_dat_2 %>%
  ggplot(aes(x =  fa, y = hit)) +
  geom_abline(slope = -1, intercept = 1, linetype = 2) +
  annotate(geom = "polygon", 
           x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "grey")+
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(aes(group = 1), linewidth = lsize) +
  geom_point(size = psize, aes(shape = "Data", colour = "Data")) +
  geom_point(aes(x = fa_gumbel, y = hit_gumbel, 
                 shape = "Gumbel", colour = "Gumbel"), size = psize) +
  geom_point(aes(x = fa_uvsd, y = hit_uvsd, 
                 shape = "UVSD", colour = "UVSD"), size = psize) + 
  geom_label(mapping = aes(x = 0.55, y = 0.45, label = n_text), 
             data = bin_comp, hjust = "left", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  geom_label(mapping = aes(x = 0.2, y = 0.15, label = elpd_diff), 
             data = bin_comp, hjust = "left", vjust = "top", parse = TRUE,
             family = "Palatino Linotype") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
  facet_wrap(vars(exp)) + 
  scale_color_manual(
     name = '',
     breaks = c('Data', 'Gumbel', 'UVSD'),
     values = c('Data' = 'black', 'Gumbel' = 'blue', 'UVSD' = 'red')
   ) +
  scale_shape_manual(
     name = '',
     breaks = c('Data', 'Gumbel', 'UVSD'),
     values = c('Data' = 19, 'Gumbel' = 3, 'UVSD' = 2)
   ) +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.1), 
        legend.direction = "horizontal") +
  labs(x = expression(italic(p)[FA]), y = expression(italic(p)[H]))
ggsave("binroc-plot1.pdf", width = 23, height = 20, units = "cm")
