library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

source("data_from_david.R")

source("gumbelrank2-stan.R")
source("uvsdtrank-stan.R")

gumbel_priors <- prior(student_t(3, 1, 2), class = Intercept)
uvsdt_priors <- prior(student_t(3, 0.5, 1), 
                      class = Intercept, dpar = "discsignal") +
  prior(student_t(3, 1, 2), class = Intercept)

##---------------------------------------------------------------
##                            4 Ranks                           -
##---------------------------------------------------------------

##------------
##  KKS 2012  
##------------

gumbel_formula_kks <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

# stancode(gumbel_formula_kks, data = kks12, 
#          stanvars = sv_gumbelrank, prior = gumbel_priors)

fit_kks_gumbel <- brm(
    gumbel_formula_kks, data = kks12, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kks <- brmsformula(
  V1 | vint(V2, V3, V4) ~ 1 + (1|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
)

# stancode(uvsdt_formula_kks, data = kks12, 
#          stanvars = sv_uvsdtrank, prior = uvsdt_priors)

fit_kks_uvsdt <- brm(
    uvsdt_formula_kks, data = kks12, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.5
  )

kks12


##-----------
##  KK14 E1  
##-----------

kk14_e1_use <- kk14_e1 %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

gumbel_formula_kke1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

fit_kke1_gumbel <- brm(
    gumbel_formula_kke1, data = kk14_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kke1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
)

fit_kke1_uvsdt <- brm(
    uvsdt_formula_kke1, data = kk14_e1_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##------------
##  MHE22 E1  
##------------


mhe_e1_use <- mhe_e1 %>% 
  pivot_longer(cols = rank1.w:rank4.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

gumbel_formula_mhe1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  family = gumbelrank_family, cmc = FALSE
)

fit_mhe1_gumbel <- brm(
    gumbel_formula_mhe1, data = mhe_e1_use, 
    stanvars = sv_gumbelrank, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_mhe1 <- brmsformula(
  rank1 | vint(rank2, rank3, rank4) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank_family, cmc = FALSE
)

fit_mhe1_uvsdt <- brm(
    uvsdt_formula_mhe1, data = mhe_e1_use, 
    stanvars = sv_uvsdtrank, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

save(fit_kke1_gumbel, fit_kke1_uvsdt, fit_kks_gumbel, fit_kks_uvsdt, 
     fit_mhe1_gumbel, fit_mhe1_uvsdt, file = "fit-4rank.rda", compress = "xz")

##---------------------------------------------------------------
##                            3 Ranks                           -
##---------------------------------------------------------------
source("gumbelrank2-3r-stan.R")
source("uvsdtrank3-stan.R")

##-----------
##  KK14 E2  
##-----------

kk14_e2

kk14_e2_use <- kk14_e2 %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))


gumbel_formula_kke2 <- brmsformula(
  rank1 | vint(rank2, rank3) ~ strength + (strength|p|id), 
  family = gumbelrank3_family, cmc = FALSE
)

fit_kke2_gumbel <- brm(
    gumbel_formula_kke2, data = kk14_e2_use, 
    stanvars = sv_gumbelrank3, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

uvsdt_formula_kke2 <- brmsformula(
  rank1 | vint(rank2, rank3) ~ strength + (strength|p|id), 
  discsignal ~ 1 + (1|p|id),
  family = uvsdtrank3_family, cmc = FALSE
)

fit_kke2_uvsdt <- brm(
    uvsdt_formula_kke2, data = kk14_e2_use, 
    stanvars = sv_uvsdtrank3, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##-----------
##  MG16 E1  
##-----------
#McAdoo and Gronlund (2016)

mg16_e1
mg16_e1_use <- mg16_e1 %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

fit_mge1_gumbel <- brm(
    gumbel_formula_kke2, data = mg16_e1_use, 
    stanvars = sv_gumbelrank3, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )
fit_mge1_uvsdt <- brm(
    uvsdt_formula_kke2, data = mg16_e1_use, 
    stanvars = sv_uvsdtrank3, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )


##-----------
##  MG16 E2  
##-----------

mg16_e2
mg16_e2_use <- mg16_e2 %>% 
  pivot_longer(cols = rank1.w:rank3.s, 
               names_to = c("rank", "strength"), names_sep = "\\.") %>% 
  pivot_wider(names_from = rank, values_from = value) %>% 
  mutate(strength = factor(strength, levels = c("w", "s")))

fit_mge2_gumbel <- brm(
    gumbel_formula_kke2, data = mg16_e2_use, 
    stanvars = sv_gumbelrank3, 
    prior = gumbel_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )
fit_mge2_uvsdt <- brm(
    uvsdt_formula_kke2, data = mg16_e2_use, 
    stanvars = sv_uvsdtrank3, 
    prior = uvsdt_priors,
    init_r = 0.25
    #control = list(adapt_delta = 0.99)
  )

save(fit_kke2_gumbel, fit_kke2_uvsdt, 
     fit_mge1_gumbel, fit_mge1_uvsdt,
     fit_mge2_gumbel, fit_mge2_uvsdt, file = "fit-3rank.rda", compress = "xz")
load("fit-3rank.rda")



##----------------------------------------------------------------
##                              Plot                             -
##----------------------------------------------------------------

kks12_agg <- kks12 %>% 
  mutate(exp = "Kellen (2012)") %>% 
  group_by(exp) %>% 
  summarise(across(V1:V4, sum)) %>% 
  mutate(total = V1 + V2 + V3 + V4) %>% 
  mutate(across(V1:V4, ~./total)) %>% 
  select(-total) %>% 
  pivot_longer(cols = -exp, names_to = "rank", values_to = "observed") %>% 
  mutate(rank = str_extract(rank, "\\d"))

kks12_gumbel <- posterior_epred(fit_kks_gumbel)
kks12_agg$gumbel <- apply(kks12_gumbel, 3, mean)
kks12_uvsdt <- posterior_epred(fit_kks_uvsdt)
kks12_agg$uvsdt <- apply(kks12_uvsdt, 3, mean)

d_4r_strength <- bind_rows(
  mutate(kk14_e1_use, exp = "Klauer (2014, E1)"),
  mutate(mhe_e1_use)
)


