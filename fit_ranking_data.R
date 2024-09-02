library("tidyverse")
library("brms")
options(mc.cores = parallel::detectCores())
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))

source("data_from_david.R")

source("gumbel6agg-stan.R")
source("uvsdt6agg-stan.R")

fit_dob231_uvsd <- xxx
