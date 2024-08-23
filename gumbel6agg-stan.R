gumbel6agg_stanvars <- "
   real gumbelmin(real x, real mu, real disc){
     //return 1- exp(-exp(-(-x-mu)/disc));
     //return exp(gumbel_lccdf(-x | mu,disc));
     return 1 - gumbel_cdf(-x|mu,disc);
   }
   real gumbel6p_lpmf(int y, real mu, 
                   real crc, real crlm, real crll, real crhm, real crhh, 
                   int itemtype; array dvec) {
     int nthres = 5;
     vector[nthres] p;
     real disc = 1;
     vector[nthres] thres;
     
     // calculate thresholds
     thres[1] = crc - (exp(crlm) + exp(crll));
     thres[2] = crc - (exp(crlm));
     thres[3] = crc;
     thres[4] = crc + (exp(crhm));
     thres[5] = crc + (exp(crhm) + exp(crhh));
     
     // calculate probabilities
     if (itemtype == 1) { // old items
      p[1] = gumbelmin(thres[1], mu, disc);
      for (i in 2:nthres) {
        p[i] = gumbelmin(thres[i], mu, disc) - gumbelmin(thres[i-1], mu, disc);
      }
     } else { // new items
      p[1] = gumbelmin(thres[1], 0, disc);
      for (i in 2:nthres) {
        p[i] = gumbelmin(thres[i], 0, disc) - gumbelmin(thres[i-1], 0, disc);
      }
     }

     return multinomial_lpmf(dvec | p);
   }
"

gumbel6agg_family <- custom_family(
  name = "gumbel6agg", 
  dpars = c("mu", "crc", "crlm", "crll", "crhm", "crhh"), 
  links = c("identity", rep("identity", 5)), lb = c(NA, rep(NA, 5)),
  type = "int", vars = c("vint1[n]", "dmat[n]")
)
sv_gumbel6agg <- stanvar(scode = gumbel6agg_stanvars, block = "functions")

# log_lik_gumbel6p <- function(i, prep) {
#   mu <- brms::get_dpar(prep, "mu", i = i)
#   #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
#   crc <- brms::get_dpar(prep, "crc", i = i)
#   crlm <- brms::get_dpar(prep, "crlm", i = i)
#   crll <- brms::get_dpar(prep, "crll", i = i)
#   crhm <- brms::get_dpar(prep, "crhm", i = i)
#   crhh <- brms::get_dpar(prep, "crhh", i = i)
#   
#   itemtype <- prep$data$vint1[i]
#   y <- prep$data$Y[i]
#   gumbel6p_lpmf(y, mu, crc, crlm, crll, crhm, crhh, itemtype)
# }
# posterior_predict_gumbel6p <- function(i, prep, ...) {
#   mu <- brms::get_dpar(prep, "mu", i = i)
#   #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
#   crc <- brms::get_dpar(prep, "crc", i = i)
#   crlm <- brms::get_dpar(prep, "crlm", i = i)
#   crll <- brms::get_dpar(prep, "crll", i = i)
#   crhm <- brms::get_dpar(prep, "crhm", i = i)
#   crhh <- brms::get_dpar(prep, "crhh", i = i)
#   itemtype <- prep$data$vint1[i]
#   nsamples <- length(mu)
#   nthres <- 5
#   thres <- matrix(NA_real_, nrow = nsamples, ncol = nthres)
#   p <- matrix(NA_real_, nrow = nsamples, ncol = nthres+1) 
#   disc = rep(1, nsamples)
#   thres[,1] = crc - (exp(crlm) + exp(crll));
#   thres[,2] = crc - (exp(crlm));
#   thres[,3] = crc;
#   thres[,4] = crc + (exp(crhm));
#   thres[,5] = crc + (exp(crhm) + exp(crhh));
#   for (y in 1:(nthres+1)) {
#     if (y == 1) {
#       p[,y] = ordinal::pgumbel(thres[,1], mu, disc, max = FALSE)
#     } else if (y == nthres + 1) {
#       p[,y] = 1 - ordinal::pgumbel(thres[,nthres], mu, disc, max = FALSE)
#     } else {
#       p[,y] = ordinal::pgumbel(thres[,y], mu, disc, max = FALSE) -
#           ordinal::pgumbel(thres[,y-1], mu, disc, max = FALSE)
#     }
#   }
#   apply(p, 1, extraDistr::rcat, n = 1)
# }