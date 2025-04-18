gumbelrank_stanvars <- "
  real gumbelrank_lpmf(int y, real mu, int r, int maxrank) {
  real log_p;
  real g = mu;
  real m = maxrank;
  
  if (y == 0) {
    return 0;
  } else {
    real e_neg_g = exp(-g);
    log_p = (-g + lgamma(m) + lgamma(r - 1 + e_neg_g) - lgamma(r) - lgamma(m + e_neg_g));
    return y * log_p;
    ///p = (exp(-g)*tgamma(m)*tgamma(r-1 + exp(-g))) / (tgamma(r)*tgamma(m + exp(-g)));
    ///return y * log(p);
  }
  }
"

gumbelrank_family <- custom_family(
  name = "gumbelrank", 
  dpars = c("mu"), 
  links = c("identity"), lb = c(NA),
  type = "int", vars = c(paste0("vint", 1:2, "[n]"))
)
sv_gumbelrank <- stanvar(scode = gumbelrank_stanvars, block = "functions") 

calc_posterior_predictions_gumbelrank <- function(i, prep) {
  Rr_gumb = function(g, m, r){
    (exp(-g)*gamma(m)*gamma(r-1 + exp(-g))) / (gamma(r)*gamma(m + exp(-g)))
  }
  
  mu <- brms::get_dpar(prep, "mu", i = i)
  OUTLEN <- length(mu)
  
  r <- prep$data$vint1[i]
  m <- prep$data$vint2[i]
  
  p <- Rr_gumb(mu, m, r)
  return(p)
}

### Needs rework, see UVSDT file
# log_lik_gumbelrank <- function(i, prep) {
#   p <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
#   dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i])
#   extraDistr::dmnom(x = dvec, size = sum(dvec), prob = p, log = TRUE)
# }

posterior_epred_gumbelrank <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, prep$ndraws, prep$nobs,
               dimnames = list(seq(prep$ndraws), seq(prep$nobs)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

## likely not possible to fix given data structure
# posterior_predict_gumbelrank <- function(i, prep, ...) {
#   p <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
#   dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i])
#   
#   lout <- length(p)
#   out <- extraDistr::rmnom(n = rep(1, lout), size = sum(dvec), prob = p)  
#   colnames(out) <- c("R1", "R2", "R3", "R4")
#   #browser()
#   lapply(seq_len(nrow(out)), function(i) out[i,])
#   #apply(out, 1, function(x) list(x))
#   #out[,1]
# }
# 
