gumbelrank_stanvars <- "
  real gumbelrank_lpmf(int y, real mu,
                   int y1, int y2, int y3) {
  vector[4] p;
  array[4] int respvec = { y, y1, y2, y3 };
  real expg = exp(mu);
  real expmg = exp(-mu);
  array[3] real fixgamma = { 1, 1, 2 };
  
  for (i in 1:3) {
    p[i] = expg * ( ( 6.0 * tgamma(i + expmg - 1) ) / ( fixgamma[i] * tgamma(4 + expmg) ) ); 
  }

  p[4] = 1-p[1]-p[2]-p[3];

  return multinomial_lpmf(respvec | p);
  }
"

gumbelrank_family <- custom_family(
  name = "gumbelrank", 
  dpars = c("mu"), 
  links = c("identity"), lb = c(NA),
  type = "int", vars = c(paste0("vint", 1:3, "[n]"))
)
sv_gumbelrank <- stanvar(scode = gumbelrank_stanvars, block = "functions") 

calc_posterior_predictions_gumbelrank <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- brms::get_dpar(prep, "discsignal", i = i)

  OUTLEN <- length(mu)
  
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = 4)
  
  G1<-function(x, mu, sd){
    ((pnorm(x)^3)*dnorm(x,mean=mu,sd=ss))
  }
  
  G2<-function(x, mu, sd){
    ((pnorm(x)^2)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x)))*3
  }
  
  G3<-function(x, mu, sd){
    (pnorm(x)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x))^2)*3
  }

  for (j in seq_len(OUTLEN)) {
    p[j, 1] <- integrate(G1,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value    
    p[j, 2] <- integrate(G2,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value
    p[j, 3] <- integrate(G3,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value
  }
  p[, 4] <- 1-p[,1]-p[,2]-p[,3]
  return(p)
}

log_lik_gumbelrank <- function(i, prep) {
  p <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
  dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i])
  extraDistr::dmnom(x = dvec, size = sum(dvec), prob = p, log = TRUE)
}

posterior_epred_gumbelrank <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 4), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("R1", "R2", "R3", "R4")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}

posterior_predict_gumbelrank <- function(i, prep, ...) {
  p <- calc_posterior_predictions_gumbelrank(i = i, prep = prep)
  dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i])
  
  lout <- length(p)
  out <- extraDistr::rmnom(n = rep(1, lout), size = sum(dvec), prob = p)  
  colnames(out) <- c("R1", "R2", "R3", "R4")
  #browser()
  lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}

