gumbelmin_dist <- "
  real gumbelmin_cdf(real x, real mu){
    return 1 - gumbel_cdf(-x | -mu, 1);
  }
  real gumbelmin_pdf(real x, real mu){
    return exp(gumbel_lpdf(-x | -mu, 1));
  }
  real gumbelmin_lcdf(real x, real mu){
    return log1m_exp(gumbel_lcdf(-x | -mu, 1));
  }
  real gumbelmin_lccdf(real x, real mu){
    return gumbel_lcdf(-x | -mu, 1);
  }
"
