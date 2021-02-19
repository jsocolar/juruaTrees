data{
  int N; // number of stems
  vector[N] height; // tree heights
  matrix[N,5] covs; // 3 environmental covariates, log_dbh, log_dbh^2
}
transformed data{
  vector[N] log_height = log(height);
}
parameters{
  real intercept;
  vector[5] beta;
  real<lower=0, upper=30> sigma;
}
model{
  vector[N] mu = intercept + covs*beta;
  log_height ~ normal(mu, sigma);
  
  // Priors
  intercept ~ normal(0,10);
  beta ~ normal(0,10);
  // Implicit uniform 0,30 prior on sigma
}