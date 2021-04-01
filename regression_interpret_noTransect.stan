//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> n_stem;
  int<lower=1> n_transect;
  int<lower=1> n_plot;
  int<lower=1> n_species;
  int<lower=1> n_genus;
  int<lower=1> n_family;
  
  vector[n_stem] agb_diff;
  vector[n_stem] dbh;
  vector[n_stem] wd_ratio;
  int<lower=1, upper=n_transect> transect[n_stem];
  int<lower=1, upper=n_plot> plot[n_stem];
  int<lower=1, upper=n_species> species[n_stem];
  int<lower=1, upper=n_genus> genus[n_stem];
  int<lower=1, upper=n_family> family[n_stem];
}

parameters {
  real a;
  real b1_dbh;
  real b2_wdRatio;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_taxonomic;
  real<lower=0, upper=1> p_taxonomic_species;
  real<lower=0, upper=1-p_taxonomic_species> p_taxonomic_genus;
  vector[n_plot] a_plot_raw;
  vector[n_species] a_species_raw;
  vector[n_genus] a_genus_raw;
  vector[n_family] a_family_raw;
}

transformed parameters{
  real<lower = 0, upper = 1> sigma_taxonomic_species = sqrt(sigma_taxonomic^2 * p_taxonomic_species);
  real<lower = 0, upper = 1> sigma_taxonomic_genus = sqrt(sigma_taxonomic^2 * p_taxonomic_genus);
  real<lower = 0, upper = 1> sigma_taxonomic_family = sqrt(sigma_taxonomic^2 * (1 - p_taxonomic_genus - p_taxonomic_species));
  
  vector[n_plot] a_plot = a_plot_raw * sigma_plot;
  vector[n_species] a_species = a_species_raw * sigma_taxonomic_species;
  vector[n_genus] a_genus = a_genus_raw * sigma_taxonomic_genus;
  vector[n_family] a_family = a_family_raw * sigma_taxonomic_family;
}

model {
  
  vector[n_stem] mu = a + a_plot[plot] + 
            a_species[species] + a_genus[genus] + a_family[family] +
            b1_dbh * dbh + b2_wdRatio * wd_ratio;
  agb_diff ~ normal(mu, sigma_resid);
  
  // priors
  a ~ normal(0,4);
  b1_dbh ~ normal(0,4);
  b2_wdRatio ~ normal(0,4);
  a_plot_raw ~ std_normal();
  a_species_raw ~ std_normal();
  a_genus_raw ~ std_normal();
  a_family_raw ~ std_normal();
  sigma_resid ~ normal(0, 3);
  sigma_plot ~ normal(0, 3);
  sigma_taxonomic ~ normal(0, 3);
}

