carbon <- readRDS('/Users/jacobsocolar/Downloads/carbdat.RDS')
newdata <- list(n_habitat = carbon$nhab,
                n_point = carbon$nsite,
                n_cluster = carbon$nclus,
                lcarbon_point = scale(log(carbon$carb))[,1],
                elevation = scale(carbon$elev)[,1],
                cluster = carbon$clus,
                habitat = carbon$habitat
)
library(cmdstanr)
carbon_model <- cmdstan_model("/Users/jacobsocolar/Desktop/jorgen_model.stan")
carbon_samples <- carbon_model$sample(data = newdata, parallel_chains = 3, chains = 3)
carbon_samples$summary()