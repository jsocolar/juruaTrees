library(BIOMASS)
library(cmdstanr)
library(posterior)

chave_sites <- read.csv("/Users/jacobsocolar/Desktop/chave_sites.csv")

chave_sites$bcp <- getBioclimParam(cbind(chave_sites$Lon, chave_sites$Lat))

chave <- read.csv("/Users/jacobsocolar/Downloads/Chave_GCB_Direct_Harvest_Data.csv")
chave$bcp <- chave_sites$bcp[match(chave$Site, chave_sites$Site),]

covariates <- as.matrix(cbind(chave$bcp, log(chave$DBH.cm.), log(chave$DBH.cm.)^2))

height_coord_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/height_bioclim.stan")
height_coord_data <- list(N = nrow(chave), height = chave$Total.height.m., covs = covariates)

height_coord_samples <- height_coord_model$sample(data = height_coord_data, parallel_chains = 4)
height_coord_samples$summary() # Strong match to the values from Chave et al 2014

height_coord_posterior <- as.data.frame(as_draws_df(height_coord_samples$draws()))[4*(1:1000),]

saveRDS(height_coord_posterior, "/Users/JacobSocolar/Dropbox/Work/JuruaTrees/height_coord_posterior.RDS")
