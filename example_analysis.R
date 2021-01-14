library(BIOMASS)
library(cmdstanr)

trees <- readxl::read_excel("/Users/jacobsocolar/Downloads/treeDF.xlsx", sheet = 2)
withheight <- trees[!is.na(trees$H_M1) & !is.na(trees$H_M2) & !is.na(trees$H_M3), ]

# Here we fit a naive error model to estimate the error in the field-based tree height measurements.
# This model will require substantial refinement if we actually want to use it in the paper.
mod <- cmdstan_model("/Users/jacobsocolar/Desktop/treeheighterr.stan")
samps <- mod$sample(data = list(y1 = withheight$H_M1, y2 = withheight$H_M2, y3 = withheight$H_M3, N=nrow(withheight)), cores = 4)
print(posterior::summarise_draws(samps$draws()[,,"sigma"]))

# Get trees with available height measurements
hwd <- withheight[!is.na(withheight$wWD),]

# Analyze terra firme first
hwd_tf <- hwd[hwd$Habitat == "tf", ]
# predict "gold-standard" biomass estimates and their full uncertianty
true_biomass_tf <- AGBmonteCarlo(D = hwd_tf$DBH_cm, WD = hwd_tf$wWD, errWD = rep(.01, nrow(hwd_tf)), H = hwd_tf$H_M, errH = rep(3,nrow(hwd_tf)))
# predict biomass estimates (and their full uncertainty) based on a DBH-height model based on the local coordinates
no_height_tf <- AGBmonteCarlo(D = hwd_tf$DBH_cm, WD = hwd_tf$wWD, errWD = rep(.01, nrow(hwd_tf)), coord = cbind(hwd_tf$long, hwd_tf$lat))

# get full posterior distribution for the stem-specific differences based on the two models
differences_tf <- no_height_tf[[5]]-true_biomass_tf[[5]]
# measure the bias
bias_tf <- colMeans(differences_tf)
hist(bias_tf, main = "tf posterior bias (per stem)")  # bias is essentially nonexistent
abline(v=bias_tf[order(bias_tf)][round(length(bias_tf)*.05)]) # 90% CI bound
# measure the coverage
coverage_tf <- apply(differences_tf, 1, function(x){sum(x > 1)})/1000
hist(coverage_tf, breaks = 60)   # despite low bias, coverage is quite poor.
abline(v=.025)  # 95% CI bound

# repeat for varzea
hwd_vz <- hwd[hwd$Habitat == "vz", ]
true_biomass_vz <- AGBmonteCarlo(D = hwd_vz$DBH_cm, WD = hwd_vz$wWD, errWD = rep(.01, nrow(hwd_vz)), H = hwd_vz$H_M, errH = rep(3,nrow(hwd_vz)))
no_height_vz <- AGBmonteCarlo(D = hwd_vz$DBH_cm, WD = hwd_vz$wWD, errWD = rep(.01, nrow(hwd_vz)), coord = cbind(hwd_vz$long, hwd_vz$lat))
differences_vz <- no_height_vz[[5]]-true_biomass_vz[[5]]
bias_vz <- colMeans(differences_vz)
hist(bias_vz, main = "vz posterior bias (per stem)") # there is clear bias
abline(v=bias_vz[order(bias_vz)][round(length(bias_vz)*.05)])
hist(bias_vz/mean(true_biomass_vz[[5]])) # bias is on the order of 15%
coverage_vz <- apply(differences_vz, 1, function(x){sum(x > 1)})/1000
hist(coverage_vz, breaks = 60) # in addition to bias, coverage is poor.
abline(v=.025)

