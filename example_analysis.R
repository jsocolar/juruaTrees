library(BIOMASS)
library(cmdstanr)

trees <- readxl::read_excel("/Users/jacobsocolar/Downloads/treeDF.xlsx", sheet = 2)
withheight <- trees[!is.na(trees$H_M1) & !is.na(trees$H_M2) & !is.na(trees$H_M3), ]

# # Here we fit a naive error model to estimate the error in the field-based tree height measurements.
# # This model will require substantial refinement if we actually want to use it in the paper.
# mod <- cmdstan_model("/Users/jacobsocolar/Desktop/treeheighterr.stan")
# samps <- mod$sample(data = list(y1 = withheight$H_M1, y2 = withheight$H_M2, y3 = withheight$H_M3, N=nrow(withheight)), cores = 4)
# print(posterior::summarise_draws(samps$draws()[,,"sigma"]))

# Get trees with available height measurements
hwd <- withheight[!is.na(withheight$wWD),]

# Analyze terra firme first
hwd_tf <- hwd[hwd$Habitat == "tf", ]
tf <- AGBmc2(D = hwd_tf$DBH_cm, WD1 = hwd_tf$wWD, H1 = hwd_tf$H_M, 
       coord2 = cbind(hwd_tf$long, hwd_tf$lat))
plot(rowMeans(tf$agb_1) ~ rowMeans(tf$agb_2))
abline(a=0, b=1)
mean(tf$agb_1-tf$agb_2)

# Now varzea
hwd_vz <- hwd[hwd$Habitat == "vz", ]
vz <- AGBmc2(D = hwd_vz$DBH_cm, WD1 = hwd_vz$wWD, H1 = hwd_vz$H_M, 
             coord2 = cbind(hwd_vz$long, hwd_vz$lat))
plot(rowMeans(vz$agb_1) ~ rowMeans(vz$agb_2))
abline(a=0, b=1)
mean(vz$agb_1-vz$agb_2)



vz <- AGBmc2(D = hwd_vz$DBH_cm, WD1 = hwd_vz$wWD, H1 = hwd_vz$H_M, WD2 = hwd_vz$uWD,
             coord2 = cbind(hwd_vz$long, hwd_vz$lat))
plot(rowMeans(vz$agb_1) ~ rowMeans(vz$agb_2))
abline(a=0, b=1)
hist((rowMeans(vz$agb_1) - rowMeans(vz$agb_2))/rowMeans(vz$agb_1))
sum(rowMeans(vz$agb_1))
sum(rowMeans(vz$agb_2))

vt1 <- AGBmonteCarlo(D=hwd_vz$DBH_cm, WD=hwd_vz$wWD, errWD = rep(.001, nrow(hwd_vz)), H = hwd_vz$H_M, errH=0)
vt2 <- AGBmonteCarlo(D=hwd_vz$DBH_cm, WD=hwd_vz$wWD, errWD = rep(.001, nrow(hwd_vz)), coord = cbind(hwd_vz$long, hwd_vz$lat))
vt1$meanAGB
sum(rowMeans(vz$agb_1))
vt2$meanAGB
sum(rowMeans(vz$agb_2))
sum(rowMeans(vt2$AGB_simu))
