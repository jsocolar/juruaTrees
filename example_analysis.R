library(BIOMASS)
library(brms)
data_for_AGBmc2 <- readRDS("/Users/jacobsocolar/Dropbox/Work/juruaTrees/data_for_AGBmc2.RDS")
attach(data_for_AGBmc2)

source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/AGBmc2.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/predict_heights.R")

# Example analysis: Varzea
stem_data_v <- stem_data[stem_data$Habitat == "vz", ]
stem_data_t <- stem_data[stem_data$Habitat == "tf", ]

# The workhorse function to generate the comparisons is AGBmc2, which takes the following arguments
# D: A vector of diameters. Generally stem_data$DBH_cm
# WD1: A vector of gold-standard wood densities.  Generally should use stem_data$wWD
# WD2:  NULL if comparison should use gold-standard densities.
#       Otherwise, a vector of wood densities or (if you want to incorporate errors in
#       the not-gold-standard WD measurements) a matrix.  To use wood densities based
#       on the taxonomy, use stem_data$taxonomic_WD_matrix
# H1: A vector or matrix of gold-standard heights. Generally should use the matrix
#       stem_data$gold_std_height
# errH1: If H1 is a vector, an optional vector of standard deviations representing
#       the error in the height measurements.
# H2: A vector or matrix of heights.  If using gold-standard heights in the comparison,
#       leave this as NULL.  If using biogeographic height allometries, leave this NULL as well.
#       If using modeled height-diameter data based on local allometry, use the matrix 
#       stem_data$local_allom_height
# errH2: If H2 is a vector, optionally errH2 can be a vector of associated errors.
# coord2: A 2-column matrix of coordinates for biogeographic height allometry

# Here's an example of how it works
v <- AGBmc2(D = stem_data_v$DBH_cm, 
            WD1 = stem_data_v$wWD, WD2 = stem_data_v$taxonomic_WD_matrix, 
            H1 = stem_data_v$gold_std_height, coord2 = cbind(stem_data_v$long, stem_data_v$lat))
tf <- AGBmc2(D = stem_data_t$DBH_cm, 
                  WD1 = stem_data_t$wWD, WD2 = stem_data_t$taxonomic_WD_matrix, 
                  H1 = stem_data_t$gold_std_height, coord2 = cbind(stem_data_t$long, stem_data_t$lat))

plot(rowMeans(v$agb_1) ~ rowMeans(v$agb_2))
abline(a=0, b=1)
mean(v$agb_1-v$agb_2)

agb_diff_v <- log(v$agb_2/v$agb_1)
agb_diff_t <- log(tf$agb_2/tf$agb_1)

hist(agb_diff_v[,1])
plot(agb_diff_v[,1] ~ stem_data_v$DBH_cm)

# Here's an example of how to set up the regression
stem_data2 <- rbind(stem_data_v, stem_data_t)

regression_data <- data.frame(agb_diff = c(agb_diff_v[,1], agb_diff_t[,1]), 
                              DBH = stem_data2$DBH_cm, rch = stem_data2$rch,
                              wd = stem_data2$wWD, genus = stem_data2$gen, 
                              transect = stem_data2$Transect, 
                              plot = stem_data2$plot_ID,
                              hab = stem_data2$Habitat)
rd2 <- regression_data[!is.na(regression_data$genus), ]

mymod <- brm(agb_diff ~ hab + DBH + rch + wd + (1|genus + transect + plot), data = rd2)
summary(mymod)

