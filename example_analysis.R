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
            WD1 = stem_data_v$wWD, WD2 = stem_data_v$wWD_s, 
            H1 = stem_data_v$gold_std_height, H2=NULL)

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

# Here's an example of how to set up the regression for one column.  Ultimately this regression
# needs to be iterated over columns (maybe a sample of 50 or 100 columns rather than all 1000),
# and the resulting posteriors need to be combined for inference.
stem_data2 <- rbind(stem_data_v, stem_data_t)
stem_data2$wd_ratio_outer_inner <- stem_data2$wd_50/stem_data2$wd_M

regression_data <- data.frame(agb_diff = c(agb_diff_v[,1], agb_diff_t[,1]), 
                              transect = stem_data2$Transect, 
                              plot = stem_data2$plot_ID,
                              hab = stem_data2$Habitat,
                              species = stem_data2$sp,
                              genus = stem_data2$gen, 
                              family = stem_data2$fam,
                              DBH = stem_data2$DBH_cm, 
                              wd_radially_weighted = stem_data2$wWD,
                              wd_unweighted = stem_data2$uWD,
                              wd_ratio_outer_inner = stem_data2$wd_ratio_outer_inner
)
sum(is.na(regression_data$family))
sum(is.na(regression_data$genus))
sum(is.na(regression_data$species))
sum(is.na(regression_data$rch) | is.infinite(regression_data$rch))
sum(is.infinite(regression_data$rch) | is.na(regression_data$species))
rd2 <- regression_data[!is.na(regression_data$species) & !is.infinite(regression_data$wd_ratio_outer_inner), ]

pairs(rd2[,c("agb_diff", "DBH", "wd_radially_weighted", "wd_unweighted", "wd_ratio_outer_inner")])
plot(log(rd2$wd_unweighted/rd2$wd_radially_weighted) ~ log(rd2$wd_ratio_outer_inner))
plot(log(rd2$wd_unweighted/rd2$wd_radially_weighted) ~ log(rd2$wd_ratio_outer_inner),
     xlim = c(-.5, .5), ylim = c(-.2,.1))
summary(lm(log(rd2$wd_unweighted/rd2$wd_radially_weighted) ~ log(rd2$wd_ratio_outer_inner)))

make_stancode(agb_diff ~ hab * DBH + wd_ratio_outer_inner * DBH + (1| species + genus + family + transect + plot), data = rd2)
make_standata(agb_diff ~ hab * DBH + wd_ratio_outer_inner * DBH + (1| species + genus + family + transect + plot), data = rd2)

mymod <- brm(agb_diff ~ hab * DBH + wd_ratio_outer_inner * DBH + (1| species + genus + family + transect + plot), data = rd2)
summary(mymod)

