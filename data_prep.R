library(BIOMASS)
library(cmdstanr)

set.seed(123)

myrtruncnorm <- function(n, lower = -1, upper = 1, mean = 0, 
                         sd = 1) {
  qnorm(runif(n, pnorm(lower, mean = mean, sd = sd), pnorm(upper, 
                                                           mean = mean, sd = sd)), mean = mean, sd = sd)
}

# get coordinate height model
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/height_coord_posterior.R")

# Read in the tree data
trees <- readxl::read_excel("/Users/jacobsocolar/Downloads/treeDF.xlsx", sheet = 2)
# Rename problem transects
trees$Transect[trees$Transect == "OUP"] <- "OUP_A"
trees$Transect[trees$Transect == "XNZ"] <- "XZN_A"

##### Insert Yennie's code here #####
# CORRECT wd_M, wWD & uWD for Tree_Tag = "6809"
stem_data2<-trees
# wd_M [869,30]
stem_data2[869,30]<-(10*(stem_data2[869,30])) 
# wWD [869,27]
stem_data2[869,27]<-0.461534032828081 
# uWD [869,34]
stem_data2[869,34]<-0.446929148153457 

# subset for Habitat:
stem_data2tf<-subset(stem_data2, Habitat=="tf")
stem_data2vz<-subset(stem_data2, Habitat=="vz") # Tree_Tag = 6809 is in vz

# SP:
stem_data3vz<-subset(stem_data2vz, sp=="Leguminosae_Pterocarpus_rohrii")
# CORRECT wWD_s for sp="Leguminosae_Pterocarpus_rohrii":
vz_wWDs<- aggregate(data=stem_data3vz, wWD~sp, FUN = mean)
stem_data2vz$wWD_s[
  stem_data2vz$sp == "Leguminosae_Pterocarpus_rohrii"] <- vz_wWDs$wWD
# CORRECT uWD_s for sp="Leguminosae_Pterocarpus_rohrii":
vz_uWDs<- aggregate(data=stem_data3vz, uWD~sp, FUN = mean)
stem_data2vz$uWD_s[
  stem_data2vz$sp == "Leguminosae_Pterocarpus_rohrii"] <- vz_uWDs$uWD

# GEN:
stem_data3vz<-subset(stem_data2vz, gen=="Leguminosae_Pterocarpus")
# CORRECT wWD_g for gen="Leguminosae_Pterocarpus":
vz_wWDg<- aggregate(data=stem_data3vz, wWD~gen, FUN = mean)
stem_data2vz$wWD_g[
  stem_data2vz$gen == "Leguminosae_Pterocarpus"] <- vz_wWDg$wWD
# CORRECT uWD_g for gen="Leguminosae_Pterocarpus":
vz_uWDg<- aggregate(data=stem_data3vz, uWD~gen, FUN = mean)
stem_data2vz$uWD_g[
  stem_data2vz$gen == "Leguminosae_Pterocarpus"] <- vz_uWDg$uWD

# FAM:
stem_data3vz<-subset(stem_data2vz, fam=="Leguminosae")
# CORRECT wWD_f for fam="Leguminosae":
vz_wWDf<- aggregate(data=stem_data3vz, wWD~fam, FUN = mean)
stem_data2vz$wWD_f[
  stem_data2vz$fam == "Leguminosae"] <- vz_wWDf$wWD
# CORRECT uWD_f for fam="Leguminosae":
vz_uWDf<- aggregate(data=stem_data3vz, uWD~fam, FUN = mean)
stem_data2vz$uWD_f[
  stem_data2vz$fam == "Leguminosae"] <- vz_uWDf$uWD

# JOIN stem_data2vz & stem_data2tf:
trees <- rbind(stem_data2vz, stem_data2tf)
##### End Yennie's insertion #####




# Get just trees with heights
withheight <- trees[!is.na(trees$H_M1) & !is.na(trees$H_M2) & !is.na(trees$H_M3), ]
# Create a stable stem identifier (that will remain stable even as we subset the data)
withheight$stem_id <- paste0(c(1:nrow(withheight)), "_")

# Extract the families, genera, and species
fgen <- do.call(rbind, strsplit(withheight$gen, "_"))
fsp <- do.call(rbind, strsplit(withheight$sp, "_"))
withheight$Family <- withheight$fam
withheight$Genus <- fgen[,2]
withheight$Species <- fsp[,3]
withheight$taxID <- paste(withheight$Family, withheight$Genus, withheight$Species, sep = "_")

# Get the "taxonomic wood densities"
twd <- getWoodDensity(genus = withheight$Genus, species = withheight$Species, family = withheight$Family, region = "World", stand = withheight$plot_ID)

# There are multiple ways ways to handle these taxonomic wood densities. We could assume that the errors are:
# 1) uncorrelated
# 2) correlated within species (particularly if the estimate is derived at the genus or family level), or 
# 3) correlated within species-by-habitat

# (3) probably is overkill (within-species sds don't appreciably differ when varzea and tf are lumped)
sigma(lm(withheight$wWD ~ withheight$sp))
sigma(lm(withheight$wWD ~ withheight$sp * withheight$Habitat))

# (2) probably isn't overkill; the within-species sds from the data are smaller than the genus- and family-level sds.
# So here we go. We decompose the genus and family level variance into a higher-order part and the species part.
twd$taxID <- withheight$taxID
twd$stem_id <- withheight$stem_id
twd2 <- twd[!duplicated(twd$taxID),]

sp_sd <- unique(twd2$sdWD[twd2$levelWD == "species"])
gen_sd <- (unique(twd2$sdWD[twd2$levelWD == "genus"])^2 - sp_sd^2)^.5
fam_sd <- (unique(twd2$sdWD[twd2$levelWD == "family"])^2 - sp_sd^2)^.5
tax_sds <- c(fam_sd, gen_sd, 0, 0)

tax_level <- rep(0, nrow(twd2))
tax_level[twd2$levelWD == "family"] <- 1
tax_level[twd2$levelWD == "genus"] <- 2
tax_level[twd2$levelWD == "species"] <- 3
tax_level[!(twd2$levelWD %in% c("family", "genus", "species"))] <- 4
wd_list <- list()
for(i in 1:1000){
  species_means <- myrtruncnorm(nrow(twd2), lower = .1, upper = 1.4, mean = twd2$meanWD, sd = tax_sds[tax_level]) # get samples of the species-level means
  w <- species_means[match(withheight$taxID, twd2$taxID)]  # index into the taxa in withheight to assign the means
  w2 <- myrtruncnorm(nrow(withheight), lower = .1, upper = 1.4, mean = w, sd = sp_sd) # sample the species sds
  # deal with the taxa with no family designation
  w2[is.na(withheight$fam)] <- myrtruncnorm(sum(is.na(withheight$fam)), lower = .1, upper = 1.4,
                                                       mean = twd$meanWD[is.na(withheight$fam)], 
                                                       sd = twd$sdWD[is.na(withheight$fam)])
  wd_list[[i]] <- w2
}
WD2 <- do.call(cbind, wd_list)
withheight$taxonomic_WD_matrix <- WD2

##### Estimate the precision of tree height measurements and get matrix of gold-standard heights #####
withheight$mean_height <- rowMeans(withheight[, c("H_M1", "H_M2", "H_M3")])
withheight$h1_diff <- withheight$H_M1 - withheight$mean_height
# obvious mean-variance relationship (as expected)
plot(withheight$h1_diff ~ withheight$mean_height)

mod <- cmdstan_model("/Users/jacobsocolar/Desktop/treeheighterr.stan")
# Here, we model the height measurements as being drawn from student-t distributions with 3 degrees of freedom
# (to allow for very heavy-tailed shapes compared to normal distributions)
samps <- mod$sample(data = list(y1 = withheight$H_M1, y2 = withheight$H_M2, y3 = withheight$H_M3, N=nrow(withheight)))
a <- samps$summary()
max(a$rhat)
min(a$ess_bulk)
min(a$ess_tail)
height_draws <- posterior::as_draws_df(samps$draws())[4*(1:1000),]
height_draws <- height_draws[,grep("mu\\[", names(height_draws))]
withheight$gold_std_height <- t(height_draws)

# Local DBH-height allometries in Varzea and Terra Firme separately
# ten largest trees plus stratified sample of 40 additional trees
wh_v <- withheight[withheight$Habitat == "vz", ]
wh_v <- wh_v[order(wh_v$DBH_cm, decreasing = T),]
v_sample <- wh_v[c(1:10, 11 + floor(c(1:40)*(nrow(wh_v)-11)/40)), ]
wh_t <- withheight[withheight$Habitat == "tf", ]
wh_t <- wh_t[order(wh_t$DBH_cm, decreasing = T),]
t_sample <- wh_t[c(1:10, 11 + floor(c(1:40)*(nrow(wh_t)-11)/40)), ]

modelHD(v_sample$DBH_cm, H = v_sample$H_M)
modelHD(t_sample$DBH_cm, H = t_sample$H_M)

v_model <- modelHD(v_sample$DBH_cm, H = v_sample$H_M, method = "michaelis")
v_height_list <- list()
for(i in 1:1000){
  v_height_list[[i]] <- BIOMASS:::predictHeight(wh_v$DBH_cm, v_model, err = T)
}
v_height <- do.call(cbind, v_height_list)
for(i in 1:50){
  v_height[c(1:10, 11 + floor(c(1:40)*(nrow(wh_v)-11)/40))[i], ] <- wh_v$gold_std_height[c(1:10, 11 + floor(c(1:40)*(nrow(wh_v)-11)/40))[i], ]
}

t_model <- modelHD(t_sample$DBH_cm, H = t_sample$H_M, method = "michaelis")
t_height_list <- list()
for(i in 1:1000){
  t_height_list[[i]] <- BIOMASS:::predictHeight(wh_t$DBH_cm, t_model, err = T)
}
t_height <- do.call(cbind, t_height_list)
for(i in 1:50){
  t_height[c(1:10, 11 + floor(c(1:40)*(nrow(wh_t)-11)/40))[i], ] <- wh_t$gold_std_height[c(1:10, 11 + floor(c(1:40)*(nrow(wh_t)-11)/40))[i], ]
}

wh_v$local_allom_height <- v_height
wh_t$local_allom_height <- t_height

wh_new <- rbind(wh_v, wh_t)

# Get trees with gold-standard wood density measures
hwd <- wh_new[!is.na(wh_new$wWD),]

data_for_AGBmc2 <- list(stem_data = hwd, height_coord_posterior = height_coord_posterior)

saveRDS(data_for_AGBmc2, file = "/Users/jacobsocolar/Dropbox/Work/juruaTrees/data_for_AGBmc2.RDS")
