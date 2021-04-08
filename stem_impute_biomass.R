# Load up the necessary functions
library(BIOMASS)
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/AGBmc2.R")

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


# Get modeled heights for all trees
v_model <- modelHD(trees$DBH_cm[trees$Habitat == "vz"], H = trees$H_M[trees$Habitat == "vz"], method = "michaelis")
t_model <- modelHD(trees$DBH_cm[trees$Habitat == "tf"], H = trees$H_M[trees$Habitat == "tf"], method = "michaelis")
trees$height_v_model <- BIOMASS:::predictHeight(trees$DBH_cm, v_model, err = F)
trees$height_t_model <- BIOMASS:::predictHeight(trees$DBH_cm, t_model, err = F)


dim(trees)
sum(trees$Habitat=="vz")
trees$wd_best <- trees$height_best <- NA
for(i in 1:nrow(trees)){
  if(!is.na(trees$wWD[i])){
    trees$wd_best[i] <- trees$wWD[i]
  }else if(trees$Habitat[i] == "vz"){
    
    if(!is.na(trees$wWD_g[i])){
      trees$wd_best[i] <- trees$wWD_g[i]
    }else if(!is.na(trees$uWD_s[i])){
      trees$wd_best[i] <- trees$uWD_s[i]
    }else if(!is.na(trees$uWD_f[i])){
      trees$wd_best[i] <- trees$uWD_f[i]
    }else{
      trees$wd_best[i] <- mean(trees$wWD[trees$plot_ID == trees$plot_ID[i]], na.rm = T)
    }
  }else if(trees$Habitat[i] == "tf"){
    if(!is.na(trees$wWD_s[i])){
      trees$wd_best[i] <- trees$wWD_s[i]
    }else if(!is.na(trees$wWD_g[i])){
      trees$wd_best[i] <- trees$wWD_g[i]
    }else if(!is.na(trees$wWD_f[i])){
      trees$wd_best[i] <- trees$wWD_f[i]
    }else{
      trees$wd_best[i] <- mean(trees$wWD[trees$plot_ID == trees$plot_ID[i]], na.rm = T)
    }
  }
  
  if(!is.na(trees$H_M[i])){
    trees$height_best[i] <- trees$H_M[i] 
  }else if(trees$Habitat[i] == "vz"){
    trees$height_best[i] <- trees$height_v_model[i]
  }else if(trees$Habitat[i] == "tf"){
    trees$height_best[i] <- trees$height_t_model[i]
  }
  
}

varzea <- trees[trees$Habitat == "vz", c("DBH_cm", "wd_best", "height_best")]
tf <- trees[trees$Habitat == "tf", c("DBH_cm", "wd_best", "height_best")]

varzea_agb <- BIOMASS::AGBmonteCarlo(D = varzea$DBH_cm, WD = varzea$wd_best, errWD = rep(0, nrow(varzea)), H = varzea$height_best, errH = 0)
tf_agb <- BIOMASS::AGBmonteCarlo(D = tf$DBH_cm, WD = tf$wd_best, errWD = rep(0, nrow(tf)), H = tf$height_best, errH = 0)

varzea_agb$meanAGB
tf_agb$meanAGB
