# Load up the necessary functions
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/AGBmc2.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/predict_heights.R")

height_coord_posterior <- readRDS("/Users/JacobSocolar/Dropbox/Work/JuruaTrees/Data_objects/height_coord_posterior.RDS")

# Import the data list and edit a couple of typos
data_list <- readxl::read_xlsx("/Users/JacobSocolar/Dropbox/Work/JuruaTrees/Data_objects/percent_table.xlsx")
unique(data_list$data_set)
data_list$data_set[data_list$data_set == "Hestem_data_v"] <- "HEstem_data_v"
data_list$data_set[data_list$data_set == "Hefam_v"] <- "HEfam_v"
data_list$data_set[data_list$data_set == "Hefam_t"] <- "HEfam_t"
# remove the rows that imply comparisons of the gold standard to itself
data_list <- data_list[-c(65,66),]

# Import a list of all the necessary datasets
dataset_list <- readRDS("/Users/JacobSocolar/Dropbox/Work/JuruaTrees/Data_objects/dataset_list.RDS")

# Do the agb calculations
agb <- list()
a <- txtProgressBar(min = 1, max = 70, style=3)
for(i in 1:nrow(data_list)){
  setTxtProgressBar(a, i)
  dataset <- dataset_list[[data_list$data_set[i]]]
    if(data_list$H[i] == "coord"){
      agb[[i]] <- AGBmc2(D = dataset$DBH_cm,
                       WD1 = dataset$wWD, 
                       WD2 = eval(parse(text=paste0("dataset$", data_list$WD[i]))),
                       H1 = dataset$gold_std_height,
                       coord2 = as.matrix(cbind(dataset$long, dataset$lat))   ###NEED TO CONFIRM THAT COORD SHOULD BE LONG-LAT AND NOT THE REVERSE
                      )
  }else{
    agb[[i]] <- AGBmc2(D = dataset$DBH_cm,
                       WD1 = dataset$wWD, 
                       WD2 = eval(parse(text=paste0("dataset$", data_list$WD[i]))),
                       H1 = dataset$gold_std_height,
                       H2 = eval(parse(text=paste0("dataset$", data_list$H[i])))
                      )
  }
}

# Subset the agb calculations by habitat
agb_tf <- agb[data_list$Habitat=="tf"]
agb_vz <- agb[data_list$Habitat=="vz"]

# Plot results in Terra Firme
plot(density(colSums(100*agb_tf[[1]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[1]]$agb_1))), 
     xlim =c(-60,35), ylim = c(0,3), ylab = "probability density", yaxt = "n",
     xlab = "AGB difference (percent)", main = "Terra firme")
for(i in 2:length(agb_tf)){
  lines(density(colSums(100*agb_tf[[i]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[i]]$agb_1))))
}

# Plot results in Varzea
plot(density(colSums(100*agb_vz[[1]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[1]]$agb_1))), 
     xlim =c(-60,35), ylim = c(0,3), ylab = "probability density", yaxt = "n", 
     xlab = "AGB difference (percent)", main = "Varzea")
for(i in 2:length(agb_vz)){
  lines(density(colSums(100*agb_vz[[i]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[i]]$agb_1))))
}

# Examples for how to make the plots pretty.
# The super tall peak is the uWD vs wWD comparison with everything else held constant.  There's no 
# uncertainty in the comparison, so the posterior distribution is a point mass concentrated
# at that exact point.  To make things a bit prettier, we can remove it, for example:
plot(density(colSums(100*agb_vz[[1]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[1]]$agb_1))), 
     xlim =c(-60,35), ylim = c(0,1.7), ylab = "probability density", yaxt = "n", 
     xlab = "AGB difference (percent)", main = "Varzea")
for(i in 2:length(agb_vz)){
  if(i != 32){
    lines(density(colSums(100*agb_vz[[i]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[i]]$agb_1))))
  }
}

# If we want to highlight a particular subset of methods, we could do something like
# Plot results in Varzea
plot(density(colSums(100*agb_vz[[1]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[1]]$agb_1))), 
     xlim =c(-60,35), ylim = c(0,1.7), ylab = "probability density", yaxt = "n", 
     xlab = "AGB difference (percent)", main = "Varzea", col = "gray")
for(i in c(2:31, 33:length(agb_vz))){
  lines(density(colSums(100*agb_vz[[i]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[i]]$agb_1))), col = "gray")
}
lines(density(colSums(100*agb_vz[[1]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[1]]$agb_1))), col = "firebrick2", lwd = 2)
lines(density(colSums(100*agb_vz[[2]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[2]]$agb_1))), col = "blue", lwd = 2)
lines(density(colSums(100*agb_vz[[20]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[20]]$agb_1))), col = "orange", lwd = 2)
lines(density(colSums(100*agb_vz[[30]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[30]]$agb_1))), col = "purple", lwd = 2)
lines(density(colSums(100*agb_vz[[33]]$stemwise_AGB_diff)/mean(colSums(agb_vz[[33]]$agb_1))), col = "forestgreen", lwd = 2)

# And compare to the same methods in terra firme
plot(density(colSums(100*agb_tf[[1]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[1]]$agb_1))), 
     xlim =c(-60,35), ylim = c(0,1.7), ylab = "probability density", yaxt = "n", 
     xlab = "AGB difference (percent)", main = "Terra firme", col = "gray")
for(i in c(2:31, 33:length(agb_tf))){
  lines(density(colSums(100*agb_tf[[i]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[i]]$agb_1))), col = "gray")
}
lines(density(colSums(100*agb_tf[[1]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[1]]$agb_1))), col = "firebrick2", lwd = 2)
lines(density(colSums(100*agb_tf[[2]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[2]]$agb_1))), col = "blue", lwd = 2)
lines(density(colSums(100*agb_tf[[20]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[20]]$agb_1))), col = "orange", lwd = 2)
lines(density(colSums(100*agb_tf[[30]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[30]]$agb_1))), col = "purple", lwd = 2)
lines(density(colSums(100*agb_tf[[33]]$stemwise_AGB_diff)/mean(colSums(agb_tf[[33]]$agb_1))), col = "forestgreen", lwd = 2)



