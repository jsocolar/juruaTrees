library(brms)
library(cmdstanr)
library(posterior)


regression1_code <- make_stancode(agb_diff ~ wd_ratio_outer_inner * DBH + (1 | species + genus + family + transect + plot), data = rd2)

fileConn<-file("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression1.stan")
writeLines(regression1_code, fileConn)
close(fileConn)

regression1_data <- make_standata(agb_diff ~ wd_ratio_outer_inner * DBH + (1 | species + genus + family + transect + plot), data = rd2)
class(regression1_data) <- "list"

regression_interpret_data <- list(n_stem = nrow(rd2),
                                  n_transect = length(unique(rd2$transect)),
                                  n_plot = length(unique(rd2$plot)),
                                  n_species = length(unique(rd2$species)),
                                  n_genus = length(unique(rd2$genus)),
                                  n_family = length(unique(rd2$fam)),
                                  agb_diff = (rd2$agb_diff - mean(rd2$agb_diff))/sd(rd2$agb_diff),
                                  dbh = (rd2$DBH - mean(rd2$DBH))/sd(rd2$agb_diff),
                                  wd_ratio = (rd2$wd_ratio_outer_inner - mean(rd2$wd_ratio_outer_inner))/sd(rd2$wd_ratio_outer_inner),
                                  transect = as.integer(as.factor(rd2$transect)),
                                  plot = as.integer(as.factor(rd2$plot)),
                                  species = as.integer(as.factor(rd2$species)),
                                  genus = as.integer(as.factor(rd2$genus)),
                                  family = as.integer(as.factor(rd2$family)))

regression1_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression1.stan")
regression_interpret_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression_interpret.stan")

regression1_fit_full_warmup <- regression1_model$sample(data = as.list(regression1_data))
regression_1_fit_short_warmup <- regression1_model$sample(data = regression1_data, iter_warmup = 500, iter_sampling = 1000)

regression1_summary <- regression1_fit_full_warmup$summary()
regression1_short_summary <- regression_1_fit_short_warmup$summary()
regression1_summary[regression1_summary$variable == "b_Intercept",]

regression1_summary
regression1_short_summary


for(j in 1:nrow(data_list)){
  print(c(nrow(data_list), j))
  
  dataset <- dataset_list[[data_list$data_set[j]]]
  
  agb_diff <- log(agb[[j]]$agb_2/agb[[j]]$agb_1)
  
  
  regression_data <- data.frame(agb_diff = agb_diff[,20], 
                                transect = dataset$Transect, 
                                plot = dataset$plot_ID,
                                species = dataset$sp,
                                genus = dataset$gen, 
                                family = dataset$fam,
                                DBH = dataset$DBH_cm, 
                                wd_radially_weighted = dataset$wWD,
                                wd_unweighted = dataset$uWD,
                                wd_ratio_outer_inner = dataset$wd_50/dataset$wd_M
  )
  rd2 <- regression_data[!is.na(regression_data$species) & !is.infinite(regression_data$wd_ratio_outer_inner), ]
  
  
posterior_list <- list()
for(i in 1:50){
  print(i)
  regression_data_i <- data.frame(agb_diff = agb_diff[,i*20], 
                                transect = dataset$Transect, 
                                plot = dataset$plot_ID,
                                species = dataset$sp,
                                genus = dataset$gen, 
                                family = dataset$fam,
                                DBH = dataset$DBH_cm, 
                                wd_radially_weighted = dataset$wWD,
                                wd_unweighted = dataset$uWD,
                                wd_ratio_outer_inner = dataset$wd_50/dataset$wd_M
                              )
  rd2_i <- regression_data_i[!is.na(regression_data_i$species) & !is.infinite(regression_data_i$wd_ratio_outer_inner), ]
  regression1_data_i <- make_standata(agb_diff ~ wd_ratio_outer_inner * DBH + (1| species + genus + family + transect + plot), data = rd2_i)
  class(regression1_data_i) <- "list"
  regression1_fit_i <- regression1_model$sample(data = as.list(regression1_data_i), chains = 1, iter_warmup = 600, 
                                                iter_sampling = 100, adapt_delta = .95)
  
  posterior_list[[i]] <- as_draws_df(regression1_fit_i$draws())
}

combined_posterior[[j]] <- do.call(rbind, posterior_list)
}
