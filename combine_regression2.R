library(brms)
library(cmdstanr)
library(posterior)
source("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/example_figure.R")

# regression1_code <- make_stancode(agb_diff ~ wd_ratio_outer_inner * DBH + (1 | species + genus + family + transect + plot), data = rd2)
# 
# fileConn<-file("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression1.stan")
# writeLines(regression1_code, fileConn)
# close(fileConn)
# 
# regression1_data <- make_standata(agb_diff ~ wd_ratio_outer_inner * DBH + (1 | species + genus + family + transect + plot), data = rd2)
# class(regression1_data) <- "list"
# 
# regression_interpret_data <- list(n_stem = nrow(rd2),
#                                   n_transect = length(unique(rd2$transect)),
#                                   n_plot = length(unique(rd2$plot)),
#                                   n_species = length(unique(rd2$species)),
#                                   n_genus = length(unique(rd2$genus)),
#                                   n_family = length(unique(rd2$fam)),
#                                   agb_diff = (rd2$agb_diff - mean(rd2$agb_diff))/sd(rd2$agb_diff),
#                                   dbh = (rd2$DBH - mean(rd2$DBH))/sd(rd2$agb_diff),
#                                   wd_ratio = (rd2$wd_ratio_outer_inner - mean(rd2$wd_ratio_outer_inner))/sd(rd2$wd_ratio_outer_inner),
#                                   transect = as.integer(as.factor(rd2$transect)),
#                                   plot = as.integer(as.factor(rd2$plot)),
#                                   species = as.integer(as.factor(rd2$species)),
#                                   genus = as.integer(as.factor(rd2$genus)),
#                                   family = as.integer(as.factor(rd2$family)))
# 
# regression1_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression1.stan")
regression_transect_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression_transect.stan")
regression_noTransect_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/juruaTrees/regression_noTransect.stan")

# regression1_fit_full_warmup <- regression1_model$sample(data = as.list(regression1_data))
# regression_1_fit_short_warmup <- regression1_model$sample(data = regression1_data, iter_warmup = 500, iter_sampling = 1000)

# regression_interpret_fit_full_warmup <- regression_interpret_model$sample(data = as.list(regression_interpret_data))
# regression_interpret_fit_short_warmup <- regression_interpret_model$sample(data = regression_interpret_data, iter_warmup = 500, iter_sampling = 1000)
# 
# 
# regression_interpret_summary <- regression_interpret_fit_full_warmup$summary()
# regression_interpret_short_summary <- regression_interpret_fit_short_warmup$summary()
# 
# regression_interpret_summary
# regression_interpret_short_summary


combined_posterior <- list()
combined_diagnostics <- list()
for(j in 1:nrow(data_list)){
  print(c(nrow(data_list), j))
  
  dataset <- dataset_list[[data_list$data_set[j]]]
  
  agb_diff <- log(agb[[j]]$agb_2/agb[[j]]$agb_1)
  
  # 
  # regression_data <- data.frame(agb_diff = agb_diff[,20], 
  #                               transect = dataset$Transect, 
  #                               plot = dataset$plot_ID,
  #                               species = dataset$sp,
  #                               genus = dataset$gen, 
  #                               family = dataset$fam,
  #                               DBH = dataset$DBH_cm, 
  #                               wd_radially_weighted = dataset$wWD,
  #                               wd_unweighted = dataset$uWD,
  #                               wd_ratio_outer_inner = dataset$wd_50/dataset$wd_M
  # )
  # rd2 <- regression_data[!is.na(regression_data$species) & !is.infinite(regression_data$wd_ratio_outer_inner), ]
  # 
  
  posterior_list <- list()
  diagnostics_list <- list()
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
    rd2_i <- regression_data_i[!is.na(regression_data_i$genus) & !is.infinite(regression_data_i$wd_ratio_outer_inner), ]
    
    regression_interpret_data_i <- list(n_stem = nrow(rd2_i),
                                        n_transect = length(unique(rd2_i$transect)),
                                        n_plot = length(unique(rd2_i$plot)),
                                        n_species = length(unique(rd2_i$species)),
                                        n_genus = length(unique(rd2_i$genus)),
                                        n_family = length(unique(rd2_i$fam)),
                                        agb_diff = rd2_i$agb_diff,
                                        dbh = rd2_i$DBH,
                                        wd_ratio = rd2_i$wd_ratio_outer_inner,
                                        transect = as.integer(as.factor(rd2_i$transect)),
                                        plot = as.integer(as.factor(rd2_i$plot)),
                                        genus = as.integer(as.factor(rd2_i$genus)),
                                        family = as.integer(as.factor(rd2_i$family)))
    
    if(regression_interpret_data_i$n_plot > (regression_interpret_data_i$n_transect + 3)){
      regression_interpret_fit_i <- regression_transect_model$sample(data = as.list(regression_interpret_data_i), chains = 1, iter_warmup = 600, 
                                                                     iter_sampling = 100, adapt_delta = .95)
    }else{
      regression_interpret_fit_i <- regression_noTransect_model$sample(data = as.list(regression_interpret_data_i), chains = 1, iter_warmup = 600, 
                                                                       iter_sampling = 100, adapt_delta = .95)
    }
    
    
    
    posterior_list[[i]] <- as_draws_df(regression_interpret_fit_i$draws())
    diagnostics_list[[i]] <- as_draws_df(regression_interpret_fit_i$sampler_diagnostics())
  }
  
  combined_posterior[[j]] <- do.call(rbind, posterior_list)
  combined_diagnostics[[j]] <- do.call(rbind, diagnostics_list)
  
}

saveRDS(combined_posterior, "/Users/JacobSocolar/Desktop/combined_posterior_v3.RDS")

post_summ <- posterior::summarise_draws(combined_posterior[[1]])
print(post_summ, n=1000)

post_summ_frame <- data.frame(a_mean = rep(NA, length(combined_posterior)), a_q5 = NA, a_q95 = NA,
                              b1_dbh_mean = NA, b1_dbh_q5 = NA, b1_dbh_q95 = NA, 
                              b2_wdRatio_mean = NA, b2_wdRatio_q5 = NA, b2_wdRatio_q95 = NA, 
                              b3_interaction_mean = NA, b3_interaction_q5 = NA, b3_interaction_q95 = NA,
                              sigma_resid_mean = NA, sigma_resid_q5 = NA, sigma_resid_q95 = NA, 
                              sigma_spatial_plot_mean = NA, sigma_spatial_plot_q5 = NA,sigma_spatial_plot_q95 = NA,
                              sigma_spatial_transect_mean = NA,sigma_spatial_transect_q5 = NA,sigma_spatial_transect_q95 = NA,
                              sigma_taxonomic_genus_mean = NA, sigma_taxonomic_genus_q5 = NA,sigma_taxonomic_genus_q95 = NA,
                              sigma_taxonomic_family_mean = NA, sigma_taxonomic_family_q5 = NA, sigma_taxonomic_family_q95 = NA)

for(i in 1:length(combined_posterior)){
  print(i)
  
  d1 <- posterior::summarise_draws(combined_posterior[[i]])
  post_summ_frame$a_mean[i] <- d1$mean[d1$variable == "a"]
  post_summ_frame$a_q5[i] <- d1$q5[d1$variable == "a"]
  post_summ_frame$a_q95[i] <- d1$q95[d1$variable == "a"]
  
  post_summ_frame$b1_dbh_mean[i] <- d1$mean[d1$variable == "b1_dbh"]
  post_summ_frame$b1_dbh_q5[i] <- d1$q5[d1$variable == "b1_dbh"]
  post_summ_frame$b1_dbh_q95[i] <- d1$q95[d1$variable == "b1_dbh"]
  
  post_summ_frame$b2_wdRatio_mean[i] <- d1$mean[d1$variable == "b2_wdRatio"]
  post_summ_frame$b2_wdRatio_q5[i] <- d1$q5[d1$variable == "b2_wdRatio"]
  post_summ_frame$b2_wdRatio_q95[i] <- d1$q95[d1$variable == "b2_wdRatio"]
  
  post_summ_frame$b3_interaction_mean[i] <- d1$mean[d1$variable == "b3_interaction"]
  post_summ_frame$b3_interaction_q5[i] <- d1$q5[d1$variable == "b3_interaction"]
  post_summ_frame$b3_interaction_q95[i] <- d1$q95[d1$variable == "b3_interaction"]

  post_summ_frame$sigma_resid_mean[i] <- d1$mean[d1$variable == "sigma_resid"]
  post_summ_frame$sigma_resid_q5[i] <- d1$q5[d1$variable == "sigma_resid"]
  post_summ_frame$sigma_resid_q95[i] <- d1$q95[d1$variable == "sigma_resid"]
  
  post_summ_frame$sigma_taxonomic_genus_mean[i] <- d1$mean[d1$variable == "sigma_taxonomic_genus"]
  post_summ_frame$sigma_taxonomic_genus_q5[i] <- d1$q5[d1$variable == "sigma_taxonomic_genus"]
  post_summ_frame$sigma_taxonomic_genus_q95[i] <- d1$q95[d1$variable == "sigma_taxonomic_genus"]
  
  post_summ_frame$sigma_taxonomic_family_mean[i] <- d1$mean[d1$variable == "sigma_taxonomic_family"]
  post_summ_frame$sigma_taxonomic_family_q5[i] <- d1$q5[d1$variable == "sigma_taxonomic_family"]
  post_summ_frame$sigma_taxonomic_family_q95[i] <- d1$q95[d1$variable == "sigma_taxonomic_family"]
  
  
  ## repeat the above three lines for each parameter of interest, except for sigma_spatial_transect
  if("sigma_spatial_transect" %in% d1$variable){
    post_summ_frame$sigma_spatial_transect_mean[i] <- d1$mean[d1$variable == "sigma_spatial_transect"]
    post_summ_frame$sigma_spatial_transect_q5[i] <- d1$q5[d1$variable == "sigma_spatial_transect"]
    post_summ_frame$sigma_spatial_transect_q95[i] <- d1$q95[d1$variable == "sigma_spatial_transect"]
    
    post_summ_frame$sigma_spatial_plot_mean[i] <- d1$mean[d1$variable == "sigma_spatial_plot"]
    post_summ_frame$sigma_spatial_plot_q5[i] <- d1$q5[d1$variable == "sigma_spatial_plot"]
    post_summ_frame$sigma_spatial_plot_q95[i] <- d1$q95[d1$variable == "sigma_spatial_plot"]
  }else{
    post_summ_frame$sigma_spatial_plot_mean[i] <- d1$mean[d1$variable == "sigma_plot"]
    post_summ_frame$sigma_spatial_plot_q5[i] <- d1$q5[d1$variable == "sigma_plot"]
    post_summ_frame$sigma_spatial_plot_q95[i] <- d1$q95[d1$variable == "sigma_plot"]
  }
}


View(post_summ_frame)

class(agb)
