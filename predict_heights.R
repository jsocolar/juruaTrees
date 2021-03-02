predict_heights <- function(D, coord){
  len = length(D)
  posterior_E <- as.matrix(height_coord_posterior[,c("beta[1]", "beta[2]", "beta[3]")]) %*% t(as.matrix(BIOMASS::getBioclimParam(coord)))
  lD <- replicate(1000,log(D))
  hmeans <- t(replicate(len, height_coord_posterior$intercept)) + t(posterior_E) + 
    sweep(lD, MARGIN=2, height_coord_posterior$`beta[4]`, `*`) +
    sweep(lD^2, MARGIN = 2, height_coord_posterior$`beta[5]`, `*`)
  err_list <- list()
  for(i in 1:len){
    err_list[[i]] <- rnorm(1000, 0, height_coord_posterior$sigma)
  }
  err_mat <- do.call(rbind, err_list)
  return(exp(hmeans + err_mat))
}

