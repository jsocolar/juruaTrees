test_sites <- getBioclimParam(cbind(rep(-70,30),rep(-8,30)))
posterior_E <- as.matrix(height_coord_posterior[,c("beta[1]", "beta[2]", "beta[3]")]) %*% t(as.matrix(test_sites))


predict_heights <- function(D, model=NULL, coord=NULL, plot=NULL){
  len = length(D)
  if(is.null(model)+is.null(coord) != 1){stop("must specify exactly one of model and coord")}
  if(!is.null(model)){
    return(BIOMASS:::predictHeight(D, model, err = T, plot))
  }else{
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
}

