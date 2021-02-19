test_sites <- getBioclimParam(cbind(-70,-8))
posterior_E <- as.matrix(height_coord_posterior[,c("beta[1]", "beta[2]", "beta[3]")]) %*% t(as.matrix(test_sites))


predict_heights <- function(D, model=NULL, coord=NULL, plot=NULL){
  len = length(D)
  if(is.null(model)+is.null(coord) != 1){stop("must specify exactly one of model and coord")}
  if(!is.null(model)){
    return(predictHeight(D, model, err = T, plot))
  }else{
    posterior_E <- as.matrix(height_coord_posterior[,c("beta[1]", "beta[2]", "beta[3]")]) %*% t(as.matrix(getBioclimParam(coord)))
    lD <- replicate(1000,log(D))
    hmeans <- t(replicate(len, height_coord_posterior$intercept)) + t(replicate(len, as.vector(posterior_E))) + 
      sweep(lD, MARGIN=2, height_coord_posterior$`beta[4]`, `*`) +
      sweep(lD^2, MARGIN = 2, height_coord_posterior$`beta[5]`, `*`)
    hsds <- as.matrix(rnorm(1000*len, 0, rep(height_coord_posterior$sigma, each = len)) # make sure to fill this thing in the right order
    return()
  }
}

