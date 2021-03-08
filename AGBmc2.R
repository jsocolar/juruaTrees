
stemwise_AGB_diff <- function (D, WD1, WD2 = NULL, H1, H2 = NULL) {
  # D: a column of diameters
  # WD1: a column of gold-standard wood density values
  # WD2: an optional second column of wood density values
  # H1: a 1000-column matrix of gold-standard height values, each column 
  #     representing a different posterior iteration
  # H2: an optional second matrix of height values
  len <- length(D)
  param_4 <- BIOMASS::param_4[2:1001,]
  
  if(is.null(H2)){H2 <- H1}
  
  D_expanded <- replicate(1000, D)
  WD1_expanded <- replicate(1000, WD1)
  if(is.null(WD2)){
    WD2_expanded <- WD1_expanded
  }else if(is.vector(WD2)){
    WD2_expanded <- replicate(1000, WD2)
  }else if(is.matrix(WD2 & ncol(WD2)==1000)){
    WD2_expanded <- WD2
  }else{stop('improper size or type for WD2')}
  
  
  RSE <- param_4[, "sd"]
  matRSE <- mapply(function(y) {
    rnorm(sd = y, n = len)
  }, y = RSE)
  Ealpha <- param_4[, "intercept"]
  Ebeta <- param_4[, "logagbt"]
  Comp1 <- t(log(WD1_expanded * H1 * D_expanded^2)) * Ebeta + 
    Ealpha
  Comp1 <- t(Comp1) + matRSE
  AGB_simu1 <- exp(Comp1)/1000
  
  Comp2 <- t(log(WD2_expanded * H2 * D_expanded^2)) * Ebeta + 
    Ealpha
  Comp2 <- t(Comp2) + matRSE
  AGB_simu2 <- exp(Comp2)/1000
  
  return(list(agb_1 = AGB_simu1, agb_2 = AGB_simu2, 
              stemwise_AGB_diff = AGB_simu1-AGB_simu2,
              stemwise_AGB_logratio = log(AGB_simu2/AGB_simu1)))
}


AGBmc2 <- function(D, WD1, WD2 = NULL, H1, errH1 = NULL, H2 = NULL, errH2 = NULL, coord2 = NULL){
  myrtruncnorm <- function(n, lower = -1, upper = 1, mean = 0, 
                           sd = 1) {
    qnorm(runif(n, pnorm(lower, mean = mean, sd = sd), pnorm(upper, 
                                                             mean = mean, sd = sd)), mean = mean, sd = sd)
  }
    
  if(is.vector(H1)){
    if(!is.null(errH1)){
      upper <- max(H1, na.rm = T) + 15
      h1l <- list()
      for(i in 1:1000){
        h1l[[i]] <- myrtruncnorm(length(D), mean = H1, sd = errH1, lower=1.3, upper=upper)
      }
      H_1 <- do.call(cbind, h1l)
    }else{
      H_1 <- replicate(1000, H1)
    }
  }else if(is.matrix(H1) & ncol(H1) == 1000){
    H_1 <- H1
  }else{stop("problem with H1")}

  if(is.vector(H2)){
    if(!is.null(errH2)){
      upper <- max(H2, na.rm = T) + 15
      h2l <- list()
      for(i in 1:1000){
        h2l[[i]] <- myrtruncnorm(length(D), H2, errH2, lower=1.3, upper=upper)
      }
      H_2 <- do.call(cbind, h2l)
    }else{
      H_2 <- replicate(1000, H2)
    }
  }else if(is.null(H2)){
    if(!is.null(coord2)){
      H_2 <- predict_heights(D = D, coord=coord2)
    }else{
      H_2 <- H_1
    }
  }else if(is.matrix(H2) & ncol(H2) == 1000){
    H_2 <- H2
  }
  else{stop("problem with H2")}
  
  return(stemwise_AGB_diff(D = D, WD1 = WD1, WD2 = WD2, H1 = H_1, H2 = H_2))
}
