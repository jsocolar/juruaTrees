
stemwise_AGB_diff <- function (D1, D2 = NULL, WD1, WD2 = NULL, H1, H2 = NULL) {
  len <- length(D1)
  param_4 <- BIOMASS::param_4[2:1001,]
  
  if(is.null(D2)){D2 <- D1}
  if(is.null(WD2)){WD2 <- WD1}
  if(is.null(H2)){H2 <- H1}
  
  D1_expanded <- replicate(1000, D1)
  D2_expanded <- replicate(1000, D2)
  WD1_expanded <- replicate(1000, WD1)
  WD2_expanded <- replicate(1000, WD2)
  H1_expanded <- replicate(1000, H1)
  H2_expanded <- replicate(1000, H2)
  
  RSE <- param_4[, "sd"]
  matRSE <- mapply(function(y) {
    rnorm(sd = y, n = len)
  }, y = RSE)
  Ealpha <- param_4[, "intercept"]
  Ebeta <- param_4[, "logagbt"]
  Comp1 <- t(log(WD1_expanded * H1_expanded * D1_expanded^2)) * Ebeta + 
    Ealpha
  Comp1 <- t(Comp1) + matRSE
  AGB_simu1 <- exp(Comp1)/1000
  
  Comp2 <- t(log(WD2_expanded * H2_expanded * D2_expanded^2)) * Ebeta + 
    Ealpha
  Comp2 <- t(Comp2) + matRSE
  AGB_simu2 <- exp(Comp2)/1000
  
  return(list(agb_1 = AGB_simu1, agb_2 = AGB_simu2, 
              stemwise_AGB_diff = AGB_simu1-AGB_simu2,
              stemwise_AGB_diff_log_prop = log((AGB_simu1-AGB_simu2)/AGB_simu1)))
}
