################################################################################
# FUNCTION:              Bimodal GEV distribution proposed in Cira EG Otiniano,
#                        Bianca S Paiva, Roberto Vila and Marcelo Bourguignon (2021)
#  dbgev                Density for the bimodal generalized extreme value distribution.
#  qbgev                Quantile function for the bimodal GEV distribution. 
#  rbgev                Random generation for the bimodal GEV distribution.
#  likbgev               maximum likelihood (ML) estimators for the parameters of a BGEV
#                        distribution
#  ebgev                Estimate the parameters of a BGEV and optionally construct 
#                        a confidence interval for the parameters.
################################################################################


library(EnvStats)


#----------------------------------------------------------------------
dbgev <- function(y, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
    # Description:
    #Compute the density for the  bimodal generalized extreme value distribution.
    #Reference: Cira EG Otiniano et al (2021). A Bimodal Model for Extremes Data.
    #1Department of Statistics, University of Brasılia, Darcy Ribeiro,
    #Brasilia, 70910-900, DF, Brazil.
    #Department of Statistics, Federal University of Rio Grande do
    # Norte, Natal, 59078-970, RN, Brazil.
    #   Parameters: y in R; mu in R; sigma > 0; xi in R ;  delta > -1;
    
    # FUNCTION:
    
    # Error treatment of input parameters
    if(sigma <= 0  || delta <= -1 )
      stop("Failed to verify condition:
             sigma <= 0  || delta <= -1")
    
    # Compute auxiliary variables:
    T      <- (y-mu)*(abs(y-mu)^delta)
    derivate_T <- (delta + 1)*(abs(y-mu)^delta)
    
    # Compute density points
    pdf    <- dgevd(T, 0, scale=sigma, shape=xi)*derivate_T
    
    # Return Value
    pdf
}
#----------------------------------------------------------------------




#----------------------------------------------------------------------
pbgev <- function(y, mu = 1, sigma = 1, xi = 0.3, delta = 2){ 
  #Distribution Function
  # FUNCTION:
  
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  
  # Compute auxiliary variables:
  Ti      <- (y-mu)*(abs(y-mu)^delta)
  # Compute 
  cdf    <- pgevd(Ti, loc=0, scale=sigma, shape=xi)
  # Return Value
  return(cdf)
}
#----------------------------------------------------------------------




#----------------------------------------------------------------------
qbgev   <- function(p, mu = 1, sigma = 1, xi = 0.3, delta = 2){
  # Description:
  #   Compute the quantile for the 
  #   Bimodal GEV distribution.
  #   Parameters: p in [0;1];  mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  
  # Compute distribution points according to their sign
  quantile <- sign(qgevd(p, 0, sigma, xi))*(abs(qgevd(p, 0, sigma, xi)))^(1/(delta + 1)) + mu
  # Return Value
  return(quantile)
}
#----------------------------------------------------------------------


#----------------------------------------------------------------------
rbgev <- function(n, mu = 1, sigma = 1, xi = 0.3, delta = 2){
  # Description:
  #   random generator for the 
  #   Bimodal GEV distribution.
  #   Parameters: p in [0;1];  mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  
  # Error treatment of input parameters
  if(sigma <= 0  || delta <= -1 )
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")
  # Compute auxiliary variables:
  U <- runif(n)
  # Compute random numbers
  rnumber <- qbgevd(U, mu, sigma, xi, delta)
  # Return Value
  return(rnumber)
}
#----------------------------------------------------------------------


#----------------------------------------------------------------------
likbgev <- function(y, theta = c(1, 1, 0.3, 2)){
  # Description:
  #  maximum likelihood (ML) function for the parameters of a BGEV distribution
  #   Parameters: y in R;  theta: vector with mu, sigma, xi and delta, respectively.
  #   mu in R; sigma > 0; xi in R ;  delta > -1;
  
  # FUNCTION:
  
  # Parameters:
  mu      <- theta[1]
  sigma   <- theta[2]
  xi      <- theta[3]
  delta   <- theta[4]
  
  # Error treatment of input parameters
  if(length(theta)!=4){
    stop("vector of parameters needs to be of length 4.")}
  if(sigma <= 0  || delta <= -1 ){
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")}
  
  # Compute auxiliary variables:
  T      <- (y-mu)*(abs(y-mu)^delta)
  derivate_T <- (delta + 1)*(abs(y-mu)^delta)
  # Compute density points
  dbgevd <- dgevd(T, mu, sigma, xi)*derivate_T
  # Log:
  logl <- sum(log(dbgevd(y,mu, sigma, xi, delta)))
  # Return negative Value for maximization
  return(-logl)
}
#----------------------------------------------------------------------



# auxiliar function
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )}

ebgev <- function( x, brep, mu = 1, sigma = 1, xi = 0.3, delta = 2,conf.level = 0.95){
  # Description:
  #   Estimate the parameters of a BGEV using MLE and optionally construct a confidence interval for
  #   the parameters.
  #   Parameters: x a numeric vector of a BGEV distribution sample;
  #   brep > 0 number of bootstrap replicas;
  #   mu in R; sigma > 0; xi in R ;  delta > -1;
  #   conf.level is the confidence level for the confidence interval.
  
  # FUNCTION:
  # Error treatment of input parameters
  if(brep<0){
    stop("Number of breplicas must be greater than 0.")}
  if(sigma <= 0  || delta <= -1 ){
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")}
  
  # Create auxiliary variables:
  muhat      = rep(NA, times = brep)
  sigmahat   = rep(NA, times = brep)
  xihat      = rep(NA, times = brep)
  deltahat   = rep(NA, times = brep)
  
  starts=c(mu, sigma, xi, delta)
  n<-length(x)
  
  # Resampling with reposition
  for(k in (1:brep)){
    i <- sample(1:n, size = n, replace = TRUE)
    Z <- x[i]
    while(show_condition(suppressWarnings(optim(par= starts, fn = likbgev, y=Z, method="BFGS")))[1]=="error"){
      Z <- x[i]
    }
    
    esti <- optim(par= starts, fn = likbgev, y=Z, method="BFGS")
    
    
    muhat[k]       <- esti$par[1]
    sigmahat[k]    <- esti$par[2]
    xihat[k]       <- esti$par[3]
    deltahat[k]    <- esti$par[4]
  }
  
  estimate_mu     <-   mean(muhat,na.rm=TRUE)
  estimate_sigma  <-   mean(sigmahat,na.rm=TRUE)
  estimate_xi     <-   mean(xihat,na.rm = TRUE)
  estimate_delta  <-   mean(deltahat,na.rm=TRUE)
  
  
  CI_infmu<- estimate_mu-qnorm((1-conf.level/2),mean=0,sd=1)*sd(muhat,na.rm = TRUE)
  CI_supmu<-estimate_mu+qnorm((1-conf.level/2),mean=0,sd=1)*sd(muhat,na.rm = TRUE)
  
  CI_infsigma<- estimate_sigma-qnorm((1-conf.level/2),mean=0,sd=1)*sd(sigmahat,na.rm = TRUE)
  CI_supsigma<-estimate_sigma+qnorm((1-conf.level/2),mean=0,sd=1)*sd(sigmahat,na.rm = TRUE)
  
  CI_infxi<- estimate_xi-qnorm((1-conf.level/2),mean=0,sd=1)*sd(xihat,na.rm = TRUE)
  CI_supxi<-estimate_xi+qnorm((1-conf.level/2),mean=0,sd=1)*sd(xihat,na.rm = TRUE)
  
  CI_infdelta<- estimate_delta-qnorm((1-conf.level/2),mean=0,sd=1)*sd(deltahat,na.rm = TRUE)
  CI_supdelta<-estimate_delta+qnorm((1-conf.level/2),mean=0,sd=1)*sd(deltahat,na.rm = TRUE)
  
  cimu<-paste("(",CI_infmu,";",CI_supmu,")")
  cisigma<-paste("(",CI_infsigma,";",CI_supsigma,")")
  cixi<-paste("(",CI_infxi,";",CI_supxi,")")
  cidelta<-paste("(",CI_infdelta,";",CI_supdelta,")")
  
  tab <- as.table(rbind(c(estimate_mu , estimate_sigma , estimate_xi ,estimate_delta ), 
                        c(cimu, cisigma, cixi,cidelta)))
  dimnames(tab) <- list(Estimative = c("Ponctual", "Conf. Interval"),
                        parameter = c("mu", "sigma", "xi","delta"))
  print(tab)
  
  
  
}






