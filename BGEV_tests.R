#auxiliar function
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )}

ebgevd<- function( x,Rep,mu,sigma,xi,delta,conf.level = 0.95){
  # Description:
  #   Estimate the parameters of a BGEV and optionally construct a confidence interval for
  #   the parameters.
  #   Parameters: x a numeric vector of a BGEV distribution sample;
  #   Rep > 0 number of bootstrap replicas;
  #   mu in R; sigma > 0; xi in R ;  delta > -1;
  #   conf.level is the confidence level for the confidence interval.
  
  # FUNCTION:
  # Error treatment of input parameters
  if(Rep<0){
    stop("Number of Replicas must be greater than 0.")}
  if(sigma <= 0  || delta <= -1 ){
    stop("Failed to verify condition:
           sigma <= 0  || delta <= -1")}
  
  # Create auxiliary variables:
  muhat      = rep(NA, times = Rep)
  sigmahat   = rep(NA, times = Rep)
  xihat      = rep(NA, times = Rep)
  deltahat   = rep(NA, times = Rep)
  
  starts=c(mu, sigma, xi, delta)
  n<-length(x)
  
  #Resampling with reposition
  for(k in (1:Rep)){
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
################################
#Examples

n=50
Z <- rdgevd(n, mu, sigma, xi, delta)
ebgevd(Z,R=10,mu=1,sigma=10,xi=0.3,delta=2,conf.level = 0.95)
#                parameter
#Estimative       mu                                     
#Ponctual       1.02279116583337                       
#Conf. Interval ( 1.01168448475556 ; 1.03389784691119 )
#parameter
#Estimative       sigma                                  
#Ponctual       16.9963296503858                       
#Conf. Interval ( 16.3867158696069 ; 17.6059434311648 )
#parameter
#Estimative       xi                                       
#Ponctual       0.139581984753769                        
#Conf. Interval ( 0.124364259639566 ; 0.154799709867972 )
#parameter
#Estimative       delta                                  
#Ponctual       2.49602823332784                       
#Conf. Interval ( 2.45331666498226 ; 2.53873980167342 )

