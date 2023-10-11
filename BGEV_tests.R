########################
# TESTS
#######################
source("BGEV_functions.R")
library(fBasics)

#-------------------------------------------
# dbgev
# mu in R; sigma > 0; xi in R ;  delta > -1;
mu = runif(1,-5,5)
sigma = runif(1,0.001,10)
#xi = 0.1
xi = runif(1,0.001,10)
delta = runif(1,-1,10)
support = bgev.support(mu, sigma, xi, delta)
distCheck(fun="bgev", mu, sigma, xi, delta,  n = 2000, 
          support.lower = support[1], support.upper = support[2], subdivisions = 5000 )
fBasics::distCheck(fun="bgev", mu, sigma, xi, delta,  n = 2000)




# WHY IS THIS GIVING STRANGE VALUES ?
fun = function(x, ...) {
  dbgev(x, ...)
}

fun2 = function(x) {
  dbgev(x, mu = mu, sigma = sigma, xi = xi, delta = delta)
}

integrate(f = fun, lower =  0.1, upper = 0.5, mu = mu, sigma = sigma, xi = xi, delta = delta,
          subdivisions = 5)

integrate(f = fun2, lower =  0.1, upper =1, subdivisions = 5)

integrate(fun, lower = 0.1, upper = 0.5, subdivisions = 5000, 
          mu = mu, sigma = sigma, xi = xi, delta = delta)

curve(expr = fun2,  from = support[1], to = support[2])
curve(fun, from = support[1], to = 3,mu = mu, sigma = sigma, xi = xi, delta = delta,)

fun2(0)


distCheck = function (fun = "norm", n = 1000, robust = TRUE, subdivisions = 1500, 
                      support.lower = -Inf, support.upper = Inf,
          ...) 
{
  #  ... : parameters passed to the tested function, e.g., mean and variance. 
  cat("\nDistribution Check for:", fun, "\n ")
  CALL = match.call()
  cat("Call: ")
  cat(paste(deparse(CALL), sep = "\n", collapse = "\n"), "\n", 
      sep = "")
  dfun = match.fun(paste("d", fun, sep = ""))
  pfun = match.fun(paste("p", fun, sep = ""))
  qfun = match.fun(paste("q", fun, sep = ""))
  rfun = match.fun(paste("r", fun, sep = ""))
  xmin = qfun(p = 0.01, ...)
  xmax = qfun(p = 0.99, ...)
  NORM = integrate(dfun, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)
  
  
  cat("\n====================================
      TEST 1. Normalization Check of density 'dfun': NORM  
====================================
     ")
  cat("test if the density function integrates out to 1. Notice here 
      that it would be advisable giving better limiting lower and upper bounds 
      for integration in case the function is not defined in some region.
      RESULT OF NUMERICAL INTEGRATION: \n")
  print(NORM)
  normCheck = (abs(NORM[[1]] - 1) < 0.01)
  
  
  cat("\n====================================
      TEST 2. Quantile and CDF function [p-pfun(qfun(p))]^2 Check: 
====================================
      ")
  p = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  P = pfun(qfun(p, ...), ...)
  pP = round(rbind(p, P), 3)
  rownames(pP) = c("fixed quantiles","computed quantiles with pfun and qfun")
  print(pP)
  RMSE = sd(p - P)
  print(c("Computed RMSE of the difference",RMSE))
  rmseCheck = (abs(RMSE) < 1e-04)
  
  
  
  cat("\n====================================
      TEST 3. Sample moments comparison using rfun and dfun
====================================
      Sample moments using the RNG of the variable and comparing with 
      the value obtained from integrating the density. 
      CAREFULL FOR DISTRIBUTIONS WHICH DOES NOT HAVE WELL DEFINED FIRST AND/OR SECOND MOMENTS.
r(", n, ") Check sample moments:\n", sep = "")
  r = rfun(n = n, ...)
  if (!robust) {
    SAMPLE.MEAN = mean(r)
    SAMPLE.VAR = var(r)
  }
  else {
    robustSample = MASS::cov.mcd(r, quantile.used = floor(0.95 * 
                                                            n))
    SAMPLE.MEAN = robustSample$center
    SAMPLE.VAR = robustSample$cov[1, 1]
  }
  SAMPLE = data.frame(t(c(MEAN = SAMPLE.MEAN, VAR = SAMPLE.VAR)), 
                      row.names = "SAMPLE")
  print(signif(SAMPLE, 3))
  fun1 = function(x, ...) {
    x * dfun(x, ...)
  }
  fun2 = function(x, ...) {
    x^2 * dfun(x, ...)
  }
  MEAN = integrate(fun1, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                   stop.on.error = FALSE, ...)

  VAR = integrate(fun2, lower = support.lower, upper = support.upper, subdivisions = subdivisions, 
                  stop.on.error = FALSE, ...)
  
  EXACT = data.frame(t(c(MEAN = MEAN[[1]], VAR = VAR[[1]] - 
                           MEAN[[1]]^2)), row.names = "EXACT ")
  cat("\nCheck moments by numerical integration\n")
  print(signif(EXACT, 3))
  meanvarCheck = (abs(SAMPLE.VAR - EXACT$VAR)/EXACT$VAR < 0.1)
  
  cat("\nPrecision of moment computation using integration \n")
  cat("   X   ")
  print(MEAN)
  cat("   X^2 ")
  print(VAR)
  
  
  cat("\n====================================
      REPORT OF ALL 3 TESTS. SHOULD BE ALL TRUE TO PASS 
====================================\n")
  ans = list(normCheck = normCheck, rmseCheck = rmseCheck, 
             meanvarCheck = meanvarCheck)
  unlist(ans)
}











# LIMIT TESTS 
# y in R; mu in R; sigma > 0; xi in R ;  delta > -1;
size = 1e5
limitvalues = matrix(c(-size, size,
                       -size, size,
                       size, size,
                       -size, size,
                        0, size), 5, 2, byrow = TRUE)
rownames(limitvalues) = c("y","mu","sigma","xi","delta")
colnames(limitvalues) = c("liminf","limsup")
for(y in limitvalues[1,]){
  for(mu in limitvalues[2,]){
    for(sigma in limitvalues[3,]){
      for(xi in limitvalues[4,]){
        for(delta in limitvalues[5,]){
          print(paste("(y,mu,sigma,xi,delta) = (",y,mu,sigma,xi,delta,") and dbgev = ",try(dbgev(y, mu, sigma, xi, delta)), sep = " " ))
        }
      }
    }
  }
}


dbgev(0, mu = 0, sigma = 0.0001, xi = 0, delta = -0.9) # INVESTIGAR PQ TEM ESSE NA AQUI... ACHO 
# QUE TEM UM INTERTRAVAMENTO ENTRE Y e Delta !!

dbgev(0, mu = 0, sigma = 0.0001, xi = 0, delta = 0)

dbgev(-1, mu = 0, sigma = 0.0001, xi = 0, delta = -0.999)

dbgev(Inf, mu = 0, sigma = 0.0001, xi = 1e5, delta = 0)




# qbgevd
# y in R; mu in R; sigma > 0; xi in R ;  delta > -1;
size = 1e5
limitvalues = matrix(c(1/size, 1-1/size,
                       -size, size,
                       size, size,
                       -size, size,
                       -1 + 1/size, size), 5, 2, byrow = TRUE)
rownames(limitvalues) = c("p","mu","sigma","xi","delta")
colnames(limitvalues) = c("liminf","limsup")
for(p in limitvalues[1,]){
  for(mu in limitvalues[2,]){
    for(sigma in limitvalues[3,]){
      for(xi in limitvalues[4,]){
        for(delta in limitvalues[5,]){
          print(paste("(p,mu,sigma,xi,delta) = (",p,mu,sigma,xi,delta,") and qbgevd = ",try(qbgevd(p, mu, sigma, xi, delta)), sep = " " ))
        }
      }
    }
  }
}

