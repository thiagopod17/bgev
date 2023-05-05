########################
# TESTS
#######################
source("BGEV_functions.R")

# dbgevd
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
          print(paste("(y,mu,sigma,xi,delta) = (",y,mu,sigma,xi,delta,") and dbgevd = ",try(dbgevd(y, mu, sigma, xi, delta)), sep = " " ))
        }
      }
    }
  }
}


dbgevd(0, mu = 0, sigma = 0.0001, xi = 0, delta = -0.9) # INVESTIGAR PQ TEM ESSE NA AQUI... ACHO 
# QUE TEM UM INTERTRAVAMENTO ENTRE Y e Delta !!

dbgevd(0, mu = 0, sigma = 0.0001, xi = 0, delta = 0)

dbgevd(-1, mu = 0, sigma = 0.0001, xi = 0, delta = -0.999)

dbgevd(Inf, mu = 0, sigma = 0.0001, xi = 1e5, delta = 0)




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

