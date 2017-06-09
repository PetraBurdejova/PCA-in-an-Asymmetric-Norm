library(mvtnorm)
#library(expectreg)
library(sn)
rho_x    = 0.9
mean_x = c(0,0)

sigma_x  = matrix(c(1,rho_x,rho_x,1), nrow = 2, ncol = 2)

n=300

x      = rmvnorm(n=n, mean = mean_x, sigma = sigma_x)

xi <- c(0,0)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,-10)
#
x <- rmsn(500, xi, Omega, alpha)

#xx <- seq(-3,3,length=15)
#pdf <- dmsn(cbind(xx, 2*xx-1), xi, Omega, alpha)
#cdf <- pmsn(cbind(xx, 2*xx-1), xi, Omega, alpha)

plot(x)

e1 <- 2*eigen(cov(x))$vectors[,1]
e2 <- 2*eigen(cov(x))$vectors[,2]

segments(0,0,e1[1],e1[2],lwd=2)
segments(0,0,e2[1],e2[2],lwd=2)


e1 <- -2*eigen(cov(x))$vectors[,1]
e2 <- -2*eigen(cov(x))$vectors[,2]

segments(0,0,e1[1],e1[2],lwd=2)
segments(0,0,e2[1],e2[2],lwd=2)

expectilea<- pec.k(t(x), alpha=0, nk=2,reset.tol=10)
expectilea[[5]]

e1_pec <- 2*expectilea[[2]][,1]
segments(0,0,e1_pec[1],(e1_pec[2]), col="green",lwd=2)

e2_pec <- 2*expectilea[[2]][,2]
segments(0,0,e2_pec[1],(e2_pec[2]), col="green",lwd=2)

expectilea<- pec.k(t(x), alpha=0.4, nk=2,reset.tol=10)
expectilea[[5]]

e1_pec <- 2* expectilea[[2]][,1]
segments(0,0,e1_pec[1],(e1_pec[2]), col="red",lwd=2)

e2_pec <- 2*expectilea[[2]][,2]
segments(0,0,e2_pec[1],(e2_pec[2]), col="red",lwd=2)


expectilea<- pec.k(t(x), alpha=-0.4, nk=2,reset.tol=10)
expectilea[[5]]

e1_pec <- 2*expectilea[[2]][,1]
segments(0,0,e1_pec[1],(e1_pec[2]), col="blue",lwd=2)

e2_pec <- 2*expectilea[[2]][,2]
segments(0,0,e2_pec[1],(e2_pec[2]), col="blue",lwd=2)



###################
#cond. distr.
mean_y.x <- function(x,mean_joint, sigma_joint){
  output <- mean_joint[2] 
  S21 <- sigma_joint[2,1]
  S11_inv <- 1/sigma_joint[1,1]
  output <- output +  (S21*S11_inv)*(x-mean_joint[1])
    return(output)
}
  
var_y.x <- function(sigma_joint){
  S11_inv <- 1/sigma_joint[1,1]
  S12 <- sigma_joint[1,2]
  S21 <- sigma_joint[2,1]
  output <- sigma_joint[2,2] - (S21*S11_inv*S12)
  return(output)
}
###################

w_alpha <- function(alpha_lev, Ymean, Yvariance){
  q_alpha <- qnorm(alpha_lev, mean =Ymean, sd = sqrt(Yvariance))
  q_alpha_std <- (q_alpha - Ymean) / sqrt(Yvariance)
  G_q <- - sqrt(Yvariance)* dnorm(q_alpha_std) + Ymean*pnorm(q_alpha_std)
  output <- (- alpha_lev*q_alpha + G_q) / (-Ymean + 2*G_q + (1-2*alpha_lev)*q_alpha)
  return(output)
}


rho_xy    = 0.9
mean_xy = c(0,0)
sigma_xy  = matrix(c(3,rho_xy,rho_xy,1), nrow = 2, ncol = 2)

n=1000
data      = rmvnorm(n=n, mean = mean_xy, sigma = sigma_xy)
plot(data)


#  ALPHA=0.5  TAU=0.5
alpha=0.5
qs <- rep (0, dim(data)[1])
ws <- rep (0, dim(data)[1])
for (i in 1:dim(data)[1]){
  x <- data[i,1]
  mu_y.x<- mean_y.x(x,mean_xy,sigma_xy)
  variance_y.x <- var_y.x(sigma_xy) 
   w <- w_alpha(alpha_lev=alpha, mu_y.x,variance_y.x )
   q_alpha <- qnorm(alpha, mean =mu_y.x, sd = sqrt(variance_y.x))
   qs[i] <- q_alpha
   ws[i] <- w
}
cat("alpha is ", alpha, "tau is", ws[i])
points(data[,1], qs, type="l", col="black")


#  ALPHA=0.806  TAU=0.9004
alpha=0.806
qs <- rep (0, dim(data)[1])
ws <- rep (0, dim(data)[1])
for (i in 1:dim(data)[1]){
  x <- data[i,1]
  mu_y.x<- mean_y.x(x,mean_xy,sigma_xy)
  variance_y.x <- var_y.x(sigma_xy) 
  w <- w_alpha(alpha_lev=alpha, mu_y.x,variance_y.x )
  q_alpha <- qnorm(alpha, mean =mu_y.x, sd = sqrt(variance_y.x))
  qs[i] <- q_alpha
  ws[i] <- w
}
cat("alpha is ", alpha, "tau is", ws[i])
points(data[,1], qs, type="l", col="red")

#  ALPHA=0.195  TAU=0.1005
alpha=0.195
qs <- rep (0, dim(data)[1])
ws <- rep (0, dim(data)[1])
for (i in 1:dim(data)[1]){
  x <- data[i,1]
  mu_y.x<- mean_y.x(x,mean_xy,sigma_xy)
  variance_y.x <- var_y.x(sigma_xy) 
  w <- w_alpha(alpha_lev=alpha, mu_y.x,variance_y.x )
  q_alpha <- qnorm(alpha, mean =mu_y.x, sd = sqrt(variance_y.x))
  qs[i] <- q_alpha
  ws[i] <- w
}
cat("alpha is ", alpha, "tau is", ws[i])
points(data[,1], qs, type="l", col="blue")

expectilea<- pec.k(t(data), alpha=0, nk=2,reset.tol=10)
expectilea[[5]]


e1_pec <- 2*expectilea[[2]][,1]
segments(0,0,e1_pec[1],(e1_pec[2]), col="green",lwd=3)

e2_pec <- 2*expectilea[[2]][,2]
segments(0,0,e2_pec[1],(e2_pec[2]), col="green",lwd=3)






  
  
  
