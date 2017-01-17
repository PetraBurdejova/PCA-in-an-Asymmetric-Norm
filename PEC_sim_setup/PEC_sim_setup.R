# ---------------------------------------------------------------------
# Paper:       	 N.M. Tran, P. Burdejova, M. Osipenko, and W.K. Härdle
#				"Principal Component Analysis in an Asymmetric Norm"
# ---------------------------------------------------------------------
# Quantlet:    	 PEC_sim_setup
# ---------------------------------------------------------------------
# Description: 	 Example of curves used for simulation.
# ---------------------------------------------------------------------
# Author:      	 Petra Burdejova, Maria Osipenko
# ---------------------------------------------------------------------

#------------------parameters
n= 20
p= 100
t=seq(0,1,length=p)
sigma_e = 0.5 
sigma_a1 = 16
sigma_a2 = 9

error_type = 1 #or 3 or 5

#----------------mean and component functions
m<-function(t){(mu=1+t+exp(-(t-0.6)^2/0.05))
		return(mu)
}
ff<-function(t){f1=sqrt(2)*sin(2*pi*t)
		f2=sqrt(2)*cos(2*pi*t)
		a1=rnorm(length(t), sigma_a1)	
		a2=rnorm(length(t),0, sigma_a2)
		return(cbind(a1*f1,a2*f2))
}

yi <- function(t,e_type){
		if (e_type==1) {e <- rnorm(length(t),0,sigma_e)}
		if (e_type==3) {e <- rt(length(t),1)}
		if (e_type==5) {e <- log(rnorm(length(t),0,sigma_e))}
		y_i <- m(t)+apply(ff(t),1,sum)
		return(y_i)
}

#----------------generate the curves
Y_mtx <- matrix(nrow=n, ncol=p)
for (i in 1:n){
	Y_mtx[i,] <- yi(t,error_type) 
}

#----------------plot the curves
matplot(t(Y_mtx),type="l", col="black", ylab="")





