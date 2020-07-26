### FB 24/04/2019 -- Dominant lyapunov exponent code, following Dennis et al. 2001
### [Dennis, B., Desharnais, R. A., Cushing, J. M., Henson, S. M., & Costantino, R. F. (2001). Estimating chaos and complex dynamics in an insect population. Ecological Monographs, 71(2), 277-303.]
### Computes Lyapunov exponents for our stochastic two-species competition model

### In that case, the Jacobian is calculated on the log-scale so the model output should be on the logscale too

rm(list=ls())
graphics.off()

# Parameters
r = c(3,2.1) #max growth rate  
a = matrix(c(-4,-2,-0.31,-3.1),2,2,byrow=TRUE) #interaction parameters
sigma = 0.1 #Note: we have sigma=0.1 and not sigma^2 which is 0.01

### 1. Simulation of the model for a long time
tmax = 10000
x=matrix(0,nrow=2,ncol=tmax)
xd=matrix(0,nrow=2,ncol=tmax) # deterministic sim. 


x[,1] = c(0,0.1)
xd[,1] = c(0,0.1)


for (t in 1:(tmax-1)){
  x[,t+1] = x[,t] + r + a %*% exp(x[,t]) + rnorm(2,0,sigma)# Ricker compet model on the log-scale
  xd[,t+1] = xd[,t] + r + a %*% exp(xd[,t])  # Ricker compet model on the log-scale
  #x[,t+1] = x[,t] + r + a %*% exp(x[,t]) + rnorm(2,0,sqrt(0.1))# More noise
}

# Check time series by plotting a subset
time_interval = 500:800
plot(time_interval,x[1,time_interval],type="b",col="blue")
lines(time_interval,x[2,time_interval],type="b",col="red")
#OK

# Plotting deterministic version
plot(time_interval,xd[1,time_interval],type="b",col="blue",ylim=c(-3,1))
lines(time_interval,xd[2,time_interval],type="b",col="red")
# Zoom
time_interval = 780:800
plot(time_interval,xd[1,time_interval],type="b",col="blue",ylim=c(-3,1))
lines(time_interval,xd[2,time_interval],type="b",col="red")
# The attractor is a two cycle -> with a negative deterministic LE

### 2. Definition of the Jacobian
Jacobian <-function(x,a){
  J = matrix(0,2,2)
  J[1,1] = 1 + a[1,1]*exp(x[1])
  J[1,2] = a[1,2]*exp(x[2])
  J[2,1] = a[2,1]*exp(x[1])
  J[2,2] = 1 + a[2,2]*exp(x[2])
  return(J)
}


### 3. LE computation

#Initializing
s=rep(1,tmax)
J1 = Jacobian(x[,1],a)
s[1]=norm(J1,type="2");
S=J1/s[1];
productStemp=S;

for (t in 2:tmax){
  J = Jacobian(x[,t],a)
  s[t]=norm(J%*%productStemp,type="2");
  #s[t]=norm(J%*%productStemp,type="f"); #other possibility. 
  S=J/s[t];
  productStemp=S%*%productStemp;
}

LE=(1/tmax)*sum(log(s));

LE

# -0.1950553 (tmax = 10000)
# -0.1896428 (tmax = 10000)
# -0.1927537 (tmax = 10000)

### Note: with more noise we go closer to noise-induced chaos (sensu Ellner and Turchin) without fully reaching it. 
