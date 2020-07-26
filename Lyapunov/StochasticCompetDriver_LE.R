### FB 24/04/2019 -- Dominant lyapunov exponent code, following Dennis et al. 2001
### [Dennis, B., Desharnais, R. A., Cushing, J. M., Henson, S. M., & Costantino, R. F. (2001). Estimating chaos and complex dynamics in an insect population. Ecological Monographs, 71(2), 277-303.]
### Computes Lyapunov exponents for our stochastic two-species competition model
### CP 13/06/2019 Add the driver
### In that case, the Jacobian is calculated on the log-scale so the model output should be on the logscale too

rm(list=ls())
graphics.off()

# Parameters
r = c(3,2.1) #max growth rate  
#a = matrix(c(-4,-2,-0.31,-3.1),2,2,byrow=TRUE) #interaction parameters
a = matrix(c(-4,0,0,-3.1),2,2,byrow=TRUE) #interaction parameters
sigma = 0.1 #Note: we have sigma=0.1 and not sigma^2 which is 0.01

### 1. Simulation of the model for a long time
tmax = 10000
x=matrix(0,nrow=2,ncol=tmax)

seasonality<-2*sin(2*pi*(1:tmax)/24)
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
y1<-seasonality+y1noise
y1=rbind(y1,y1)

x[,1] = c(0,0.1)

for (t in 1:(tmax-1)){
  x[,t+1] = x[,t] + r + a %*% exp(x[,t]) + 0.5*y1[,t]+rnorm(2,0,sigma)# Ricker compet model on the log-scale
}

# Check time series by plotting a subset
time_interval = 500:800
plot(time_interval,x[1,time_interval],type="b",col="blue")
lines(time_interval,x[2,time_interval],type="b",col="red")
#OK


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
  S=J/s[t];
  productStemp=S%*%productStemp;
}

LE=(1/tmax)*sum(log(s));

LE

### Note: with more noise we go closer to noise-induced chaos (sensu Ellner and Turchin) without fully reaching it. 
### FB: Above sentence was for the previous model without the added forcing, but it does looks like we do reach noise-induced chaos here though, since the LE is positive. 