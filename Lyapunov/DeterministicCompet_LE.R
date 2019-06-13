### FB 24/04/2019 -- Dominant lyapunov exponent code, following Dennis et al. 2001
### [Dennis, B., Desharnais, R. A., Cushing, J. M., Henson, S. M., & Costantino, R. F. (2001). Estimating chaos and complex dynamics in an insect population. Ecological Monographs, 71(2), 277-303.]
### Computes Lyapunov exponents for the deterministic two-species competition model by Sugihara et al. 2012
### Should demonstrate positive LE, due to chaos 

rm(list=ls())
graphics.off()

# The model is x[,t+1] = x[,t]*(r-A%*%x[,t]) in R notation
# Parameters
r = c(3.8,3.5)
a = matrix(c(-3.8,-0.02,-0.1,-3.5),2,2,byrow=TRUE)
#a = matrix(c(-3.8,-0.0,-0.,-3.5),2,2,byrow=TRUE)

### 1. Simulation of the model for a long time
tmax = 10000
x=matrix(0,nrow=2,ncol=tmax)

x[,1] = c(0.1,0.1)

for (t in 1:(tmax-1)){
  x[,t+1] = x[,t]*(r+a%*%x[,t])
}

# Check time series by plotting a subset
time_interval = 500:800
plot(time_interval,x[1,time_interval],type="b",col="blue")
lines(time_interval,x[2,time_interval],type="b",col="red")
#OK

### 2. Definition of the Jacobian
Jacobian <-function(x,r,a){
  J = matrix(0,2,2)
  J[1,1] = r[1] + 2*a[1,1]*x[1]+a[1,2]*x[2]
  J[1,2] = a[1,2]*x[1]
  J[2,1] = a[2,1]*x[2]
  J[2,2] = r[2] + a[2,1]*x[1] + 2*a[2,2]*x[2]
  return(J)
}

### 3. LE computation

#Initializing
s=rep(1,tmax)
J1 = Jacobian(x[,1],r,a)
s[1]=norm(J1,type="2");
S=J1/s[1];
productStemp=S;

for (t in 2:tmax){
  J = Jacobian(x[,t],r,a)
  s[t]=norm(J%*%productStemp,type="2");
  S=J/s[t];
  productStemp=S%*%productStemp;
}

LE=(1/tmax)*sum(log(s));

LE

# 0.410414 (tmax=5000)
# 0.4092075 (tmax=6000)
# 0.4112785 (tmax=7000)
# 0.4142969 (tmax=10000)


