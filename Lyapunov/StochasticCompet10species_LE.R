### FB 24/04/2019 -- Dominant lyapunov exponent code, following Dennis et al. 2001
### [Dennis, B., Desharnais, R. A., Cushing, J. M., Henson, S. M., & Costantino, R. F. (2001). Estimating chaos and complex dynamics in an insect population. Ecological Monographs, 71(2), 277-303.]
### Computes Lyapunov exponents for our stochastic 10-species (mostly) competition model

### In that case, the Jacobian is calculated on the log-scale so the model output should be on the logscale too

rm(list=ls())
graphics.off()

### Same exact network structure used 
interaction_matrix = rbind(c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,1,0,0,0,0,0,0),
                           c(0,0,0,1,1,1,1,0,0,0),
                           c(0,0,0,0,1,1,1,0,0,0),
                           c(0,0,0,0,1,1,1,0,0,0),
                           c(0,0,0,0,0,0,0,1,1,1),
                           c(0,0,0,0,0,0,0,1,1,1),
                           c(0,0,0,0,0,0,0,1,1,1))
nspecies=nrow(interaction_matrix)

### We use a reference parameter set
r = c(4,3.1,0.12,3,3,3.1,0.8,4,3.1,0.12)
alpha = rbind (c(-4,-2,-0.4,0,0,0,0,0,0,0),
               c(-0.31,-3.1,-0.93,0,0,0,0,0,0,0),
               c(0.636,0.636,-2.12,0,0,0,0,0,0,0),
               c(-0.111,-0.111,0.131,-3.8,0,0,0,0,0,0),
               c(0,0,0,0.5,-2,-2,-0.4,0,0,0),
               c(0,0,0,0,-0.31,-3.1,-0.93,0,0,0),
               c(0,0,0,0,0.636,-0.636,-2.12,0,0,0), #instead +0.636,4.12 for deterministic model
               c(0,0,0,0,0,0,0,-4,-2,-0.4),
               c(0,0,0,0,0,0,0,-0.31,-3.1,-0.93),
               c(0,0,0,0,0,0,0,0.636,0.636,-2.12)
)
### Check correspondance 
((abs(alpha)>0)==interaction_matrix)

sigma=0.3 # Sigma^2= 0.1
#sigma = 0.1 

### 1. Simulation of the model for a long time
tmax = 3000
x=matrix(0,nrow=nspecies,ncol=tmax)
x[,1]=abs(rnorm(nspecies,1,1)) #### Initial condition - most changing part across repeats

for (t in 1:(tmax-1)){
  x[,t+1] = x[,t] + r + alpha %*% exp(x[,t]) + rnorm(nspecies,0,sigma)# Ricker compet model on the log-scale
}

# Check time series by plotting a subset
time_interval = 500:800
matplot(as.matrix(time_interval),t(x[,time_interval]),type="b",ylim=c(-5,5))
#OK

### 2. Definition of the Jacobian

### Equilibrium, mean values - some algebra
N_star<-solve(alpha) %*% (-r) ### check the equilibrium is feasible
N_star # good -- should always be stable in theory
N_star_mat = matrix(rep(as.vector(N_star),nspecies),nspecies,nspecies,byrow=T)
Jacobian_eq = diag(nspecies) + alpha*N_star_mat ## Jacobian at the equilibrium values

### Now we need it around each point
Jacobian <-function(z,alph){
  nspecies = length(z)
  N_mat = matrix(rep(as.vector(exp(z)),nspecies),nspecies,nspecies,byrow=T)
  J = diag(nspecies) + alph*N_mat
  return(J)
}

### 3. LE computation

#Initializing
s=rep(1,tmax)
J1 = Jacobian(x[,1],alpha)
s[1]=norm(J1,type="2");
S=J1/s[1];
productStemp=S;

for (t in 2:tmax){
  J = Jacobian(x[,t],alpha)
  s[t]=norm(J%*%productStemp,type="2");
  S=J/s[t];
  productStemp=S%*%productStemp;
}

LE=(1/tmax)*sum(log(s));

LE

# 0.3289029 (tmax=3000)
# 0.328281  (tmax=3000)
# 0.3270377 (tmax=3000)

### With sigma=0.1
# 0.3402373 (tmax=1000)
# 0.3540836 (tmax=2000)
# 0.3279887 (tmax=3000)
# 0.3143527 (tmax=3000)

### Note: we might have some noise-induced chaos (sensu Ellner and Turchin) 