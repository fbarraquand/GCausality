### FB 24/04/2019 -- Dominant lyapunov exponent code, following Dennis et al. 2001
### [Dennis, B., Desharnais, R. A., Cushing, J. M., Henson, S. M., & Costantino, R. F. (2001). Estimating chaos and complex dynamics in an insect population. Ecological Monographs, 71(2), 277-303.]
### Computes Lyapunov exponents for our stochastic 20-species interaction model

### In that case, the Jacobian is calculated on the log-scale so the model output should be on the logscale too

rm(list=ls())
graphics.off()

### 1. Definition of the Jacobian

### Now we need it around each point
Jacobian <-function(z,alph){
  nsp = length(z)
  N_mat = matrix(rep(as.vector(exp(z)),nsp),nsp,nsp,byrow=T)
  J = diag(nsp) + alph*N_mat
  return(J)
}

### 2. Interaction matrix

### We use the one below, modifying a 10 species model
interactions10species = rbind(c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,1,1,0,0,0,0,0),
                              c(0,0,0,1,1,1,1,0,0,0),
                              c(0,0,0,0,1,1,1,0,0,0),
                              c(0,0,0,0,1,1,1,0,0,0),
                              c(0,0,0,0,0,0,0,1,1,1),
                              c(0,0,0,0,0,0,0,1,1,1),
                              c(0,0,0,0,0,0,0,1,1,1))
null_mat = matrix(0,10,10)
interaction_matrix = rbind(cbind(interactions10species,null_mat),cbind(null_mat,interactions10species))
### Adding some links between the two main compartments by adding a big connected cluster
interaction_matrix[8:13,8:13] = matrix(1,6,6) ## module filled with ones
interaction_matrix
image(interaction_matrix)
nspecies = nrow(interaction_matrix)

### 3. Simulation and computation of SLE

# How many repeats? 
nrepeats=25
LE = rep(0,nrepeats)
ModEig = rep(0,nrepeats) # Fixed point 
### High-level parameters
tmax=3000 ## so that the computation stays practical
index_time=1:tmax
# Temporal noise on the population growth rate in logarithmic scale 
sigma=0.3 # Sigma^2= 0.1

niter=1

while (niter<(nrepeats+1))
{
  
  ### Part 3a - Finding feasible parameter set
  
  print(paste("site/repeat =",niter))
  
  # Nonlinear LV discrete-time (Ricker) competition
  stability_crit = 1.5 #initializing stability criteria
  ### Use feasability (strictly positive interior fixed point) with min threshold for abundances
  
  while(stability_crit >1){
    print("Try out new parameter set for LV model")
    
    ### Simulate the A = (alpha_ij) matrix
    min_inter_pos = 0.05
    max_inter_pos = 0.10
    min_inter_neg = -0.1
    max_inter_neg = -0.2
    min_intra = - 0.8
    max_intra = -0.3
    frac_pos = 0.2
    ### To avoid introducing too weak interactions, we change the techniques
    ### There is a fraction of positive interactions, set to frac_pos, which is set to 20% of off-diagonal interactions
    ### If a Bernouilli draw decides positive, drawn from a Beta distrib for positive values and else from a Beta distrib for negative values
    ### The Beta distribution is scaled, so that the min_inter_pos max_inter_pos, min_inter_neg, max_inter_neg, min_intra,max_intra 
    ### specify instead the distribution. I could use Uniform even (would be easiest to compute the variance). 
    
    alpha=interaction_matrix
    for (ki in 1:nspecies){
      for (kj in 1:nspecies){
        noise_strength=rbeta(1,2,2) ### so that there's a min interaction strength
        ## noise_strength=runif(1,0,sigma_A) ## somewhat arbitrary but chosen to keep sign and magnitude of interaction in check
        if (ki==kj){alpha[ki,kj]= alpha[ki,kj]*(min_intra+(max_intra-min_intra)*noise_strength) } 
        else  { sign_inter = rbinom(1,1,frac_pos)
        if (sign_inter==1)
        {
          alpha[ki,kj]= alpha[ki,kj]*(min_inter_pos+(max_inter_pos-min_inter_pos)*noise_strength)
        }   
        else {alpha[ki,kj]= alpha[ki,kj]*(min_inter_neg+(max_inter_neg-min_inter_neg)*noise_strength)}
        }
      }
    }
    alpha
    
    ### Simulate the growth rates -- better to simulate them here to check the combination alpha - r is feasible
    r<-runif(nspecies,0.5,1) # Maximal growth rates // could covary with interaction strength later on. 
    # r<-runif(nspecies,2,3) # for more chaotic dynamics. 
    
    ### Equilibrium, mean values - some algebra
    N_star<-solve(alpha) %*% (-r) ### check the equilibrium is feasible, at least with those r values...
    N_star # 
    
    # Simulate the model
    
    Y=matrix(1,nrow=tmax,ncol=nspecies)
    epsilon<-matrix(0, nrow=tmax-1,ncol=nspecies,byrow=TRUE) # Noise on growth rates
    
    for (i in 1:nspecies){
      Y[1,i]=abs(rnorm(1,1,1))} #### Initial condition - most changing part across repeats
    for (t in 1:(tmax-1)){
      epsilon[t,]<-rnorm(nspecies,mean = 0, sd = sigma) # Residual noise on log-population growth rates
      Y[t+1,] = Y[t,]*exp(r+alpha %*% Y[t,] + epsilon[t,]) # Ricker dynamics // Note FB: Beverton-Holt would be another option, more stable
    }
    y=log(Y)
    #y=y-mean(y)
    
    matplot(y)
    matlines(y)
    
    #feasible = (min(y)>(-25))
    #if (feasible==TRUE){print("Found LV TS")}
    #found one right parameter set 
    if (min(N_star)>0.00001){stability_crit=0.5; print("Found LV param set")}

  }


### Part 3b - LE computation (within the loop on repeats)

#Initializing
s=rep(1,tmax)
J1 = Jacobian(y[1,],alpha)
s[1]=norm(J1,type="2");
S=J1/s[1];
productStemp=S;

for (t in 2:tmax){
  J = Jacobian(y[t,],alpha)
  s[t]=norm(J%*%productStemp,type="2");
  S=J/s[t];
  productStemp=S%*%productStemp;
}
### Store and print Lyapunov exponent
LE[niter]=(1/tmax)*sum(log(s));
print(paste("SLE = ",LE[niter]))

# computes whether the fixed point Nstar is stable. 
Jstar = Jacobian(log(N_star),alpha)
ModEig[niter]=max(Mod(eigen(Jstar)$values))

print(paste("Modulus eigenvalue (fixed point)= ",ModEig[niter]))

### Increment iteration counter
niter=niter+1

}

LE

## -0.020174461 -0.122244437 -0.028048958 -0.084889740 -0.021444842 -0.097381420 -0.089103743 -0.022303016 -0.016130163
## -0.008933913 -0.041934415 -0.019597191 -0.022439223 -0.038746568 -0.013317825 -0.029779366 -0.036169234 -0.005780739
## -0.019796340 -0.061151384 -0.018476075 -0.004844198 -0.069370259 -0.189545005 -0.121243528

## -0.03200816 -0.08300824 -0.17171034 -0.11553928 -0.08729588 -0.03383929 -0.09145159 -0.07055240 -0.17389386 -0.09508202
## -0.10106448 -0.10103577 -0.10846045 -0.05416752 -0.06308711 -0.08895230 -0.01030335 -0.01975854 -0.07279155 -0.08024503
## -0.01234832 -0.08004140 -0.07243068 -0.06091174 -0.01842431

ModEig


