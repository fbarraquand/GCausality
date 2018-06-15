####################################################################################
### FB 12/06/2018 -- Generating replicated datasets for 10 species #################
### Both Lotka-Volterra discrete-time (Ricker) competition and VAR(1) (Gompertz) ###
####################################################################################

rm(list = ls())
set.seed(5)

### Mainly competition dataset

### Do we use the same exact network structure? Sounds OK
### We use the one below

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
nspecies = nrow(interaction_matrix)

### Random draw of parameter values using Gaussian variables with an interaction strength mean and variance

# How many repeats? 
nrepeats=25

### High-level parameters
tmax=500 ## Number of time units. 500 timesteps total for the simulation, and we select the 300 last time steps for model fitting
index_time=1:tmax
# Variability on the matrix elements
sigma_A = sqrt(0.0025) # sigma_A^2 = 0.05
# Temporal noise on the population growth rate in logarithmic scale 
sigma=0.3 # Sigma^2= 0.1

niter=1

while (niter<nrepeats)
{

print(paste("site/repeat =",niter))
  
# Nonlinear LV discrete-time (Ricker) competition
stability_crit = 1.5 #initializing stability criteria
### Use feasability (strictly positive interior fixed point) with min threshold for abundances

while(stability_crit >1){
  print("Try out new parameter set for LV model")
  
  ### Simulate the A = (alpha_ij) matrix
  interspe_strength = 0.2
  intraspe_strength = 0.5
  
  alpha=interaction_matrix
  for (ki in 1:nspecies){
    for (kj in 1:nspecies){
      noise_strength=runif(1,0,sigma_A) ## somewhat arbitrary but chosen to keep sign and magnitude of interaction in check
      if (ki==kj){alpha[ki,kj]= - alpha[ki,kj]*intraspe_strength+noise_strength} else {alpha[ki,kj]= - alpha[ki,kj]*interspe_strength+noise_strength}
    }
  }
  alpha
  

 
  ### Simulate the growth rates -- better to simulate them here to check the combination alpha - r is feasible
  r<-runif(nspecies,0.2,0.5) # Maximal growth rates // could covary with interaction strength later on. 
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
  
# save data to file
DataTemp_wTime_abs_LV = data.frame(niter*rep(1,tmax),index_time,Y)
DataTemp_wTime_abs_LV  = format(DataTemp_wTime_abs_LV,digits=3,scientific=T)
### Write this in total dataframes
if (niter==1){ 	# Create data structures
  Data_wTime_abs_LV=DataTemp_wTime_abs_LV
} else {
  Data_wTime_abs_LV=rbind(Data_wTime_abs_LV,DataTemp_wTime_abs_LV)}
# save quantitative interaction matrix
write.csv(alpha,file=paste("../data/random_param_set/alpha/alphaMatrix",niter,".csv",sep=""))
write.csv(r,file=paste("../data/random_param_set/r/rVector",niter,".csv",sep=""))


### Now switch model 

# M/VAR(1) (Gompertz) competition model
spectral_radius = 1.5 #initializing
  while(spectral_radius>1){
    print("Try out new initial conditions VAR")
    
    # Define the interaction matrix
    B=interaction_matrix*rnorm(nspecies*nspecies,-0.05,sigma_A)+0.8*diag(nspecies)
    # Your typical MAR(1) matrix - is it stable
    spectral_radius = max(abs(eigen(B)$values)) #stable
    
    # Simulate the VAR(1) model
    ysim=matrix(0,nrow = tmax,ncol=nspecies)
    ysim[1,] = rnorm(nspecies,0.5,0.1)
    for (time_index in 1:(tmax-1)){
      ysim[time_index+1,] = B %*% ysim[time_index,] + rnorm(nspecies,0,sigma)
    }
    matplot(ysim)
    matlines(ysim) #check
    
    # Say we found a nice param set
    if (spectral_radius>1){print("Found VAR(1) param set")}
  }
# save data to file
DataTemp_wTime_abs_VAR = data.frame(niter*rep(1,tmax),index_time,Y)
DataTemp_wTime_abs_VAR  = format(DataTemp_wTime_abs_VAR,digits=3,scientific=T)
### Write this in total dataframes
if (niter==1){ 	# Create data structures
  Data_wTime_abs_VAR=DataTemp_wTime_abs_VAR
} else {
  Data_wTime_abs_VAR=rbind(Data_wTime_abs_VAR,DataTemp_wTime_abs_VAR)}
# save quantitative interaction matrix
write.csv(alpha,file=paste("../data/random_param_set/B/BMatrix",niter,".csv",sep=""))

### Increment iteration counter

niter=niter+1
}

### --------------------- Write data to file ------------------------- ### 

names(Data_wTime_abs_LV)=c("Site","Time_index","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10")
names(Data_wTime_abs_VAR)=c("Site","Time_index","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10")

#options(scipen=-500)
#Data_wTime_abs_LV = format(Data_wTime_abs_LV,digits=4)
#Data_wTime_abs_VAR = format(Data_wTime_abs_VAR,digits=4)
#sapply(Data_wTime_abs_LV, class)

write.csv(Data_wTime_abs_LV,file="../data/random_param_set/Data_wTime_abs_LV.csv")
write.csv(Data_wTime_abs_VAR,file="../data/random_param_set/Data_wTime_abs_VAR.csv")
