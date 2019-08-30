####################################################################################
### FB 12/06/2018 -- Generating replicated datasets for 10 species #################
### Both Lotka-Volterra discrete-time (Ricker) competition and VAR(1) (Gompertz) ###
####################################################################################

rm(list = ls())
set.seed(5)

##options(scipen=-500)

### Mainly competition dataset

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

### Equilibrium, mean values - some algebra
N_star<-solve(alpha) %*% (-r) ### check the equilibrium is feasible
N_star # good -- should always be stable in theory
N_star_mat = matrix(rep(as.vector(N_star),nspecies),nspecies,nspecies,byrow=T)
Jacobian = diag(nspecies) + alpha*N_star_mat

spectral_radius = max(abs(eigen(Jacobian)$values)) 
spectral_radius
# add just a tad more stability
B = Jacobian 
diag(B) = diag(B)/3 + 0.2
B ## More reasonable diagonal values!!
spectral_radius = max(abs(eigen(B)$values)) 
spectral_radius # now barely stable, perfect

### High-level parameters
tmax=500 ## Number of time units. 500 timesteps total for the simulation, and we select the 300 last time steps for model fitting
index_time=1:tmax
# How many repeats? 
nrepeats=25

niter=1 ### initializing repeats
while (niter<nrepeats+1)
{
  
print(paste("site/repeat =",niter))
# Nonlinear LV discrete-time (Ricker) competition
feasible = FALSE #initializing
### use feasability (strictly positive interior fixed point) with min threshold for abundances?? 
while(feasible == FALSE){
  print("Try out new initial conditions LV")
  
  #simulate the model
  sigma=0.3 # Sigma^2= 0.1
  Y=matrix(1,nrow=tmax,ncol=nspecies)
  
  ### Deterministic dynamics ------------------------------------------------------------------------------------------------### 
  #for (i in 1:nspecies){
  #  Y[1,i]=abs(rnorm(1,1,1))} #### Initial condition - most changing part across repeats
  #for (t in 1:(tmax-1)){
  #Y[t+1,] = Y[t,]*exp(r+alpha %*% Y[t,]) # Ricker dynamics // Note FB: Beverton-Holt would be another option, more stable
  #}
  ### Deterministic dynamics ------------------------------------------------------------------------------------------------### 
  epsilon<-matrix(0, nrow=tmax-1,ncol=nspecies,byrow=TRUE) # Noise on growth rates
  
  for (i in 1:nspecies){
    Y[1,i]=abs(rnorm(1,1,1))} #### Initial condition - most changing part across repeats
  for (t in 1:(tmax-1)){
    epsilon[t,]<-rnorm(nspecies,mean = 0, sd = sigma) # Residual noise on log-population growth rates
    Y[t+1,] = Y[t,]*exp(r+alpha %*% Y[t,] + epsilon[t,]) # Ricker dynamics
  }
  y=log(Y)
  #y=y-mean(y)
  feasible = (min(y)>(-20))&(!is.nan(min(y)))
  if (feasible==TRUE){print("Found LV TS")}
}#found one right parameter set 

matplot(y)
matlines(y)

# save data 
DataTemp_wTime_abs_LV = data.frame(niter*rep(1,tmax),index_time,y)
DataTemp_wTime_abs_LV  = format(DataTemp_wTime_abs_LV,digits=3,scientific=T)

### Write this in total dataframes
if (niter==1){ 	# Create data structures
  Data_wTime_abs_LV=DataTemp_wTime_abs_LV
} else {
  Data_wTime_abs_LV=rbind(Data_wTime_abs_LV,DataTemp_wTime_abs_LV)}

# M/VAR(1) (Gompertz) competition model
feasible = FALSE #initializing
  while(feasible == FALSE){ 
    print("Try out new initial conditions VAR")
    
    #simulate the model  
    ysim=matrix(0,nrow = tmax,ncol=nspecies)
    ysim[1,] = rnorm(nspecies,0.5,0.1)
    for (time_index in 1:(tmax-1)){
      ysim[time_index+1,] = B %*% ysim[time_index,] + rnorm(nspecies,0,sigma)
    }
    matplot(ysim)
    matlines(ysim) #check
    feasible = (min(ysim)>(-25)) ### add smthg for NAN
    if (feasible==TRUE){print("Found VAR TS")}
    
  }

# save data to file
DataTemp_wTime_abs_VAR = data.frame(niter*rep(1,tmax),index_time,ysim)
DataTemp_wTime_abs_VAR  = format(DataTemp_wTime_abs_VAR,digits=3,scientific=T)
### Write this in total dataframes
if (niter==1){ 	# Create data structures
  Data_wTime_abs_VAR=DataTemp_wTime_abs_VAR
} else {
  Data_wTime_abs_VAR=rbind(Data_wTime_abs_VAR,DataTemp_wTime_abs_VAR)}

### Increment iteration counter
niter=niter+1
}

### --------------------- Write data to file ------------------------- ### 
#Data_wTime_abs_LV = format(Data_wTime_abs_LV,digits=6)
#Data_wTime_abs_VAR = format(Data_wTime_abs_VAR,digits=6)
names(Data_wTime_abs_LV)=c("Site","Time_index","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10")
names(Data_wTime_abs_VAR)=c("Site","Time_index","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10")
write.csv(Data_wTime_abs_LV,file="../data/ref_param_set/Data_wTime_abs_LV.csv")
write.csv(Data_wTime_abs_VAR,file="../data/ref_param_set/Data_wTime_abs_VAR.csv")
