########################################################################################################################
########### FBarraquand 26/06/2018 - CCM analysis on nonlinear community dynamics, mimicking GC analysis (with pvals) ##
########################################################################################################################

#library(devtools)
#devtools::install_github("ha0ye/rEDM")
rm(list=ls())
library("vars")
library("rEDM")


### NB This code does not compute the embedding dimension (ED) which is set to 3
### On the minus side, recomputing the ED could help capture the nonlinearies 
### On the plus side, we avoid artefacts due to a too large ED when we do 500 simulations. 

########################################################################################################################
########## Sugihara two species deterministic competition model -- interactions
########################################################################################################################

tmax=300
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=0.1
Y[1]=0.2
for (t in 1:(tmax-1))
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}
X
Y
x=log(X)
y=log(Y)
# centering
x=x-mean(x)
y=y-mean(y)

z=cbind(x,y)
vec_col=c("black","red")
#pdf(file="CompetitionDeterministicNL2Species.pdf",width=12,height=6)
par(cex=1.5)
minz=min(z[151:201,])
maxz=max(z[151:201,])
for (i in 1:2){
  if (i==1){plot(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz),xlab="Time",ylab="ln(Density)")}
  lines(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz))
}
#dev.off()


### Do classic CCM analysis in the vein of sardine / anchovy / sst

species12=data.frame(1:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, 80, by = 10)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = 3, lib_column = "sp1", 
                        target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples)
                   #can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = 3, lib_column = "sp2", target_column = "sp1", 
                        lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1
Pval_1xmap2 #0.1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2
Pval_2xmap1 #0

### Let's try to see those
sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

plot(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")
legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

########################################################################################################################
########## Sugihara two species deterministic competition model -- no interactions
########################################################################################################################

tmax=300
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=0.1
Y[1]=0.2
for (t in 1:(tmax-1))
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.0*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.0*X[t])
}
X
Y
x=log(X)
y=log(Y)
# centering
x=x-mean(x)
y=y-mean(y)

z=cbind(x,y)
vec_col=c("black","red")
#pdf(file="CompetitionDeterministicNL2Species.pdf",width=12,height=6)
par(cex=1.5)
minz=min(z[151:201,])
maxz=max(z[151:201,])
for (i in 1:2){
  if (i==1){plot(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz),xlab="Time",ylab="ln(Density)")}
  lines(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz))
}
#dev.off()

species12=data.frame(1:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, 80, by = 10)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = 3, lib_column = "sp1", 
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = 3, lib_column = "sp2", target_column = "sp1", 
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1
Pval_1xmap2 #0.49

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2
Pval_2xmap1 #0.47

### Let's try and see those
sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

plot(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")
legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

# How do I know there's no increase? Granted, the numbers are low... Should I set a rho_min?

### Let's try with a different cross-validation setup, 200 for the model and 100 for validation

species12=data.frame(1:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(201, 300, by = 1)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, lib = c(1,200), pred=c(201,300), E = 3, lib_column = "sp1", 
                    target_column = "sp2", lib_sizes = libsizes,num_samples = numsamples)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, lib = c(1,200), pred=c(201,300), E = 3, lib_column = "sp2", target_column = "sp1", 
                    lib_sizes = libsizes, num_samples = numsamples) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1
Pval_1xmap2 #0.41

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2
Pval_2xmap1 #0.46

### Let's try and see those
sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

plot(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")
legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)

# How do I know there's no increase? 

### --- I guess a rho_min around 0.1 is needed to disquality situations in which there's no causality whatsoever -- 

########################################################################################################################
####################### Repeating many times with many initial conditions - yet to do ####################################
########################################################################################################################

ncond<-500 #how many initial conditions do we consider? 

#Initializing vectors 
Pval_12_inter=Pval_21_inter=Pval_12_noInter=Pval_21_noInter=rep(NA,ncond)
RhoLMax_12_inter=RhoLMax_21_inter=RhoLMax_12_noInter=RhoLMax_21_noInter=rep(NA,ncond)
index_1cause2_inter=index_2cause1_inter=index_1cause2_noInter=index_2cause1_noInter=rep(NA,ncond)
lag_order_inter=lag_order_noInter=rep(NA,ncond)

#F: NB I need to replace the lag order by the embedding dimension, so I need to really introduce a code for the simplex

########################################################################################################################
########## Sugihara two species deterministic competition model -- interactions
########################################################################################################################


for (kcond in 1:ncond){
kcond
  ### Data simulation
tmax=800
tmin=501
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=runif(1,0.1,0.7)
Y[1]=runif(1,0.1,0.7)
for (t in 1:tmax)
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}
X
Y
x=log(X[tmin:800])
y=log(Y[tmin:800])
# centering
x=x-mean(x)
y=y-mean(y)
z=cbind(x,y)

 ### determination of the lag order
 # Here we just assume this one to be 3
 lag_order_inter[kcond] = 3
 
 ### CCM Analysis 

species12=data.frame(tmin:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, 80, by = 10)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter[kcond] , lib_column = "sp1", 
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter[kcond] , lib_column = "sp2", target_column = "sp1", 
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1
Pval_1xmap2 #0.93

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2
Pval_2xmap1 #1

### Let's try to see those
sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

plot(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")
legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)


Pval_12_inter[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter[kcond]<0.1)&(RhoLMax_12_inter[kcond]>0.1))
{index_1cause2_inter[kcond]=1} else {index_1cause2_inter[kcond]=0}

if ((Pval_21_inter[kcond]<0.1)&(RhoLMax_21_inter[kcond]>0.1))
{index_2cause1_inter[kcond]=1} else {index_2cause1_inter[kcond]=0}

}


########################################################################################################################
########## Sugihara two species deterministic competition model -- no interactions
########################################################################################################################


for (kcond in 1:ncond){
  kcond
  tmax=800
  X=Y=rep(1,tmax)## Problem if I start with 1
  X[1]=runif(1,0.1,0.7)
  Y[1]=runif(1,0.1,0.7)
  for (t in 1:tmax)
  {
    X[t+1]=X[t]*(3.8-3.8*X[t]-0*Y[t])
    Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0*X[t])
  }
  X
  Y
  x=log(X[501:800])
  y=log(Y[501:800])
  # centering
  x=x-mean(x)
  y=y-mean(y)
  z=cbind(x,y)
  
  ### determination of the lag order
  # Here we just assume this one to be 3
  lag_order_noInter[kcond] = 3
  
  ### CCM Analysis 
  
  species12=data.frame(tmin:tmax,z)
  names(species12)=c("time","sp1","sp2")
  libsizes = seq(10, 80, by = 10)
  lm=length(libsizes)
  numsamples = 100
  sp1_xmap_sp2 <- ccm(species12, E = lag_order_noInter[kcond] , lib_column = "sp1", 
                      target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples)
  #can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
  sp2_xmap_sp1 <- ccm(species12, E = lag_order_noInter[kcond] , lib_column = "sp2", target_column = "sp1", 
                      lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?
  
  ### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
  rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
  rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
  # Fraction of samples for which rho(L_max)<rho(L_min)
  Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1
  Pval_1xmap2 #0.93
  
  rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
  rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]
  
  Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2
  Pval_2xmap1 #1
  
  ### Let's try to see those
  sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
  sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)
  
  plot(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), type = "l", col = "red", 
       xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
  lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")
  legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
  
  
  Pval_12_noInter[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
  Pval_21_noInter[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2
  
  RhoLMax_12_noInter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
  RhoLMax_21_noInter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2
  
  if ((Pval_12_noInter[kcond]<0.1)&(RhoLMax_12_noInter[kcond]>0.1))
  {index_1cause2_noInter[kcond]=1} else {index_1cause2_noInter[kcond]=0}
  
  if ((Pval_21_noInter[kcond]<0.1)&(RhoLMax_21_noInter[kcond]>0.1))
  {index_2cause1_noInter[kcond]=1} else {index_2cause1_noInter[kcond]=0}
  
}

#Fraction of 1->2 when present
sum(index_1cause2_inter)/ncond
#Fraction of 2->1 when present
sum(index_2cause1_inter)/ncond

#Fraction of 1->2 when absent
sum(index_1cause2_noInter)/ncond
#Fraction of 2->1 when absent
sum(index_2cause1_noInter)/ncond

#### What would have happened if we only used P-values?

#Fraction of 1->2 when present
sum(Pval_12_inter<0.1)/ncond
#Fraction of 2->1 when present
sum(Pval_21_inter<0.1)/ncond

#Fraction of 1->2 when absent
sum(Pval_12_noInter<0.1)/ncond
#Fraction of 2->1 when absent
sum(Pval_21_noInter<0.1)/ncond
Pval_21_noInter

# Almost the same thing

DataCompet_sugiharaDeterModel_CCM = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)


########################################################################################################################
########## Stochastic competition model -- interactions
########################################################################################################################

### -- to fill --- ###

########################################################################################################################
########## Stochastic competition model -- no interactions
########################################################################################################################

### -- to fill --- ###