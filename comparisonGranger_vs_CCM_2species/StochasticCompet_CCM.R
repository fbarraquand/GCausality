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

### FB 14/08/2018  - We now include a stochastic model in this code

##############################################################################################################
####################### Repeating many times with many initial conditions ####################################
##############################################################################################################

ncond<-500 #how many initial conditions we consider

#Initializing vectors 
Pval_12_inter=Pval_21_inter=Pval_12_noInter=Pval_21_noInter=rep(NA,ncond)
RhoLMax_12_inter=RhoLMax_21_inter=RhoLMax_12_noInter=RhoLMax_21_noInter=rep(NA,ncond)
index_1cause2_inter=index_2cause1_inter=index_1cause2_noInter=index_2cause1_noInter=rep(NA,ncond)
lag_order_inter=lag_order_noInter=rep(NA,ncond)

#F: NB replace the lag order by the embedding dimension, need to introduce a code for the simplex? 

######################################################################################################
########## Stochastic two species competition model -- interactions
#######################################################################################################


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
  X[t+1]=X[t]*exp(3-4*X[t]-2*Y[t]+ rnorm(1,0,0.1))
  Y[t+1]=Y[t]*exp(2.1-3.1*Y[t]-0.31*X[t]+ rnorm(1,0,0.1))
}
X
Y
x=X[tmin:800]#log(X[tmin:800]) #no log-transfo for CCM to avoid weird stuff
y=Y[tmin:800]#log(Y[tmin:800])
# centering
x=(x-mean(x))/sd(x)
y=(y-mean(y))/sd(y)
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


#########################################################################################################
########## Stochastic two species competition model -- no interactions
#########################################################################################################


for (kcond in 1:ncond){
  kcond
  tmax=800
  X=Y=rep(1,tmax)## Problem if I start with 1
  X[1]=runif(1,0.1,0.7)
  Y[1]=runif(1,0.1,0.7)
  for (t in 1:tmax)
  {
    X[t+1]=X[t]*exp(3-4*X[t]-0*Y[t]+ rnorm(1,0,0.1))
    Y[t+1]=Y[t]*exp(2.1-3.1*Y[t]-0*X[t]+ rnorm(1,0,0.1))
  }
  X
  Y
  x=X[tmin:800]#log(X[tmin:800]) #no log-transfo for CCM to avoid weird stuff
  y=Y[tmin:800]#log(Y[tmin:800])
  # centering
  x=(x-mean(x))/sd(x)
  y=(y-mean(y))/sd(y)
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

# Remark -- although it does not work at all with P-values, the coefficient rho is high whenever 
# there is an interaction and low whenever there is none. 

DataCompet_stochModel_CCM = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter,RhoLMax_12_inter,RhoLMax_21_inter,RhoLMax_12_noInter,RhoLMax_21_noInter)
#Write down results
write.csv(DataCompet_stochModel_CCM,file="results/DataCompet_stochModel_CCM.csv")
