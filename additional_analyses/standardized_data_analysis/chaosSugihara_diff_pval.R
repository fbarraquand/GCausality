########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################
### FB 14/08/2018  - We now include a stochastic model in this code
### CP 10/04/2019 - We now test GC and CCM at the same time, and save evg, as well as plot results
### CP 16/04/19 - We now keep the GC coefficients and compute CCM-pvalues on surrogates
### CP 28/05/19 - Checking standardization does not change the results of centered-only
rm(list=ls())
graphics.off()

library("vars")
library("rEDM")
########################################################################################################################
####################### Repeating many times with many initial conditions / stochastic sequences    ####################
########################################################################################################################
set.seed(42)
#Initializing vectors 
ncond<-500
Pval_12_inter_GC=Pval_21_inter_GC=Pval_12_noInter_GC=Pval_21_noInter_GC=rep(NA,ncond)
Pval_12_inter_CCM=Pval_21_inter_CCM=Pval_12_noInter_CCM=Pval_21_noInter_CCM=rep(NA,ncond)
Pval_12_inter_CCM_surr=Pval_21_inter_CCM_surr=Pval_12_noInter_CCM_surr=Pval_21_noInter_CCM_surr=rep(NA,ncond)
Pval_12_inter_CCM_surr_ebi=Pval_21_inter_CCM_surr_ebi=Pval_12_noInter_CCM_surr_ebi=Pval_21_noInter_CCM_surr_ebi=rep(NA,ncond)
Pval_12_inter_CCM_surr_twin=Pval_21_inter_CCM_surr_twin=Pval_12_noInter_CCM_surr_twin=Pval_21_noInter_CCM_surr_twin=rep(NA,ncond)

log_12_inter=log_21_inter=log_12_noInter=log_21_noInter=rep(NA,ncond)
effect_12_inter=effect_21_inter=effect_12_noInter=effect_21_noInter=rep(NA,ncond)
RhoLMax_12_inter_v1=RhoLMax_21_inter_v1=RhoLMax_12_noInter_v1=RhoLMax_21_noInter_v1=rep(NA,ncond) #v1 is from the bootstrap
RhoLMax_12_inter_v2=RhoLMax_21_inter_v2=RhoLMax_12_noInter_v2=RhoLMax_21_noInter_v2=rep(NA,ncond) #v2 is from the total time series and only it

index_1cause2_inter_GC=index_2cause1_inter_GC=index_1cause2_noInter_GC=index_2cause1_noInter_GC=rep(NA,ncond)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=index_1cause2_noInter_CCM=index_2cause1_noInter_CCM=rep(NA,ncond)

lag_order_inter_GC=lag_order_noInter_GC=rep(NA,ncond)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=lag_order_noInter_CCM_predictx=lag_order_noInter_CCM_predicty=rep(NA,ncond)

########################################################################################################################
########## Chaotic two species competition model -- interactions present
########################################################################################################################
#if(1==0){
print("Inter")
for (kcond in 1:ncond){
print(kcond)
tmax=800
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=runif(1,0.1,0.7)
Y[1]=runif(1,0.1,0.7)

for (t in 1:tmax)
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}

x=log(X[501:800])
y=log(Y[501:800])
# centering
#x=x-mean(x)
#y=y-mean(y)
x=c(scale(x))
y=c(scale(y))
z=cbind(x,y)

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[kcond] <- varcompet$p

#Let's compute log ratio
ar_x=ar(x,order=varcompet$p,AIC=F,method="ols")
ar_y=ar(y,order=varcompet$p,AIC=F,method="ols")

log_12_inter[kcond]=log(sum((ar_y$resid)^2,na.rm=T)/sum((varcompet$varresult$y$residuals)^2,na.rm=T))
log_21_inter[kcond]=log(sum((ar_x$resid)^2,na.rm=T)/sum((varcompet$varresult$x$residuals)^2,na.rm=T))

n1=names(varcompet$varresult$x$coefficients)
effect_21_inter[kcond]=sum(abs(varcompet$varresult$x$coefficients[grep("y",n1)]))
n2=names(varcompet$varresult$y$coefficients)
effect_12_inter[kcond]=sum(abs(varcompet$varresult$y$coefficients[grep("x",n2)]))


gxy = grangertest(x,y,order = lag_order_inter_GC[kcond]) #x causes y 
gyx = grangertest(y,x,order = lag_order_inter_GC[kcond]) #y causes x

Pval_12_inter_GC[kcond]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[kcond]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[kcond]<0.1)
{index_1cause2_inter_GC[kcond]=1} else {index_1cause2_inter_GC[kcond]=0}

if (Pval_21_inter_GC[kcond]<0.1)
{index_2cause1_inter_GC[kcond]=1} else {index_2cause1_inter_GC[kcond]=0}


###Let's CCM
#Chose E
simplex_output_predictx = simplex(x,E=1:10)
 lag_order_inter_CCM_predictx[kcond] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_inter_CCM_predicty[kcond] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
libsizes = c(5,8,seq(10, nrow(species12)-11, by = 10))
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs=T,num_samples = numsamples,replace=T)  #We want the real bootstrap from CB
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs =T, num_samples = numsamples,replace=T) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)


### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==max(sp1_xmap_sp2$lib_size)]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==max(sp2_xmap_sp1$lib_size)]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter_v1[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==max(sp2_xmap_sp1_means$lib_size)] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter_v1[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==max(sp1_xmap_sp2_means$lib_size)] # 2 causes 1 if 1 xmap 2

####Doing a test here, I really want RhoLMax with all points
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), num_samples = 1,replace=F)  #We want the real bootstrap from CB
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = max(sp2_xmap_sp1$lib_size), num_samples = 1,replace=F)

RhoLMax_12_inter_v2[kcond]=sp2_xmap_sp1$rho # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter_v2[kcond]=sp1_xmap_sp2$rho # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond]<0.1)&(RhoLMax_12_inter_v1[kcond]>0.1))
{index_1cause2_inter_CCM[kcond]=1} else {index_1cause2_inter_CCM[kcond]=0}

if ((Pval_21_inter_CCM[kcond]<0.1)&(RhoLMax_21_inter_v1[kcond]>0.1))
{index_2cause1_inter_CCM[kcond]=1} else {index_2cause1_inter_CCM[kcond]=0}

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp2"]=sample(species12[,"sp2"])
        sp1_xmap_sp2_random <- ccm(species_random, E = lag_order_inter_CCM_predictx[kcond], lib_column = "sp1",target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp1_xmap_sp2_random$rho
}
  Pval_21_inter_CCM_surr[kcond] = sum(rho_dist>RhoLMax_21_inter_v2[kcond])/numsamples

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_inter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_12_inter_CCM_surr[kcond] = sum(rho_dist>RhoLMax_12_inter_v2[kcond])/numsamples


####Surrogates
surr_twin_s1 <- make_surrogate_data(species12$sp1, method = "twin", num_surr = numsamples,phase_lock=F)
surr_twin_s2 <- make_surrogate_data(species12$sp2, method = "twin", num_surr = numsamples,phase_lock=F)
surr_ebi_s1 <- make_surrogate_data(species12$sp1, method = "ebisuzaki", num_surr = numsamples)
surr_ebi_s2 <- make_surrogate_data(species12$sp2, method = "ebisuzaki", num_surr = numsamples)


rho_twin<- list(species1 =  matrix(NA,nrow=numsamples,ncol=1), species2 =  matrix(NA,nrow=numsamples,ncol=1))
rho_ebi<- list(species1 =  matrix(NA,nrow=numsamples,ncol=1), species2 =  matrix(NA,nrow=numsamples,ncol=1))

for (i in 1:numsamples) {
  rho_twin$species1[i,1] <- ccm(cbind(species12$sp2, surr_twin_s1[,i]), E =  lag_order_inter_CCM_predicty[kcond], lib_column = 1, target_column = 2, lib_sizes =max(sp2_xmap_sp1$lib_size),num_samples=1,replace=F)$rho #species 1 the cause so it is the surrogated TS and we use the embedding for species2
  rho_twin$species2[i,1] <- ccm(cbind(surr_twin_s2[,i], species12$sp1), E =  lag_order_inter_CCM_predictx[kcond], lib_column = 2, target_column = 1, lib_sizes = max(sp1_xmap_sp2$lib_size),num_samples=1,replace=F)$rho #species2 is te cause so it's the target and surrogated TS, and we use the embedding for species 2
  rho_ebi$species1[i,1] <- ccm(cbind(species12$sp2, surr_ebi_s1[,i]), E =  lag_order_inter_CCM_predicty[kcond], lib_column = 1, target_column = 2, lib_sizes =max(sp2_xmap_sp1$lib_size),num_samples=1,replace=F)$rho #species 1 the cause so it is the surrogated TS and we use the embedding for species2
  rho_ebi$species2[i,1] <- ccm(cbind(surr_ebi_s2[,i], species12$sp1), E =  lag_order_inter_CCM_predictx[kcond], lib_column = 2, target_column = 1, lib_sizes = max(sp1_xmap_sp2$lib_size),num_samples=1,replace=F)$rho #species2 is te cause so it's the target and surrogated TS, and we use the embedding for species 2
}

Pval_21_inter_CCM_surr_twin[kcond]=sum(RhoLMax_21_inter_v2[kcond]<rho_twin$species2) /numsamples
Pval_12_inter_CCM_surr_twin[kcond]=sum(RhoLMax_12_inter_v2[kcond]<rho_twin$species1) /numsamples
Pval_21_inter_CCM_surr_ebi[kcond]=sum(RhoLMax_21_inter_v2[kcond]<rho_ebi$species2) /numsamples
Pval_12_inter_CCM_surr_ebi[kcond]=sum(RhoLMax_12_inter_v2[kcond]<rho_ebi$species1) /numsamples


}
#}
########################################################################################################################
########## Chaotic two species competition model -- no interactions
########################################################################################################################
print("no inter")
for (kcond in 1:ncond){
  print(kcond)
  tmax=800
  X=Y=rep(1,tmax)## Problem if I start with 1
  X[1]=runif(1,0.1,0.7)
  Y[1]=runif(1,0.1,0.7)
  for (t in 1:tmax)
  {
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.0*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.0*X[t])
  }
  x=log(X[501:800])
  y=log(Y[501:800])
  # centering
#  x=x-mean(x)
#  y=y-mean(y)
  x=c(scale(x))
  y=c(scale(y))
  z=cbind(x,y)
  
  varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
  lag_order_noInter_GC[kcond] <- varcompet$p

#Let's compute log ratio
ar_x=ar(x,order=varcompet$p,AIC=F,method="ols")
ar_y=ar(y,order=varcompet$p,AIC=F,method="ols")
log_12_noInter[kcond]=log(sum((ar_y$resid)^2,na.rm=T)/sum((varcompet$varresult$y$residuals)^2,na.rm=T))
log_21_noInter[kcond]=log(sum((ar_x$resid)^2,na.rm=T)/sum((varcompet$varresult$x$residuals)^2,na.rm=T))


n1=names(varcompet$varresult$x$coefficients)
effect_21_noInter[kcond]=sum(abs(varcompet$varresult$x$coefficients[grep("y",n1)]))
n2=names(varcompet$varresult$y$coefficients)
effect_12_noInter[kcond]=sum(abs(varcompet$varresult$y$coefficients[grep("x",n2)]))

tryCatch({  
  gxy = grangertest(x,y,order = lag_order_noInter_GC[kcond]) #x causes y 
  gyx = grangertest(y,x,order = lag_order_noInter_GC[kcond]) #y causes x
  
  Pval_12_noInter_GC[kcond]=gxy$`Pr(>F)`[2]
  Pval_21_noInter_GC[kcond]=gyx$`Pr(>F)`[2]

if (Pval_12_noInter_GC[kcond]<0.1)
{index_1cause2_noInter_GC[kcond]=1} else {index_1cause2_noInter_GC[kcond]=0}

if (Pval_21_noInter_GC[kcond]<0.1)
{index_2cause1_noInter_GC[kcond]=1} else {index_2cause1_noInter_GC[kcond]=0}
},error=function(e){
	return(0)
})

#Let's CCM
simplex_output_predictx = simplex(x,E=1:10)
 lag_order_noInter_CCM_predictx[kcond] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_noInter_CCM_predicty[kcond] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
libsizes = c(5,8,seq(10, nrow(species12)-11, by = 10))
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_noInter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_size = libsizes, num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_noInter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_size = libsizes, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==max(sp1_xmap_sp2$lib_size)]
# Fraction of samples for which rho(L_max)<rho(L_min)
tmp= na.omit(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)
Pval_1xmap2 = sum(tmp)/length(tmp) #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==max(sp2_xmap_sp1$lib_size)]
tmp=na.omit(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)
Pval_2xmap1 = sum(tmp)/length(tmp) #1 towards 2

if(Pval_2xmap1>1|Pval_1xmap2>1){stop("Euh..")}

Pval_12_noInter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_noInter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_noInter_v1[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==max(sp2_xmap_sp1_means$lib_size)] # 1 causes 2 if 2 xmap 1
RhoLMax_21_noInter_v1[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==max(sp1_xmap_sp2_means$lib_size)] # 2 causes 1 if 1 xmap 2

####Doing a test here, I really want RhoLMax with all points
sp1_xmap_sp2 <- ccm(species12, E = lag_order_noInter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), num_samples = 1,replace=F)  #We want the real bootstrap from CB
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_noInter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = max(sp2_xmap_sp1$lib_size), num_samples = 1,replace=F)

RhoLMax_12_noInter_v2[kcond]=sp2_xmap_sp1$rho # 1 causes 2 if 2 xmap 1
RhoLMax_21_noInter_v2[kcond]=sp1_xmap_sp2$rho # 2 causes 1 if 1 xmap 2


if ((Pval_12_noInter_CCM[kcond]<0.1)&(RhoLMax_12_noInter_v1[kcond]>0.1))
{index_1cause2_noInter_CCM[kcond]=1} else {index_1cause2_noInter_CCM[kcond]=0}

if ((Pval_21_noInter_CCM[kcond]<0.1)&(RhoLMax_21_noInter_v1[kcond]>0.1))
{index_2cause1_noInter_CCM[kcond]=1} else {index_2cause1_noInter_CCM[kcond]=0}

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp2"]=sample(species12[,"sp2"])
        sp1_xmap_sp2_random <- ccm(species_random, E = lag_order_noInter_CCM_predictx[kcond], lib_column = "sp1",target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp1_xmap_sp2_random$rho
}
  Pval_21_noInter_CCM_surr[kcond] = sum(rho_dist>RhoLMax_21_noInter_v2[kcond],na.rm=T)/numsamples

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_noInter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_12_noInter_CCM_surr[kcond] = sum(rho_dist>RhoLMax_12_noInter_v2[kcond])/numsamples

####Surrogates
surr_twin_s1 <- make_surrogate_data(species12$sp1, method = "twin", num_surr = numsamples,phase_lock=F)
surr_twin_s2 <- make_surrogate_data(species12$sp2, method = "twin", num_surr = numsamples,phase_lock=F)
surr_ebi_s1 <- make_surrogate_data(species12$sp1, method = "ebisuzaki", num_surr = numsamples)
surr_ebi_s2 <- make_surrogate_data(species12$sp2, method = "ebisuzaki", num_surr = numsamples)


rho_twin<- list(species1 =  matrix(NA,nrow=numsamples,ncol=1), species2 =  matrix(NA,nrow=numsamples,ncol=1))
rho_ebi<- list(species1 =  matrix(NA,nrow=numsamples,ncol=1), species2 =  matrix(NA,nrow=numsamples,ncol=1))

for (i in 1:numsamples) {
  rho_twin$species1[i,1] <- ccm(cbind(species12$sp2, surr_twin_s1[,i]), E =  lag_order_noInter_CCM_predicty[kcond], lib_column = 1, target_column = 2, lib_sizes =max(sp2_xmap_sp1$lib_size),num_samples=1,replace=F)$rho #species 1 the cause so it is the surrogated TS and we use the embedding for species2
  rho_twin$species2[i,1] <- ccm(cbind(surr_twin_s2[,i], species12$sp1), E =  lag_order_noInter_CCM_predictx[kcond], lib_column = 2, target_column = 1, lib_sizes = max(sp1_xmap_sp2$lib_size),num_samples=1,replace=F)$rho #species2 is te cause so it's the target and surrogated TS, and we use the embedding for species 2
  rho_ebi$species1[i,1] <- ccm(cbind(species12$sp2, surr_ebi_s1[,i]), E =  lag_order_noInter_CCM_predicty[kcond], lib_column = 1, target_column = 2, lib_sizes =max(sp2_xmap_sp1$lib_size),num_samples=1,replace=F)$rho #species 1 the cause so it is the surrogated TS and we use the embedding for species2
  rho_ebi$species2[i,1] <- ccm(cbind(surr_ebi_s2[,i], species12$sp1), E =  lag_order_noInter_CCM_predictx[kcond], lib_column = 2, target_column = 1, lib_sizes = max(sp1_xmap_sp2$lib_size),num_samples=1,replace=F)$rho #species2 is te cause so it's the target and surrogated TS, and we use the embedding for species 2
}

Pval_21_noInter_CCM_surr_twin[kcond]=sum(RhoLMax_21_noInter_v2[kcond]<rho_twin$species2) /numsamples
Pval_12_noInter_CCM_surr_twin[kcond]=sum(RhoLMax_12_noInter_v2[kcond]<rho_twin$species1) /numsamples
Pval_21_noInter_CCM_surr_ebi[kcond]=sum(RhoLMax_21_noInter_v2[kcond]<rho_ebi$species2) /numsamples
Pval_12_noInter_CCM_surr_ebi[kcond]=sum(RhoLMax_12_noInter_v2[kcond]<rho_ebi$species1) /numsamples



}
DataCompet_stochModel_inter = data.frame(1:ncond,lag_order_inter_GC,Pval_12_inter_GC,Pval_21_inter_GC,index_1cause2_inter_GC,index_2cause1_inter_GC,effect_12_inter,effect_21_inter,log_12_inter,log_21_inter,lag_order_inter_CCM_predictx,Pval_12_inter_CCM,lag_order_inter_CCM_predicty,Pval_21_inter_CCM,index_1cause2_inter_CCM,index_2cause1_inter_CCM,Pval_12_inter_CCM_surr,Pval_21_inter_CCM_surr,Pval_12_inter_CCM_surr_twin,Pval_21_inter_CCM_surr_twin,Pval_12_inter_CCM_surr_ebi,Pval_21_inter_CCM_surr_ebi,RhoLMax_12_inter_v1,RhoLMax_21_inter_v1,RhoLMax_12_inter_v2,RhoLMax_21_inter_v2)
DataCompet_stochModel_noInter = data.frame(1:ncond,lag_order_noInter_GC,Pval_12_noInter_GC,Pval_21_noInter_GC,index_1cause2_noInter_GC,index_2cause1_noInter_GC,effect_12_noInter,effect_21_noInter,log_12_noInter,log_21_noInter,lag_order_noInter_CCM_predictx,Pval_12_noInter_CCM,lag_order_noInter_CCM_predicty,Pval_21_noInter_CCM,index_1cause2_noInter_CCM,index_2cause1_noInter_CCM,Pval_12_noInter_CCM_surr,Pval_21_noInter_CCM_surr,Pval_12_noInter_CCM_surr_twin,Pval_21_noInter_CCM_surr_twin,Pval_12_noInter_CCM_surr_ebi,Pval_21_noInter_CCM_surr_ebi,RhoLMax_12_noInter_v1,RhoLMax_21_noInter_v1,RhoLMax_12_noInter_v2,RhoLMax_21_noInter_v2)
#Write down results
write.csv(DataCompet_stochModel_inter,file="results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
write.csv(DataCompet_stochModel_noInter,file="results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")

