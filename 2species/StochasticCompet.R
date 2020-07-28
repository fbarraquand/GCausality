########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################
### FB 14/08/2018  - We now include a stochastic model in this code
### CP 10/04/2019 - We now test GC and CCM at the same time, and save evg, as well as plot results
### CP 16/04/2019 - We now keep the GC coefficients and compute CCM-pvalues on surrogates
### CP 08/07/2019 - We now use pvalue as (r+1)/(n+1)

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

tab_simu=array(NA,dim=c(ncond,30,2,2))

########################################################################################################################
########## Stochastic two species competition model -- interactions present
########################################################################################################################
print("Inter")
for (kcond in 1:ncond){
print(kcond)
tmax=800
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=runif(1,0.1,0.7)
Y[1]=runif(1,0.1,0.7)
for (t in 1:tmax)
{
  X[t+1]=X[t]*exp(3-4*X[t]-2*Y[t]+ rnorm(1,0,0.1))
  Y[t+1]=Y[t]*exp(2.1-3.1*Y[t]-0.31*X[t]+ rnorm(1,0,0.1))
}
x=log(X[501:800])
y=log(Y[501:800])
# centering
x=x-mean(x)
y=y-mean(y)
z=cbind(x,y)

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[kcond] <- varcompet$p

#Let's compute log ratio
ar_x=ar(x,order=varcompet$p,aic=T,method="ols")
ar_y=ar(y,order=varcompet$p,aic=T,method="ols")

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


if(kcond==1){
pdf("fig/stochastic_competition_model.pdf",width=10)
par(cex=1.25,lwd=2)
plot(551:600,x[51:100],col="blue",pch=16,xlab="Time",ylab="ln(abundance)",t="o",lty=1,ylim=c(-2.75,3.0))
lines(551:600,y[51:100],col="red",pch=16,t="o")
legend("topleft",c("x inter","y inter","x no inter","y no inter"),col=c("blue","red",rgb(255/256,165/256,0,1),rgb(155/256,79/256,150/256,1)),pch=16,lty=1,bty="n")
}
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
#libsizes = c(seq(3, 90, by = 3)) #only for precise
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs=T,num_samples = numsamples,replace=T)  #We want the real bootstrap from CB
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs =T, num_samples = numsamples,replace=T) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

tab_simu[kcond,,1,1]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2,1]=sp2_xmap_sp1_means$rho

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
  Pval_21_inter_CCM_surr[kcond] = (sum(rho_dist>=RhoLMax_21_inter_v2[kcond])+1)/(numsamples+1)

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_inter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_12_inter_CCM_surr[kcond] = (sum(rho_dist>=RhoLMax_12_inter_v2[kcond])+1)/(numsamples+1)


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

Pval_21_inter_CCM_surr_twin[kcond]=(sum(RhoLMax_21_inter_v2[kcond]<=rho_twin$species2)+1) /(numsamples+1)
Pval_12_inter_CCM_surr_twin[kcond]=(sum(RhoLMax_12_inter_v2[kcond]<=rho_twin$species1)+1) /(numsamples+1)
Pval_21_inter_CCM_surr_ebi[kcond]=(sum(RhoLMax_21_inter_v2[kcond]<=rho_ebi$species2)+1) /(numsamples+1)
Pval_12_inter_CCM_surr_ebi[kcond]=(sum(RhoLMax_12_inter_v2[kcond]<=rho_ebi$species1)+1) /(numsamples+1)


}
########################################################################################################################
########## Stochastic two species competition model -- no interactions
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
    X[t+1]=X[t]*exp(3-4*X[t]-0*Y[t]+ rnorm(1,0,0.1))
    Y[t+1]=Y[t]*exp(2.1-3.1*Y[t]-0*X[t]+ rnorm(1,0,0.1))
  }
  x=log(X[501:800])
  y=log(Y[501:800])
  # centering
  x=x-mean(x)
  y=y-mean(y)
  z=cbind(x,y)
  
  varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
  lag_order_noInter_GC[kcond] <- varcompet$p


#Let's compute log ratio
ar_x=ar(x,order=varcompet$p,aic=F,method="ols")
ar_y=ar(y,order=varcompet$p,aic=F,method="ols")

log_12_noInter[kcond]=log(sum((ar_y$resid)^2,na.rm=T)/sum((varcompet$varresult$y$residuals)^2,na.rm=T))
log_21_noInter[kcond]=log(sum((ar_x$resid)^2,na.rm=T)/sum((varcompet$varresult$x$residuals)^2,na.rm=T))


n1=names(varcompet$varresult$x$coefficients)
effect_21_noInter[kcond]=sum(abs(varcompet$varresult$x$coefficients[grep("y",n1)]))
n2=names(varcompet$varresult$y$coefficients)
effect_12_noInter[kcond]=sum(abs(varcompet$varresult$y$coefficients[grep("x",n2)]))

  
  gxy = grangertest(x,y,order = lag_order_noInter_GC[kcond]) #x causes y 
  gyx = grangertest(y,x,order = lag_order_noInter_GC[kcond]) #y causes x
  
  Pval_12_noInter_GC[kcond]=gxy$`Pr(>F)`[2]
  Pval_21_noInter_GC[kcond]=gyx$`Pr(>F)`[2]

if (Pval_12_noInter_GC[kcond]<0.1)
{index_1cause2_noInter_GC[kcond]=1} else {index_1cause2_noInter_GC[kcond]=0}

if (Pval_21_noInter_GC[kcond]<0.1)
{index_2cause1_noInter_GC[kcond]=1} else {index_2cause1_noInter_GC[kcond]=0}


if(kcond==1){
lines(551:600,y[51:100],col=rgb(255/256,165/256,0,1),pch=16,lty=1,t="o")
lines(551:600,x[51:100],col= rgb(155/256,79/256,150/256,1),pch=16,t="o")
dev.off()

}

#Let's CCM
simplex_output_predictx = simplex(x,E=1:10)
 lag_order_noInter_CCM_predictx[kcond] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_noInter_CCM_predicty[kcond] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
#libsizes = c(seq(3, 90, by = 3))
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_noInter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_size = libsizes, num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_noInter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_size = libsizes, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

tab_simu[kcond,,1,2]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2,2]=sp2_xmap_sp1_means$rho

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==max(sp1_xmap_sp2$lib_size)]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==max(sp2_xmap_sp1$lib_size)]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

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
  Pval_21_noInter_CCM_surr[kcond] = (sum(rho_dist>=RhoLMax_21_noInter_v2[kcond])+1)/(numsamples+1)

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_noInter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_12_noInter_CCM_surr[kcond] = (sum(rho_dist>=RhoLMax_12_noInter_v2[kcond])+1)/(numsamples+1)

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

Pval_21_noInter_CCM_surr_twin[kcond]=(sum(RhoLMax_21_noInter_v2[kcond]<rho_twin$species2)+1) /(numsamples+1)
Pval_12_noInter_CCM_surr_twin[kcond]=(sum(RhoLMax_12_noInter_v2[kcond]<rho_twin$species1)+1)/(numsamples+1)
Pval_21_noInter_CCM_surr_ebi[kcond]=(sum(RhoLMax_21_noInter_v2[kcond]<rho_ebi$species2)+1) /(numsamples+1)
Pval_12_noInter_CCM_surr_ebi[kcond]=(sum(RhoLMax_12_noInter_v2[kcond]<rho_ebi$species1)+1) /(numsamples+1)


}
DataCompet_stochModel_inter = data.frame(1:ncond,lag_order_inter_GC,Pval_12_inter_GC,Pval_21_inter_GC,index_1cause2_inter_GC,index_2cause1_inter_GC,effect_12_inter,effect_21_inter,log_12_inter,log_21_inter,lag_order_inter_CCM_predictx,Pval_12_inter_CCM,lag_order_inter_CCM_predicty,Pval_21_inter_CCM,index_1cause2_inter_CCM,index_2cause1_inter_CCM,Pval_12_inter_CCM_surr,Pval_21_inter_CCM_surr,Pval_12_inter_CCM_surr_twin,Pval_21_inter_CCM_surr_twin,Pval_12_inter_CCM_surr_ebi,Pval_21_inter_CCM_surr_ebi,RhoLMax_12_inter_v1,RhoLMax_21_inter_v1,RhoLMax_12_inter_v2,RhoLMax_21_inter_v2)

DataCompet_stochModel_noInter = data.frame(1:ncond,lag_order_noInter_GC,Pval_12_noInter_GC,Pval_21_noInter_GC,index_1cause2_noInter_GC,index_2cause1_noInter_GC,effect_12_noInter,effect_21_noInter,log_12_noInter,log_21_noInter,lag_order_noInter_CCM_predictx,Pval_12_noInter_CCM,lag_order_noInter_CCM_predicty,Pval_21_noInter_CCM,index_1cause2_noInter_CCM,index_2cause1_noInter_CCM,Pval_12_noInter_CCM_surr,Pval_21_noInter_CCM_surr,Pval_12_noInter_CCM_surr_twin,Pval_21_noInter_CCM_surr_twin,Pval_12_noInter_CCM_surr_ebi,Pval_21_noInter_CCM_surr_ebi,RhoLMax_12_noInter_v1,RhoLMax_21_noInter_v1,RhoLMax_12_noInter_v2,RhoLMax_21_noInter_v2)
#Write down results
write.csv(DataCompet_stochModel_inter,file="results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
write.csv(DataCompet_stochModel_noInter,file="results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")

#pdf("stochasticcompet_ccm_precise.pdf",width=7,height=5) #Change libsizes if you want the precise version
pdf("fig/stochasticcompet_ccm.pdf",width=7,height=5)
tmp=libsizes
plot(tmp,rep(0,length(tmp)),ylim=c(0,1.1),t="n",xlab="Library size",ylab="rho")
for(kcond in 1:ncond){
if(kcond==1){
alpha=1
}else{
alpha=0.1
}
lines(tmp,tab_simu[kcond,,1,1],col=rgb(1,0,0,alpha)) #1xmap2
lines(tmp,tab_simu[kcond,,2,1],col= rgb(0,0,1,alpha)) #2xmap1
}

for(kcond in 1:ncond){
if(kcond==1){
alpha=1
}else{
alpha=0.1
}
lines(tmp,tab_simu[kcond,,1,2],col=rgb(255/256,165/256,0,alpha),lty=1) #1xmap2
lines(tmp,tab_simu[kcond,,2,2],col= rgb(155/256,79/256,150/256,alpha)) #2xmap1
}
alpha=1
legend(x = "right", legend = c("2 causes 1, inter", "1 causes 2, inter","2 causes 1, no inter", "1 causes 2, no inter"), col = c("red","blue",rgb(255/256,165/256,0,alpha),rgb(155/256,79/256,150/256,alpha)), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
dev.off()
