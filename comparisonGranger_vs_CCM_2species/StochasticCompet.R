########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################
### FB 14/08/2018  - We now include a stochastic model in this code

library("vars")
library("rEDM")
########################################################################################################################
####################### Repeating many times with many initial conditions / stochastic sequences    ####################
########################################################################################################################

#Initializing vectors 
ncond<-500
#Initializing vectors 
Pval_12_inter_GC=Pval_21_inter_GC=Pval_12_noInter_GC=Pval_21_noInter_GC=rep(NA,ncond)
Pval_12_inter_CCM=Pval_21_inter_CCM=Pval_12_noInter_CCM=Pval_21_noInter_CCM=rep(NA,ncond)

RhoLMax_12_inter=RhoLMax_21_inter=RhoLMax_12_noInter=RhoLMax_21_noInter=rep(NA,ncond)

index_1cause2_inter_GC=index_2cause1_inter_GC=index_1cause2_noInter_GC=index_2cause1_noInter_GC=rep(NA,ncond)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=index_1cause2_noInter_CCM=index_2cause1_noInter_CCM=rep(NA,ncond)

lag_order_inter_GC=lag_order_noInter_GC=rep(NA,ncond)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=lag_order_noInter_CCM_predictx=lag_order_noInter_CCM_predicty=rep(NA,ncond)

tab_simu=array(NA,dim=c(ncond,15,2,2))

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

gxy = grangertest(x,y,order = lag_order_inter_GC[kcond]) #x causes y 
gyx = grangertest(y,x,order = lag_order_inter_GC[kcond]) #y causes x

Pval_12_inter_GC[kcond]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[kcond]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[kcond]<0.1)
{index_1cause2_inter_GC[kcond]=1} else {index_1cause2_inter_GC[kcond]=0}

if (Pval_21_inter_GC[kcond]<0.1)
{index_2cause1_inter_GC[kcond]=1} else {index_2cause1_inter_GC[kcond]=0}


if(kcond==1){
pdf("stochastic_competition_model.pdf",width=10)
par(cex=1.25,lwd=2)
plot(551:600,x[51:100],col="blue",pch=16,xlab="Time",ylab="ln(abundance)",t="o",lty=1,ylim=c(-2,2))
lines(551:600,y[51:100],col="red",pch=16,t="o")
}
###Let's CCM
#Chose E
smap_output_predictx = simplex(x,E=1:10)
 lag_order_inter_CCM_predictx[kcond] = smap_output_predictx$E[which(smap_output_predictx$rho==max(smap_output_predictx$rho))]

smap_output_predicty = simplex(y,E=1:10)
 lag_order_inter_CCM_predicty[kcond] = smap_output_predicty$E[which(smap_output_predicty$rho==max(smap_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(2, 30, by = 2)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

tab_simu[kcond,,1,1]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2,1]=sp2_xmap_sp1_means$rho

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond]<0.1)&(RhoLMax_12_inter[kcond]>0.1))
{index_1cause2_inter_CCM[kcond]=1} else {index_1cause2_inter_CCM[kcond]=0}

if ((Pval_21_inter_CCM[kcond]<0.1)&(RhoLMax_21_inter[kcond]>0.1))
{index_2cause1_inter_CCM[kcond]=1} else {index_2cause1_inter_CCM[kcond]=0}


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
  
  gxy = grangertest(x,y,order = lag_order_noInter_GC[kcond]) #x causes y 
  gyx = grangertest(y,x,order = lag_order_noInter_GC[kcond]) #y causes x
  
  Pval_12_noInter_GC[kcond]=gxy$`Pr(>F)`[2]
  Pval_21_noInter_GC[kcond]=gyx$`Pr(>F)`[2]

if (Pval_12_noInter_GC[kcond]<0.1)
{index_1cause2_noInter_GC[kcond]=1} else {index_1cause2_noInter_GC[kcond]=0}

if (Pval_21_inter_GC[kcond]<0.1)
{index_2cause1_noInter_GC[kcond]=1} else {index_2cause1_noInter_GC[kcond]=0}


if(kcond==1){
lines(551:600,y[51:100],col=rgb(255/256,165/256,0,1),pch=16,lty=1,t="o")
lines(551:600,x[51:100],col= rgb(155/256,79/256,150/256,1),pch=16,t="o")
dev.off()
}


#Let's CCM
smap_output_predictx = simplex(x,E=1:10)
 lag_order_noInter_CCM_predictx[kcond] = smap_output_predictx$E[which(smap_output_predictx$rho==max(smap_output_predictx$rho))]

smap_output_predicty = simplex(y,E=1:10)
 lag_order_noInter_CCM_predicty[kcond] = smap_output_predicty$E[which(smap_output_predicty$rho==max(smap_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(2, 30, by = 2)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_noInter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_noInter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

tab_simu[kcond,,1,2]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2,2]=sp2_xmap_sp1_means$rho

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_noInter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_noInter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_noInter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_noInter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_noInter_CCM[kcond]<0.1)&(RhoLMax_12_noInter[kcond]>0.1))
{index_1cause2_noInter_CCM[kcond]=1} else {index_1cause2_noInter_CCM[kcond]=0}

if ((Pval_21_noInter_CCM[kcond]<0.1)&(RhoLMax_21_noInter[kcond]>0.1))
{index_2cause1_noInter_CCM[kcond]=1} else {index_2cause1_noInter_CCM[kcond]=0}

}

#par(mfrow=c(1,2))
#sm.density.compare(c(Pval_12_inter,Pval_12_noInter), group = c(1,0), model = "none",lwd=2)
#sm.density.compare(c(Pval_21_inter,Pval_21_noInter), group = c(1,0), model = "none",lwd=2)
### Poor representation

#DataCompet_sugiharaDeterModel = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)


### Let's try the percentage of initial conditions at which the Granger tests works
#sum(Pval_12_inter<0.1)/length(Pval_12_inter) #100% at 0.1 level
#sum(Pval_21_inter<0.1)/length(Pval_21_inter) #52% at 0.1 level

#sum(Pval_12_noInter>0.1)/length(Pval_12_noInter) #90% at 0.1 level
#sum(Pval_21_noInter>0.1)/length(Pval_21_noInter) #87% at 0.1 level

### Plot 4 panel figure with (a) P-values 1->2, (b) P-values 2->1 interactions and no interactions
### Other plot with lag order - should I overlay lines??
### Or should I do biplots? 
#plot(Pval_12_inter,Pval_21_inter,xlim=c(0,1),ylim=c(0,1))
#plot(Pval_12_noInter,Pval_21_noInter,xlim=c(0,1),ylim=c(0,1))
### Second hypothesis consistent with the null... 
DataCompet_stochModel_inter = data.frame(1:ncond,lag_order_inter_GC,Pval_12_inter_GC,Pval_21_inter_GC,index_1cause2_inter_GC,index_2cause1_inter_GC,lag_order_inter_CCM_predictx,Pval_12_inter_CCM,lag_order_inter_CCM_predicty,Pval_21_inter_CCM,index_1cause2_inter_CCM,index_2cause1_inter_CCM)

DataCompet_stochModel_noInter = data.frame(1:ncond,lag_order_noInter_GC,Pval_12_noInter_GC,Pval_21_noInter_GC,index_1cause2_noInter_GC,index_2cause1_noInter_GC,lag_order_noInter_CCM_predictx,Pval_12_noInter_CCM,lag_order_noInter_CCM_predicty,Pval_21_noInter_CCM,index_1cause2_noInter_CCM,index_2cause1_noInter_CCM)
#Write down results
write.csv(DataCompet_stochModel_inter,file="results/DataCompet_stochModel_inter.csv")
write.csv(DataCompet_stochModel_noInter,file="results/DataCompet_stochModel_noInter.csv")


pdf("stochasticcompet_ccm.pdf",width=7,height=5)
tmp=libsizes
plot(0,0,xlim=c(1,30),ylim=c(0,1.1),t="n",xlab="Library size",ylab="rho")
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
legend(x = "right", legend = c("sp1_xmap_sp2, inter", "sp2_xmap_sp1, inter","sp1_xmap_sp2, no inter", "sp2_xmap_sp1, no inter"), col = c("red","blue",rgb(255/256,165/256,0,alpha),rgb(155/256,79/256,150/256,alpha)), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
dev.off()
