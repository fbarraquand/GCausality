###CP April 2019, based on FB's previous work
### Simulates a chaotic 2-species deterministic model with and without interactions, for 500 different starting points.
## Compares results from CCM and Granger causality for each simulation

rm(list=ls())
graphics.off()
set.seed(42)

library(rEDM)
library(vars)

ncond<-500 #how many initial conditions we consider
tmax=800
tmin=501
inter_list=c(TRUE,FALSE)

pval_12_perlag=array(NA,dim=c(ncond,10,2)) 
pval_21_perlag=array(NA,dim=c(ncond,10,2)) 
tab_simu=array(NA,dim=c(ncond,28,2,2))


pdf("chaos_comparison_CCM.pdf",height=10,width=10)
par(mfcol=c(2,2),cex=1.25)
for(i_inter in 1:length(inter_list)){
inter=inter_list[i_inter]
if(i_inter==1){
mt="a)"
yl="Cross Map Skill"
}else{
mt="b)"
yl=""
}
par(mar=c(2,4,2,0.5))
plot(0,0,xlim=c(1,28),ylim=c(0,1),ylab=yl,xaxt="n",xlab="",t="n")
mtext(mt,side=2.5,line=2,at=1,las=2,cex=1.2)
#Initializing vectors 
Pval_12_inter_GC=Pval_21_inter_GC=Pval_12_noInter_GC=Pval_21_noInter_GC=rep(NA,ncond)
Pval_12_inter_CCM=Pval_21_inter_CCM=Pval_12_noInter_CCM=Pval_21_noInter_CCM=rep(NA,ncond)

RhoLMax_12_inter=RhoLMax_21_inter=RhoLMax_12_noInter=RhoLMax_21_noInter=rep(NA,ncond)

index_1cause2_inter_GC=index_2cause1_inter_GC=index_1cause2_noInter_GC=index_2cause1_noInter_GC=rep(NA,ncond)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=index_1cause2_noInter_CCM=index_2cause1_noInter_CCM=rep(NA,ncond)

lag_order_inter_GC=lag_order_noInter_GC=rep(NA,ncond)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=lag_order_noInter_CCM_predictx=lag_order_noInter_CCM_predicty=rep(NA,ncond)

for (kcond in 1:ncond){
print(kcond)
  ### Data simulation
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=runif(1,0.1,0.7)
Y[1]=runif(1,0.1,0.7)
for (t in 1:tmax)
{
if(inter){
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}else{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.0*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.0*X[t])
}

}
x=log(X[tmin:800])
y=log(Y[tmin:800])
# centering
x=x-mean(x)
y=y-mean(y)
z=cbind(x,y)

#####################################CCM#####################################
 ### determination of the lag order based on simplex
#Chose E
simplex_output_predictx = simplex(x,E=1:10)
 lag_order_inter_CCM_predictx[kcond] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_inter_CCM_predicty[kcond] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(tmin:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, nrow(species12)-11, by = 10)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==max(sp1_xmap_sp2$lib_size)]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==max(sp2_xmap_sp1$lib_size)]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

### Let's try to see those
sp1_xmap_sp2_means <- data.frame(ccm_means(sp1_xmap_sp2), sd.rho = with(sp1_xmap_sp2,
                                                                            tapply(rho, lib_size, sd)))
sp2_xmap_sp1_means <- data.frame(ccm_means(sp2_xmap_sp1), sd.rho = with(sp2_xmap_sp1,
                                                                            tapply(rho, lib_size, sd)))

tab_simu[kcond,,1,i_inter]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2,i_inter]=sp2_xmap_sp1_means$rho

if(kcond==1){
	lines(1:length(libsizes),tab_simu[kcond,,1,i_inter],col=rgb(1,0,0,1))
	lines(1:length(libsizes),sp1_xmap_sp2_means$rho+2*sp1_xmap_sp2_means$sd.rho,col=rgb(1,0,0,1),lty=2)
	lines(1:length(libsizes),sp1_xmap_sp2_means$rho-2*sp1_xmap_sp2_means$sd.rho,col=rgb(1,0,0,1),lty=2)
	lines(1:length(libsizes),tab_simu[kcond,,2,i_inter],col=rgb(0,0,1,1))
	lines(1:length(libsizes),sp2_xmap_sp1_means$rho+2*sp2_xmap_sp1_means$sd.rho,col=rgb(0,0,1,1),lty=2)
	lines(1:length(libsizes),sp2_xmap_sp1_means$rho-2*sp2_xmap_sp1_means$sd.rho,col=rgb(0,0,1,1),lty=2)
if(i_inter==1){
mt="c)"
yl="Cross Map Skill"
}else{
mt="d)"
yl=""
}
par(mar=c(4,4,0,0.5))
tmp_ax=seq(0,28,5)
tmp_ax[1]=1
plot(0,0,xlim=c(1,28),ylim=c(0,1),xlab="Library size",ylab=yl,t='n',xaxt="n")
axis(1,at=tmp_ax,lab=libsizes[tmp_ax])
mtext(mt,side=2.5,line=2,at=1,las=2,cex=1.2)

}

Pval_12_inter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==max(sp2_xmap_sp1_means$lib_size)] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==max(sp1_xmap_sp2_means$lib_size)] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond]<0.1)&(RhoLMax_12_inter[kcond]>0.1))
{index_1cause2_inter_CCM[kcond]=1} else {index_1cause2_inter_CCM[kcond]=0}

if ((Pval_21_inter_CCM[kcond]<0.1)&(RhoLMax_21_inter[kcond]>0.1))
{index_2cause1_inter_CCM[kcond]=1} else {index_2cause1_inter_CCM[kcond]=0}

#####################################Granger causalité#####################################
#Here, we compute Granger causality for the best lag, chosen with BIC
#We need to catch errors because some of the matrices can be singular in the estimation process
  tryCatch({

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

},
error=function(e){
	print(paste("Could not GC for",lag_order_inter_GC[kcond]))
})


#Here, we compute GC for all lags from 1 to 10, to check the proportion of false positive and negatives
#We need to catch errors because with high lags, time series are colinear, which leads the wald.test to crash (we could still re-write it, though)
tryCatch({
for(i in 1:10){
pval_12_perlag[kcond,i,i_inter]=grangertest(x,y,order = i)$`Pr(>F)`[2]
pval_21_perlag[kcond,i,i_inter]=grangertest(y,x,order = i)$`Pr(>F)`[2]
}
},
error=function(e){
	print(paste("error for lag",i))
})

}

DataCompet_chaos = data.frame(1:ncond,lag_order_inter_GC,Pval_12_inter_GC,Pval_21_inter_GC,index_1cause2_inter_GC,index_2cause1_inter_GC,lag_order_inter_CCM_predictx,Pval_12_inter_CCM,lag_order_inter_CCM_predicty,Pval_21_inter_CCM,index_1cause2_inter_CCM,index_2cause1_inter_CCM)

for(i in 2:ncond){
	lines(1:length(libsizes),tab_simu[i,,1,i_inter],col=rgb(1,0,0,0.1))
	lines(1:length(libsizes),tab_simu[i,,2,i_inter],col=rgb(0,0,1,0.1))
}
#Write down results
if(inter){
write.csv(DataCompet_chaos,file="results/DataCompet_chaos_withinter.csv")
}else{
write.csv(DataCompet_chaos,file="results/DataCompet_chaos_withoutinter.csv")
}

}
legend("topleft",c("2 causes 1","1 causes 2"),col=c("red","blue"),lty=1,bty="n")
dev.off()

#Plot the GC assessment of causality for all lags from 1 to 10, with and without interactions
pdf("GC_per_lag.pdf",height=5,width=10)
par(mfrow=c(1,2),cex=1.25,mar=c(4,4,1,0.5))

plot(1,0,t="n",xlab="Lag",ylab="% detected causality",ylim=c(0,1),xlim=c(1,11))
prop_ok=apply(pval_12_perlag[,,1],2,function(x) sum(x<0.1))/ncond
rect(1:length(prop_ok),rep(0,length(prop_ok)),seq(1.2,length(prop_ok)+0.2,1),prop_ok,col="red")
prop_ok=apply(pval_21_perlag[,,1],2,function(x) sum(x<0.1))/ncond
rect(seq(1.3,length(prop_ok)+0.3,1),rep(0,length(prop_ok)),seq(1.5,length(prop_ok)+0.5,1),prop_ok,col="blue")

plot(1,0,t="n",xlab="Lag",ylab="",ylim=c(0,1),xlim=c(1,11))
prop_ok=apply(pval_12_perlag[,,2],2,function(x) sum(x<0.1))/ncond
rect(1:length(prop_ok),rep(0,length(prop_ok)),seq(1.2,length(prop_ok)+0.2,1),prop_ok,col="red")
prop_ok=apply(pval_21_perlag[,,2],2,function(x) sum(x<0.1))/ncond
rect(seq(1.3,length(prop_ok)+0.3,1),rep(0,length(prop_ok)),seq(1.5,length(prop_ok)+0.5,1),prop_ok,col="blue")
legend("topright",c("1 -> 2","2 -> 1"),fill=c("red","blue"),bty="n")
dev.off()
