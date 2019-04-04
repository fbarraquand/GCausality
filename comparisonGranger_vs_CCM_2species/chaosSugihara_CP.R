rm(list=ls())
graphics.off()

library(rEDM)
library(vars)

ncond<-500 #how many initial conditions we consider
tmax=800
tmin=501
inter=FALSE


#Initializing vectors 
Pval_12_inter_GC=Pval_21_inter_GC=Pval_12_noInter_GC=Pval_21_noInter_GC=rep(NA,ncond)
Pval_12_inter_CCM=Pval_21_inter_CCM=Pval_12_noInter_CCM=Pval_21_noInter_CCM=rep(NA,ncond)

RhoLMax_12_inter=RhoLMax_21_inter=RhoLMax_12_noInter=RhoLMax_21_noInter=rep(NA,ncond)

index_1cause2_inter_GC=index_2cause1_inter_GC=index_1cause2_noInter_GC=index_2cause1_noInter_GC=rep(NA,ncond)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=index_1cause2_noInter_CCM=index_2cause1_noInter_CCM=rep(NA,ncond)

lag_order_inter_GC=lag_order_noInter_GC=rep(NA,ncond)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=lag_order_noInter_CCM_predictx=lag_order_noInter_CCM_predicty=rep(NA,ncond)

tab_simu=array(NA,dim=c(ncond,8,2))

#pdf("chaos_with_interactions_all_lines.pdf")
#plot(0,0,t="n",xlim=c(0,8),ylim=c(0,1.0),ylab="rho",xlab="theta")


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
smap_output_predictx = s_map(x,E=1:10,theta=2)
 lag_order_inter_CCM_predictx[kcond] = smap_output_predictx$E[which(smap_output_predictx$rho==max(smap_output_predictx$rho))]
#Chose theta
#plou=s_map(x,E= lag_order_inter_CCM_predictx[kcond])
#lines(plou$theta,plou$rho)


smap_output_predicty = s_map(y,E=1:10,theta=2)
 lag_order_inter_CCM_predicty[kcond] = smap_output_predicty$E[which(smap_output_predicty$rho==max(smap_output_predicty$rho))]
#Chose theta
#plou=s_map(y,E= lag_order_inter_CCM_predicty[kcond])
#lines(plou$theta,plou$rho,col="blue")
#CP : it seems that theta=2 is indeed good enough

 ### CCM Analysis 
species12=data.frame(tmin:tmax,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, 80, by = 10)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

### Let's try to see those
sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

tab_simu[kcond,,1]=sp1_xmap_sp2_means$rho
tab_simu[kcond,,2]=sp2_xmap_sp1_means$rho

#lines(sp1_xmap_sp2_means$lib_size, pmax(0, sp1_xmap_sp2_means$rho), col = "red")
#lines(sp2_xmap_sp1_means$lib_size, pmax(0, sp2_xmap_sp1_means$rho), col = "blue")


Pval_12_inter_CCM[kcond]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond]<0.1)&(RhoLMax_12_inter[kcond]>0.1))
{index_1cause2_inter_CCM[kcond]=1} else {index_1cause2_inter_CCM[kcond]=0}

if ((Pval_21_inter_CCM[kcond]<0.1)&(RhoLMax_21_inter[kcond]>0.1))
{index_2cause1_inter_CCM[kcond]=1} else {index_2cause1_inter_CCM[kcond]=0}

#####################################Granger causalité#####################################

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

})



#legend(x = "topleft", legend = c("sp1_xmap_sp2", "sp2_xmap_sp1"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
#dev.off()

}

DataCompet_chaos = data.frame(1:ncond,lag_order_inter_GC,Pval_12_inter_GC,Pval_21_inter_GC,index_1cause2_inter_GC,index_2cause1_inter_GC,lag_order_inter_CCM_predictx,Pval_12_inter_CCM,lag_order_inter_CCM_predicty,Pval_21_inter_CCM,index_1cause2_inter_CCM,index_2cause1_inter_CCM)

#Write down results
if(inter){
write.csv(DataCompet_chaos,file="results/DataCompet_chaos_withinter.csv")
}else{
write.csv(DataCompet_chaos,file="results/DataCompet_chaos_withoutinter.csv")
}

min1xmap2=apply(tab_simu[,,1],2,min)
min2xmap1=apply(tab_simu[,,2],2,min)
max1xmap2=apply(tab_simu[,,1],2,max)
max2xmap1=apply(tab_simu[,,2],2,max)

if(inter){
pdf("Chaos_shaded_area_withinter.pdf")
}else{
pdf("Chaos_shaded_area_withoutinter.pdf")
}
plot(0,0,xlim=c(1,8),ylim=c(0,1))
polygon(c(1:8, rev(1:8)), c(max1xmap2, rev(min1xmap2)),
     col=rgb(0, 0, 1,0.5),border = NA)
polygon(c(1:8, rev(1:8)), c(max2xmap1, rev(min2xmap1)),
     col=rgb(1, 0, 0,0.5),border=NA)
legend("topleft",c("1 xmap 2","2 xmap 1"),fill=c("blue","red"))
dev.off()

if(inter){
pdf("Chaos_lines_withinter.pdf")
}else{
pdf("Chaos_lines_withoutinter.pdf")
}
plot(0,0,xlim=c(1,8),ylim=c(0,1))
for(i in 1:ncond){
	lines(1:8,tab_simu[i,,1],col="blue")
	lines(1:8,tab_simu[i,,2],col="red")
}
legend("topleft",c("1 xmap 2","2 xmap 1"),col=c("blue","red"),lty=1,bty="n")
dev.off()
