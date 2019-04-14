### CP April 2019, based on FB's previous work
### Compares GC and CCM performance on the 2-species-and-a-driver model
### NOTE I think we should do it again, using the surrogates for the p-values

rm(list=ls())
graphics.off()

library(rEDM)
library(vars)


### Initialize
ncond=500

libsizes = seq(5, 100, by = 5)
lm=length(libsizes)
numsamples=100
set.seed(42)
tmax=300
seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates

Y_with=array(1,dim=c(tmax,2,ncond))
Y_without=array(1,dim=c(tmax,2,ncond))
y1=array(1,dim=c(tmax,ncond))

# 3 because we compare : 1-sp1 and temp, 2-sp1 and sp2, 3-sp2 and temp
Pval_12_inter_GC=Pval_21_inter_GC=matrix(NA,ncond,3)
Pval_12_inter_CCM=Pval_21_inter_CCM=matrix(NA,ncond,3)

RhoLMax_12_inter=RhoLMax_21_inter=matrix(NA,ncond,3)

index_1cause2_inter_GC=index_2cause1_inter_GC=matrix(NA,ncond,3)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=matrix(NA,ncond,3)

lag_order_inter_GC=matrix(NA,ncond,3)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=matrix(NA,ncond,3)

end_id=c("sp1temp","sp1sp2","sp2temp")

for(kcond in 1:ncond){
Y_with[1,1,kcond]=abs(rnorm(1,1,1))
Y_with[1,2,kcond]=abs(rnorm(1,1,1))
Y_without[1,1,kcond]=Y_with[1,1,kcond]
Y_without[1,2,kcond]=Y_with[1,2,kcond]

### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1[,kcond]<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y_with[t+1,1,kcond] = Y_with[t,1,kcond]*exp(3+0.5*y1[t,kcond] - 4*Y_with[t,1,kcond]-2*Y_with[t,2,kcond] + rnorm(1,0,0.1))
  Y_with[t+1,2,kcond] = Y_with[t,2,kcond]*exp(2.1+0.5*y1[t,kcond] -0.31*Y_with[t,1,kcond]-3.1*Y_with[t,2,kcond] + rnorm(1,0,0.1))
}

for (t in 1:(tmax-1)){
  Y_without[t+1,1,kcond] = Y_without[t,1,kcond]*exp(3+0.5*y1[t,kcond] - 4*Y_without[t,1,kcond]-0*Y_without[t,2,kcond] + rnorm(1,0,0.1))
  Y_without[t+1,2,kcond] = Y_without[t,2,kcond]*exp(2.1+0.5*y1[t,kcond] -0*Y_without[t,1,kcond]-3.1*Y_without[t,2,kcond] + rnorm(1,0,0.1))
}

}

for (i in 1:2){
for (kcond in 1:ncond){
	print(kcond)
	y1_tmp=y1[,kcond]-mean(y1[,kcond])
	if(i==1){
		y_with=log(Y_with[,,kcond])
		y_with[,1]=y_with[,1]-mean(y_with[,1])
		y_with[,2]=y_with[,2]-mean(y_with[,2])
		species2_species1_temp=data.frame(1:tmax,y_with,y1_tmp)
	}else{
		y_without=log(Y_without[,,kcond])
		y_without[,1]=y_without[,1]-mean(y_without[,1])
		y_without[,2]=y_without[,2]-mean(y_without[,2])
		species2_species1_temp=data.frame(1:tmax,y_without,y1_tmp)

	}

names(species2_species1_temp)=c("time","species1","species2","temp")

##########SP1 and TEMP
#### GC ####

varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,"species1"],species2_species1_temp[,"temp"])), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[kcond,1] <- varcompet$p

gxy = grangertest(species2_species1_temp[,"species1"],species2_species1_temp[,"temp"],order = lag_order_inter_GC[kcond]) #x causes y 
gyx = grangertest(species2_species1_temp[,"temp"],species2_species1_temp[,"species1"],order = lag_order_inter_GC[kcond]) #y causes x

Pval_12_inter_GC[kcond,1]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[kcond,1]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[kcond,1]<0.1)
{index_1cause2_inter_GC[kcond,1]=1} else {index_1cause2_inter_GC[kcond,1]=0}

if (Pval_21_inter_GC[kcond,1]<0.1)
{index_2cause1_inter_GC[kcond,1]=1} else {index_2cause1_inter_GC[kcond,1]=0}


### CCM ###
simplex_output_predictx = simplex(species2_species1_temp[,"species1"],E=1:10)
 lag_order_inter_CCM_predictx[kcond,1] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp[,"temp"],E=1:10)
 lag_order_inter_CCM_predicty[kcond,1] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_temp <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predictx[kcond,1], lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predicty[kcond,1], lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t <- ccm_means(species1_xmap_temp)
t_xmap_a <- ccm_means(temp_xmap_species1)

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[kcond,1]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond,1]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond,1]=a_xmap_t$rho[a_xmap_t$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond,1]=t_xmap_a$rho[t_xmap_a$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond,1]<0.1)&(RhoLMax_12_inter[kcond,1]>0.1))
{index_1cause2_inter_CCM[kcond,1]=1} else {index_1cause2_inter_CCM[kcond,1]=0}

if ((Pval_21_inter_CCM[kcond,1]<0.1)&(RhoLMax_21_inter[kcond,1]>0.1))
{index_2cause1_inter_CCM[kcond,1]=1} else {index_2cause1_inter_CCM[kcond,1]=0}



############SP1 and SP2
### GC ####
varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,"species1"],species2_species1_temp[,"species2"])), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[kcond,2] <- varcompet$p

gxy = grangertest(species2_species1_temp[,"species1"],species2_species1_temp[,"species2"],order = lag_order_inter_GC[kcond]) #x causes y 
gyx = grangertest(species2_species1_temp[,"species2"],species2_species1_temp[,"species1"],order = lag_order_inter_GC[kcond]) #y causes x

Pval_12_inter_GC[kcond,2]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[kcond,2]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[kcond,2]<0.1)
{index_1cause2_inter_GC[kcond,2]=1} else {index_1cause2_inter_GC[kcond,2]=0}

if (Pval_21_inter_GC[kcond,2]<0.1)
{index_2cause1_inter_GC[kcond,2]=1} else {index_2cause1_inter_GC[kcond,2]=0}


### CCM ###
 lag_order_inter_CCM_predictx[kcond,2] = lag_order_inter_CCM_predictx[kcond,1] #It's the same time series (sp1), so it's the same embedding

simplex_output_predicty = simplex(species2_species1_temp[,"species2"],E=1:10)
 lag_order_inter_CCM_predicty[kcond,2] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_species2 <- ccm(species2_species1_temp, E = lag_order_inter_CCM_predictx[kcond,2], lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp, E = lag_order_inter_CCM_predicty[kcond,2], lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t <- ccm_means(species1_xmap_species2)
t_xmap_a <- ccm_means(species2_xmap_species1)


### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[kcond,2]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond,2]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond,2]=a_xmap_t$rho[a_xmap_t$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond,2]=t_xmap_a$rho[t_xmap_a$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[kcond,2]<0.1)&(RhoLMax_12_inter[kcond,2]>0.1))
{index_1cause2_inter_CCM[kcond,2]=1} else {index_1cause2_inter_CCM[kcond,2]=0}

if ((Pval_21_inter_CCM[kcond,2]<0.1)&(RhoLMax_21_inter[kcond,2]>0.1))
{index_2cause1_inter_CCM[kcond,2]=1} else {index_2cause1_inter_CCM[kcond,2]=0}



##########SP2 and TEMP
### GC ####
varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,"species2"],species2_species1_temp[,"temp"])), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[kcond,3] <- varcompet$p

gxy = grangertest(species2_species1_temp[,"species2"],species2_species1_temp[,"temp"],order = lag_order_inter_GC[kcond]) #x causes y 
gyx = grangertest(species2_species1_temp[,"temp"],species2_species1_temp[,"species2"],order = lag_order_inter_GC[kcond]) #y causes x

Pval_12_inter_GC[kcond,3]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[kcond,3]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[kcond,3]<0.1)
{index_1cause2_inter_GC[kcond,3]=1} else {index_1cause2_inter_GC[kcond,3]=0}

if (Pval_21_inter_GC[kcond,3]<0.1)
{index_2cause1_inter_GC[kcond,3]=1} else {index_2cause1_inter_GC[kcond,3]=0}

### CCM ###
 lag_order_inter_CCM_predictx[kcond,3] = lag_order_inter_CCM_predicty[kcond,2] #species 2 was y in the 2nd model
 lag_order_inter_CCM_predicty[kcond,3] = lag_order_inter_CCM_predicty[kcond,1] #temp was y in the 1st model

species2_xmap_temp <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predictx[kcond,3], lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

temp_xmap_species2 <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predicty[kcond,3], lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t_means <- ccm_means(species2_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species2)
### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[kcond,3]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[kcond,3]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[kcond,3]=a_xmap_t$rho[a_xmap_t$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond,3]=t_xmap_a$rho[t_xmap_a$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2
                      
if ((Pval_12_inter_CCM[kcond,3]<0.1)&(RhoLMax_12_inter[kcond,3]>0.1))
{index_1cause2_inter_CCM[kcond,3]=1} else {index_1cause2_inter_CCM[kcond,3]=0}

if ((Pval_21_inter_CCM[kcond,3]<0.1)&(RhoLMax_21_inter[kcond,3]>0.1))
{index_2cause1_inter_CCM[kcond,3]=1} else {index_2cause1_inter_CCM[kcond,3]=0}

}
for(j in 1:3){
DataCompet_stochModel_inter = data.frame(1:ncond,lag_order_inter_GC[,j],Pval_12_inter_GC[,j],Pval_21_inter_GC[,j],lag_order_inter_CCM_predictx[kcond,j],lag_order_inter_CCM_predicty[kcond,j],index_1cause2_inter_GC[,j],index_2cause1_inter_GC[,j],Pval_12_inter_CCM[,j],Pval_21_inter_CCM[,j],index_1cause2_inter_CCM[,j],index_2cause1_inter_CCM[,j])

	if(i==1){
	#With interactions
	write.csv(DataCompet_stochModel_inter,file=paste("DataCompet_driver_inter",end_id[j],".csv",sep=""))

	}else{
	#Without interactions

	write.csv(DataCompet_stochModel_inter,file=paste("DataCompet_driver_noInter",end_id[j],".csv",sep=""))

	}
}
}
