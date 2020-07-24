#######################################################################################################################
### Compares GC and CCM performance on the 2-species-and-a-driver model                                          ######
### From 2competitors_and_climate_withSurrogates.R, FBarraquand 2018                                             ######
### CP April 2019                                                                                                ######
#######################################################################################################################

### 16/04/2019 : Changing from twin surrogate data to seasonal surrogate data because a) the former tended to crash, b) apart from phase locking, it's not handling seasonality as well (at least, as clearly), as seasonal
### Also corrected a bug on Rho (Rho12=rho(2xmap1))
### Also modified libsizes to remove burn-in AND use 100 for rho max
### July 2019 factorized the code 
### 17/06/2020 corrected variable names and commented code + added a comparison for Wald vs F-test + added a GC-test for temperature

rm(list=ls())
graphics.off()

library(vars)
library(rEDM)
library(VARtests)

standard_GC=function(species2_species1_temp,name_x,name_y,name_exo=NA){
if(is.na(name_exo)){ #Pairwise GC, no exogen variable
	varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,name_x],species2_species1_temp[,name_y])), type="none",lag.max=20,ic="SC")
	vartot=VARfit(cbind(species2_species1_temp[,name_x],species2_species1_temp[,name_y]),p=varcompet$p)
	ar_x=VARfit(species2_species1_temp[,name_x],p=varcompet$p)
	ar_y=VARfit(species2_species1_temp[,name_y],p=varcompet$p)
	effect_temp1=NA
	effect_temp2=NA
	Pval_temp1=NA
	Pval_temp2=NA
}else{ #Conditional GC, exogen variable
	varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,name_x],species2_species1_temp[,name_y])), type="none",lag.max=20,ic="SC",exogen=species2_species1_temp[,name_exo])
	summ_varcompet=summary(varcompet)
	vartot=VARfit(cbind(species2_species1_temp[,name_x],species2_species1_temp[,name_y]),p=varcompet$p,exogen=species2_species1_temp[,name_exo])
	ar_x=VARfit(species2_species1_temp[,name_x],p=varcompet$p,exogen=species2_species1_temp[,name_exo])
	ar_y=VARfit(species2_species1_temp[,name_y],p=varcompet$p,exogen=species2_species1_temp[,name_exo])
	effect_temp1=summ_varcompet$varresult$X1$coefficients["exo1",1]
	effect_temp2=summ_varcompet$varresult$X2$coefficients["exo1",1]
	Pval_temp1=summ_varcompet$varresult$X1$coefficients["exo1","Pr(>|t|)"]
	Pval_temp2=summ_varcompet$varresult$X2$coefficients["exo1","Pr(>|t|)"]
}
lag_order_inter_GC<- varcompet$p

#Let's compute log ratio #there is only 10^-1 differences for the sum of squares of residual between VARfit and VAR
log_12_inter=log(sum((ar_y$resid[,1])^2,na.rm=T)/sum((vartot$resid[,2])^2,na.rm=T))
log_21_inter=log(sum((ar_x$resid[,1])^2,na.rm=T)/sum((vartot$resid[,1])^2,na.rm=T))

n1=names(varcompet$varresult$X1$coefficients)
effect_21_inter=sum(abs(varcompet$varresult$X1$coefficients[grep("X2",n1)]))
n2=names(varcompet$varresult$X2$coefficients)
effect_12_inter=sum(abs(varcompet$varresult$X2$coefficients[grep("X1",n2)]))

gxy=causality(varcompet,cause="X1") #Ftest
gyx=causality(varcompet,cause="X2") #Ftest

Pval_12_inter_Ftest=gxy$Granger$p.value 
Pval_21_inter_Ftest=gyx$Granger$p.value 

if(is.na(name_exo)){ #Pairwise GC, no exogen variable
gxy = grangertest(species2_species1_temp[,name_x],species2_species1_temp[,name_y],order = lag_order_inter_GC) #Wald test
gyx = grangertest(species2_species1_temp[,name_y],species2_species1_temp[,name_x],order = lag_order_inter_GC) #Wald test
Pval_12_inter_Wald=gxy$`Pr(>F)`[2]
Pval_21_inter_Wald=gyx$`Pr(>F)`[2]
}else{ #We don't want to compare Wald and Ftest for the conditional case (grangertest does not allow that)
Pval_12_inter_Wald=NA
Pval_21_inter_Wald=NA
}


return(list(lag_order_inter_GC,log_12_inter,log_21_inter,effect_12_inter,effect_21_inter,Pval_12_inter_Ftest,Pval_21_inter_Ftest,Pval_12_inter_Wald,Pval_21_inter_Wald,effect_temp1,effect_temp2,Pval_temp1,Pval_temp2))
}

standard_CCM=function(species2_species1_temp,name_x,name_y){
#libsizes = seq(5, 100, by = 5)
libsizes = seq(10, nrow(species2_species1_temp)-11, by = 10)
lm=length(libsizes)
num_surr=100
numsamples=100

simplex_output_predictx = simplex(species2_species1_temp[,name_x],E=1:10)
 lag_order_inter_CCM_predictx = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))] #lagx is the embedding for species1

simplex_output_predicty = simplex(species2_species1_temp[,name_y],E=1:10)
 lag_order_inter_CCM_predicty= simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))] #lagy is the embedding for temperature

species1_xmap_temp <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predictx, lib_column = name_x,
                      target_column = name_y, lib_size=libsizes,num_samples=numsamples)

temp_xmap_species1 <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predicty, lib_column = name_y, target_column = name_x,
                      lib_sizes = libsizes,num_samples=numsamples)

a_xmap_t <- ccm_means(species1_xmap_temp) #1 xmap 2
t_xmap_a <- ccm_means(temp_xmap_species1) #2 xmap 1

### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- a_xmap_t$rho[a_xmap_t$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- a_xmap_t$rho[a_xmap_t$lib_size ==max(a_xmap_t$lib_size)]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- t_xmap_a$rho[t_xmap_a$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- t_xmap_a$rho[t_xmap_a$lib_size ==max(t_xmap_a$lib_size)]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

a_xmap_t <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predictx, lib_column = name_x,
                      target_column = name_y, lib_size=max(a_xmap_t$lib_size),num_samples=1,replace=F)

t_xmap_a <- ccm(species2_species1_temp, E =  lag_order_inter_CCM_predicty, lib_column = name_y, target_column = name_x,
                      lib_sizes =max(t_xmap_a$lib_size),replace=F,num_samples=1)

RhoLMax_21_inter=a_xmap_t$rho[a_xmap_t$lib_size==max(a_xmap_t$lib_size)] # 1 causes 2 if 2 xmap 1
RhoLMax_12_inter=t_xmap_a$rho[t_xmap_a$lib_size==max(t_xmap_a$lib_size)] # 2 causes 1 if 1 xmap 2

### pvalue from surrogates seasonal
surr_temp_twin <- make_surrogate_data(species2_species1_temp[,name_y], method = "seasonal", num_surr = num_surr,T_period=24)
surr_species1_twin <- make_surrogate_data(species2_species1_temp[,name_x], method = "seasonal", num_surr = num_surr,T_period=24)
rho_surr1_twin<- list(temp =  matrix(NA,nrow=num_surr,ncol=1), species1 =  matrix(NA,nrow=num_surr,ncol=1))

for (i in 1:num_surr) {
  rho_surr1_twin$temp[i,1] <- ccm(cbind(species2_species1_temp[,name_x], surr_temp_twin[,i]), E = lag_order_inter_CCM_predictx, lib_column = 1, target_column = 2, lib_sizes = max(a_xmap_t$lib_size),replace = FALSE,num_samples=1)$rho #temp is the cause, so it is the target, and the "surrogated" time series, and the embedding is the one for species 1, ie lagx
  rho_surr1_twin$species1[i,1] <- ccm(cbind(surr_species1_twin[,i], species2_species1_temp[,name_y]), E = lag_order_inter_CCM_predicty , lib_column = 2, target_column = 1, lib_sizes = max(t_xmap_a$lib_size),num_samples=1,replace = FALSE)$rho #sp1 is the cause, so it is the target and the "surrogated" TS, and the embedding is the one for temperature, ie lagy
}

Pval_21_inter_CCM_surr_seasonal=(sum(RhoLMax_21_inter<=rho_surr1_twin$temp)+1) /(num_surr+1)
Pval_12_inter_CCM_surr_seasonal=(sum(RhoLMax_12_inter<=rho_surr1_twin$species1)+1) /(num_surr+1)

### pvalue from surrogates
rho_surr1_twin<- list(temp =  matrix(NA,nrow=num_surr,ncol=1), species1 =  matrix(NA,nrow=num_surr,ncol=1))
for (i in 1:num_surr) {
  rho_surr1_twin$temp[i,1] <- ccm(cbind(species2_species1_temp[,name_x], sample(species2_species1_temp[,name_y])), E = lag_order_inter_CCM_predictx, lib_column = 1, target_column = 2, lib_sizes = max(a_xmap_t$lib_size),replace = FALSE,num_samples=1)$rho #temp is the cause, so it is the target, and the "surrogated" time series, and the embedding is the one for species 1, ie lagx
  rho_surr1_twin$species1[i,1] <- ccm(cbind(sample(species2_species1_temp[,name_x]), species2_species1_temp[,name_y]), E = lag_order_inter_CCM_predicty , lib_column = 2, target_column = 1, lib_sizes = max(t_xmap_a$lib_size),num_samples=1,replace = FALSE)$rho #sp1 is the cause, so it is the target and the "surrogated" TS, and the embedding is the one for temperature, ie lagy
}

Pval_21_inter_CCM_surr_sample=(sum(RhoLMax_21_inter<rho_surr1_twin$temp)+1) /(num_surr+1)
Pval_12_inter_CCM_surr_sample=(sum(RhoLMax_12_inter<rho_surr1_twin$species1)+1) /(num_surr+1)

return(list(lag_order_inter_CCM_predictx,lag_order_inter_CCM_predicty,RhoLMax_12_inter,RhoLMax_21_inter,Pval_12_inter_CCM,Pval_21_inter_CCM,Pval_12_inter_CCM_surr_seasonal,Pval_21_inter_CCM_surr_seasonal,Pval_12_inter_CCM_surr_sample,Pval_21_inter_CCM_surr_sample))
}

######################################### Let's start computing
### Initialize
ncond=500

set.seed(42)
#tmax=300
tmax=800
seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates

Y_with=array(1,dim=c(tmax,2,ncond))
Y_without=array(1,dim=c(tmax,2,ncond))
y1=array(1,dim=c(tmax,ncond))

Pval_12_inter_GC_Ftest=Pval_21_inter_GC_Ftest=matrix(NA,ncond,2) #With and without exogen variable
Pval_12_inter_GC_Wald=Pval_21_inter_GC_Wald=matrix(NA,ncond,2) #With and without exogen variable
Pval_temp1=Pval_temp2=matrix(NA,ncond,2) #Pvalue for the effect of the temperature (/exogen variable) on sp1 or sp2
Pval_12_inter_CCM=Pval_21_inter_CCM=matrix(NA,ncond,3) #1-sp1 and temp, 2-sp1 and sp2, 3-sp2 and temp
Pval_12_inter_CCM_surr_seasonal=Pval_21_inter_CCM_surr_seasonal=matrix(NA,ncond,3)
Pval_12_inter_CCM_surr_sample=Pval_21_inter_CCM_surr_sample=matrix(NA,ncond,3)

log_12_inter=log_21_inter=matrix(NA,ncond,2)
effect_12_inter=effect_21_inter=matrix(NA,ncond,2)
effect_temp1=effect_temp2=matrix(NA,ncond,2) #Effect of the temperature
RhoLMax_12_inter=RhoLMax_21_inter=matrix(NA,ncond,3)

lag_order_inter_GC=matrix(NA,ncond,2)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=matrix(NA,ncond,3)

end_id=c("sp1temp","sp1sp2","sp2temp")


names_GC=c("Time","lag_order_inter_GC_exo","Pval_12_inter_GC_exo","Pval_21_inter_GC_exo","effect_12_inter_exo","effect_21_inter_exo","log_12_inter_exo","log_21_inter_exo","lag_order_inter_GC_noexo","Pval_12_inter_GC_noexo_Ftest","Pval_21_inter_GC_no_exo_Ftest","Pval_12_inter_GC_noexo_Wald","Pval_21_inter_GC_no_exo_Wald","effect_12_inter_noexo","effect_21_inter_noexo","log_12_inter_noexo","log_21_inter_noexo","effect_exo1","effect_exo2","Pval_exo1","Pval_exo2")
for(j in 1:3){
write(names_GC,file=paste("results/DataCompet_driver_inter",end_id[j],"factorized_GC_otf_with_F_Wald_test_with800.csv",sep=""),sep=",",append=F,ncolumns=length(names_GC))
}

names_CCM=c("Time","lag_order_inter_CCM_predictx","lag_order_inter_CCM_predicty","Pval_12_inter_CCM","Pval_21_inter_CCM","Rho_12","Rho_21","Pval_12_inter_CCM_surr_season","Pval_21_inter_CCM_surr_season","Pval_12_inter_CCM_surr_sample","Pval_21_inter_CCM_surr_sample")
for(j in 1:3){
write(names_CCM,file=paste("results/DataCompet_driver_inter",end_id[j],"factorized_CCM_otf_test_with800.csv",sep=""),sep=",",append=F,ncolumns=length(names_CCM))
write(names_CCM,file=paste("results/DataCompet_driver_noInter",end_id[j],"factorized_CCM_otf_test_with800.csv",sep=""),sep=",",append=F,ncolumns=length(names_CCM))
}

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

for (type_inter in 1:2){
for (kcond in 1:ncond){
        print(kcond)
 ##       y1_tmp=y1[201:300,kcond]-mean(y1[201:300,kcond]) CP Saw a pb here
        y1_tmp=y1[501:800,kcond]-mean(y1[501:800,kcond])
        if(type_inter==1){
                ##y_with=log(Y_with[201:300,,kcond]) CP Saw a pb here
                y_with=log(Y_with[501:800,,kcond])
                y_with[,1]=y_with[,1]-mean(y_with[,1])
                y_with[,2]=y_with[,2]-mean(y_with[,2])
                #species2_species1_temp=data.frame(1:100,y_with,y1_tmp)
                species2_species1_temp=data.frame(1:nrow(y_with),y_with,y1_tmp)
        }else{
                ##y_without=log(Y_without[201:300,,kcond]) CP Saw a pb here
                y_without=log(Y_without[501:800,,kcond])
                y_without[,1]=y_without[,1]-mean(y_without[,1])
                y_without[,2]=y_without[,2]-mean(y_without[,2])
                #species2_species1_temp=data.frame(1:100,y_without,y1_tmp)
                species2_species1_temp=data.frame(1:nrow(y_without),y_without,y1_tmp)

        }

names(species2_species1_temp)=c("time","species1","species2","temp")

######## Granger cause
val_id=2 #sp1sp2
#With exogen variable
plou=standard_GC(species2_species1_temp,"species1","species2",name_exo="temp")
lag_order_inter_GC[kcond,1]=plou[[1]]
log_12_inter[kcond,1]=plou[[2]]
log_21_inter[kcond,1]=plou[[3]]
effect_12_inter[kcond,1]=plou[[4]]
effect_21_inter[kcond,1]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,1]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,1]=plou[[7]]
effect_temp1[kcond,1]=plou[[10]]
effect_temp2[kcond,1]=plou[[11]]
Pval_temp1[kcond,1]=plou[[12]]
Pval_temp2[kcond,1]=plou[[13]]

#Without exogen variables
plou=standard_GC(species2_species1_temp,"species1","species2",name_exo=NA)
lag_order_inter_GC[kcond,2]=plou[[1]]
log_12_inter[kcond,2]=plou[[2]]
log_21_inter[kcond,2]=plou[[3]]
effect_12_inter[kcond,2]=plou[[4]]
effect_21_inter[kcond,2]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,2]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,2]=plou[[7]]
Pval_12_inter_GC_Wald[kcond,2]=plou[[8]]
Pval_21_inter_GC_Wald[kcond,2]=plou[[9]]

write(c(kcond,lag_order_inter_GC[kcond,1],Pval_12_inter_GC_Ftest[kcond,1],Pval_21_inter_GC_Ftest[kcond,1],effect_12_inter[kcond,1],effect_21_inter[kcond,1],log_12_inter[kcond,1],log_21_inter[kcond,1],lag_order_inter_GC[kcond,2],Pval_12_inter_GC_Ftest[kcond,2],Pval_21_inter_GC_Ftest[kcond,2],Pval_12_inter_GC_Wald[kcond,2],Pval_21_inter_GC_Wald[kcond,2],effect_12_inter[kcond,2],effect_21_inter[kcond,2],log_12_inter[kcond,2],log_21_inter[kcond,2],effect_temp1[kcond,1],effect_temp2[kcond,1],Pval_temp1[kcond,1],Pval_temp2[kcond,1]),file=paste("results/DataCompet_driver_inter",end_id[val_id],"factorized_GC_otf_with_F_Wald_test_with800.csv",sep=""),sep=",",append=T,ncolumns=length(names_GC))

val_id=1 #sp1temp
#With exogen variable
plou=standard_GC(species2_species1_temp,"species1","temp",name_exo="species2")
lag_order_inter_GC[kcond,1]=plou[[1]]
log_12_inter[kcond,1]=plou[[2]]
log_21_inter[kcond,1]=plou[[3]]
effect_12_inter[kcond,1]=plou[[4]]
effect_21_inter[kcond,1]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,1]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,1]=plou[[7]]

#Without exogen variables
plou=standard_GC(species2_species1_temp,"species1","temp",name_exo=NA)
lag_order_inter_GC[kcond,2]=plou[[1]]
log_12_inter[kcond,2]=plou[[2]]
log_21_inter[kcond,2]=plou[[3]]
effect_12_inter[kcond,2]=plou[[4]]
effect_21_inter[kcond,2]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,2]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,2]=plou[[7]]
Pval_12_inter_GC_Wald[kcond,2]=plou[[8]]
Pval_21_inter_GC_Wald[kcond,2]=plou[[9]]

write(c(kcond,lag_order_inter_GC[kcond,1],Pval_12_inter_GC_Ftest[kcond,1],Pval_21_inter_GC_Ftest[kcond,1],effect_12_inter[kcond,1],effect_21_inter[kcond,1],log_12_inter[kcond,1],log_21_inter[kcond,1],lag_order_inter_GC[kcond,2],Pval_12_inter_GC_Ftest[kcond,2],Pval_21_inter_GC_Ftest[kcond,2],Pval_12_inter_GC_Wald[kcond,2],Pval_21_inter_GC_Wald[kcond,2],effect_12_inter[kcond,2],effect_21_inter[kcond,2],log_12_inter[kcond,2],log_21_inter[kcond,2],effect_temp1[kcond,1],effect_temp2[kcond,1],Pval_temp1[kcond,1],Pval_temp2[kcond,1]),file=paste("results/DataCompet_driver_inter",end_id[val_id],"factorized_GC_otf_with_F_Wald_test_with800.csv",sep=""),sep=",",append=T,ncolumns=length(names_GC))

val_id=3 #sp2temp
#With exogen variable
plou=standard_GC(species2_species1_temp,"species2","temp",name_exo="species1")
lag_order_inter_GC[kcond,1]=plou[[1]]
log_12_inter[kcond,1]=plou[[2]]
log_21_inter[kcond,1]=plou[[3]]
effect_12_inter[kcond,1]=plou[[4]]
effect_21_inter[kcond,1]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,1]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,1]=plou[[7]]

#Without exogen variables
plou=standard_GC(species2_species1_temp,"species2","temp",name_exo=NA)
lag_order_inter_GC[kcond,2]=plou[[1]]
log_12_inter[kcond,2]=plou[[2]]
log_21_inter[kcond,2]=plou[[3]]
effect_12_inter[kcond,2]=plou[[4]]
effect_21_inter[kcond,2]=plou[[5]]
Pval_12_inter_GC_Ftest[kcond,2]=plou[[6]]
Pval_21_inter_GC_Ftest[kcond,2]=plou[[7]]
Pval_12_inter_GC_Wald[kcond,2]=plou[[8]]
Pval_21_inter_GC_Wald[kcond,2]=plou[[9]]

write(c(kcond,lag_order_inter_GC[kcond,1],Pval_12_inter_GC_Ftest[kcond,1],Pval_21_inter_GC_Ftest[kcond,1],effect_12_inter[kcond,1],effect_21_inter[kcond,1],log_12_inter[kcond,1],log_21_inter[kcond,1],lag_order_inter_GC[kcond,2],Pval_12_inter_GC_Ftest[kcond,2],Pval_21_inter_GC_Ftest[kcond,2],Pval_12_inter_GC_Wald[kcond,2],Pval_21_inter_GC_Wald[kcond,2],effect_12_inter[kcond,2],effect_21_inter[kcond,2],log_12_inter[kcond,2],log_21_inter[kcond,2],effect_temp1[kcond,1],effect_temp2[kcond,1],Pval_temp1[kcond,1],Pval_temp2[kcond,1]),file=paste("results/DataCompet_driver_inter",end_id[val_id],"factorized_GC_otf_with_F_Wald_test_with800.csv",sep=""),sep=",",append=T,ncolumns=length(names_GC))


###############CCM 
#SP1 and temp
id_simu=1
plou=standard_CCM(species2_species1_temp,"species1","temp")
lag_order_inter_CCM_predictx[kcond,id_simu]=plou[[1]]
lag_order_inter_CCM_predicty[kcond,id_simu]=plou[[2]]
RhoLMax_12_inter[kcond,id_simu]=plou[[3]]
RhoLMax_21_inter[kcond,id_simu]=plou[[4]]
Pval_12_inter_CCM[kcond,id_simu]=plou[[5]]
Pval_21_inter_CCM[kcond,id_simu]=plou[[6]]
Pval_12_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[7]]
Pval_21_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[8]]
Pval_12_inter_CCM_surr_sample[kcond,id_simu]=plou[[9]]
Pval_21_inter_CCM_surr_sample[kcond,id_simu]=plou[[10]]

#SP1 and SP2
id_simu=2
plou=standard_CCM(species2_species1_temp,"species1","species2")
lag_order_inter_CCM_predictx[kcond,id_simu]=plou[[1]]
lag_order_inter_CCM_predicty[kcond,id_simu]=plou[[2]]
RhoLMax_12_inter[kcond,id_simu]=plou[[3]]
RhoLMax_21_inter[kcond,id_simu]=plou[[4]]
Pval_12_inter_CCM[kcond,id_simu]=plou[[5]]
Pval_21_inter_CCM[kcond,id_simu]=plou[[6]]
Pval_12_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[7]]
Pval_21_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[8]]
Pval_12_inter_CCM_surr_sample[kcond,id_simu]=plou[[9]]
Pval_21_inter_CCM_surr_sample[kcond,id_simu]=plou[[10]]

#SP2 and temp
id_simu=3
plou=standard_CCM(species2_species1_temp,"species2","temp")
lag_order_inter_CCM_predictx[kcond,id_simu]=plou[[1]]
lag_order_inter_CCM_predicty[kcond,id_simu]=plou[[2]]
RhoLMax_12_inter[kcond,id_simu]=plou[[3]]
RhoLMax_21_inter[kcond,id_simu]=plou[[4]]
Pval_12_inter_CCM[kcond,id_simu]=plou[[5]]
Pval_21_inter_CCM[kcond,id_simu]=plou[[6]]
Pval_12_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[7]]
Pval_21_inter_CCM_surr_seasonal[kcond,id_simu]=plou[[8]]
Pval_12_inter_CCM_surr_sample[kcond,id_simu]=plou[[9]]
Pval_21_inter_CCM_surr_sample[kcond,id_simu]=plou[[10]]

for(j in 1:3){
if(type_inter==1){
write(c(kcond,lag_order_inter_CCM_predictx[kcond,j],lag_order_inter_CCM_predicty[kcond,j],Pval_12_inter_CCM[kcond,j],Pval_21_inter_CCM[kcond,j],RhoLMax_12_inter[kcond,j],RhoLMax_21_inter[kcond,j],Pval_12_inter_CCM_surr_seasonal[kcond,j],Pval_21_inter_CCM_surr_seasonal[kcond,j],Pval_12_inter_CCM_surr_sample[kcond,j],Pval_21_inter_CCM_surr_sample[kcond,j]),file=paste("results/DataCompet_driver_inter",end_id[j],"factorized_CCM_otf_test_with800.csv",sep=""),sep=",",append=T,ncolumns=length(names_CCM))
}else{
write(c(kcond,lag_order_inter_CCM_predictx[kcond,j],lag_order_inter_CCM_predicty[kcond,j],Pval_12_inter_CCM[kcond,j],Pval_21_inter_CCM[kcond,j],RhoLMax_12_inter[kcond,j],RhoLMax_21_inter[kcond,j],Pval_12_inter_CCM_surr_seasonal[kcond,j],Pval_21_inter_CCM_surr_seasonal[kcond,j],Pval_12_inter_CCM_surr_sample[kcond,j],Pval_21_inter_CCM_surr_sample[kcond,j]),file=paste("results/DataCompet_driver_noInter",end_id[j],"factorized_CCM_otf_test_with800.csv",sep=""),sep=",",append=T,ncolumns=length(names_CCM))
}
}


}
}

