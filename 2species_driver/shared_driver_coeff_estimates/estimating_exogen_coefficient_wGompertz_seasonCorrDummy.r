rm(list=ls())
graphics.off()
library("vars")

#set.seed(42)
set.seed(11)
tmax=800
seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates
ncond=500

noise_seq=c("white","autocorrelated","autocorrelatedANDseason")
### different types of noise -- the latter is used in the main text
### autocorrelated in Certain et al. 2018

growth_rate_seq=c("ref","low")
### ref.  = param set main text, low = low growth rate -> stable fixed point in absence of stochasticity.  

coeff_temperature=0.5 #0.5 in the original simulation

for(noise in noise_seq){
	if(grepl("season",noise)){
		correction_season_seq=c(0,24)
	}else{
		correction_season_seq=c(0)
	}
		for(growth_rate in growth_rate_seq){
	for(correction in correction_season_seq){	
#par(mfrow=c(2,1))
#plot(1:tmax,seasonality,type="o")


Y_with=array(1,dim=c(tmax,2,ncond))
Y_without=array(1,dim=c(tmax,2,ncond))
z_with=array(1,dim=c(tmax,2,ncond)) # for Gompertz growth
z_without=array(1,dim=c(tmax,2,ncond)) # for Gompertz growth
y1=array(1,dim=c(tmax,ncond))

for(kcond in 1:ncond){
Y_with[1,1,kcond]=abs(rnorm(1,1,1))
Y_with[1,2,kcond]=abs(rnorm(1,1,1))
Y_without[1,1,kcond]=Y_with[1,1,kcond]
Y_without[1,2,kcond]=Y_with[1,2,kcond]
z_with[1,1,kcond]=rnorm(1,0,0.1)
z_with[1,1,kcond]=rnorm(1,0,0.1)
z_without[1,1,kcond]=rnorm(1,0,0.1)
z_without[1,1,kcond]=rnorm(1,0,0.1)

### Environmental variables
if(noise=="white"){
	y1noise<-rnorm(tmax)
	y1[,kcond]<-y1noise
}else if(noise=="season"){
	y1[,kcond]<-seasonality

}else if(noise=="autocorrelated"){
	y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
	y1[,kcond]<-y1noise

}else if(noise=="autocorrelatedANDseason"){
	y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
	y1[,kcond]<-seasonality+y1noise
}

for (t in 1:(tmax-1)){
	if(growth_rate=="low"){
  Y_with[t+1,1,kcond] = Y_with[t,1,kcond]*exp(1+coeff_temperature*y1[t,kcond] - 4*Y_with[t,1,kcond]-2*Y_with[t,2,kcond] + rnorm(1,0,0.1))
  Y_with[t+1,2,kcond] = Y_with[t,2,kcond]*exp(0.5+coeff_temperature*y1[t,kcond] -0.31*Y_with[t,1,kcond]-3.1*Y_with[t,2,kcond] + rnorm(1,0,0.1))
	}else if(growth_rate=="ref"){
  Y_with[t+1,1,kcond] = Y_with[t,1,kcond]*exp(3+coeff_temperature*y1[t,kcond] - 4*Y_with[t,1,kcond]-2*Y_with[t,2,kcond] + rnorm(1,0,0.1))
  Y_with[t+1,2,kcond] = Y_with[t,2,kcond]*exp(2.1+coeff_temperature*y1[t,kcond] -0.31*Y_with[t,1,kcond]-3.1*Y_with[t,2,kcond] + rnorm(1,0,0.1))
	}

  z_with[t+1,1,kcond] = z_with[t,1,kcond] + 0.2 + coeff_temperature* y1[t,kcond] - 0.4*z_with[t,1,kcond] -0.1*z_with[t,2,kcond]+ rnorm(1,0,0.1)
  z_with[t+1,2,kcond] = z_with[t,2,kcond] + 0.1 +coeff_temperature*y1[t,kcond] -0.31*z_with[t,2,kcond] - 0.05*z_with[t,1,kcond]+ rnorm(1,0,0.1)
  
	if(growth_rate=="low"){
  Y_without[t+1,1,kcond] = Y_without[t,1,kcond]*exp(1+coeff_temperature*y1[t,kcond] - 4*Y_without[t,1,kcond]-0*Y_without[t,2,kcond] + rnorm(1,0,0.1))
  Y_without[t+1,2,kcond] = Y_without[t,2,kcond]*exp(0.5+coeff_temperature*y1[t,kcond] -0*Y_without[t,1,kcond]-3.1*Y_without[t,2,kcond] + rnorm(1,0,0.1))
	}else if(growth_rate=="ref"){
  Y_without[t+1,1,kcond] = Y_without[t,1,kcond]*exp(3+coeff_temperature*y1[t,kcond] - 4*Y_without[t,1,kcond]-0*Y_without[t,2,kcond] + rnorm(1,0,0.1))
  Y_without[t+1,2,kcond] = Y_without[t,2,kcond]*exp(2.1+coeff_temperature*y1[t,kcond] -0*Y_without[t,1,kcond]-3.1*Y_without[t,2,kcond] + rnorm(1,0,0.1))
}


  z_without[t+1,1,kcond] = z_without[t,1,kcond] + 0.2 +coeff_temperature* y1[t,kcond] - 0.4*z_without[t,1,kcond] + rnorm(1,0,0.1)
  z_without[t+1,2,kcond] = z_without[t,2,kcond] + 0.1 +coeff_temperature*y1[t,kcond] -0.31*z_without[t,2,kcond] + rnorm(1,0,0.1)
  }

}

#plot(1:tmax,y1noise,type="o")

#par(mfrow=c(2,1))

#plot(1:tmax,Y_with[,1,1],type="o",col="blue")
#lines(1:tmax,Y_with[,2,1],type="o",col="red")

#plot(1:tmax,Y_without[,1,1],type="o",col="blue")
#lines(1:tmax,Y_without[,2,1],type="o",col="red")

#plot(1:tmax,z_with[,1,1],type="o",col="blue",ylim=c(-3,3))
#lines(1:tmax,z_with[,2,1],type="o",col="red")

#plot(1:tmax,z_without[,1,1],type="o",col="blue",ylim=c(-3,3))
#lines(1:tmax,z_without[,2,1],type="o",col="red")

tmp_coeff1=matrix(NA,2,ncond)
tmp_coeff2=matrix(NA,2,ncond)
inter_coeff12=matrix(NA,2,ncond)
inter_coeff21=matrix(NA,2,ncond)

#for (type_inter in 1:4){ ### Considers Gompertz dynamics as well -- for code-checking purposes. 
for (type_inter in 1:2){
for (kcond in 1:ncond){
        print(kcond)
        #y1_tmp=y1[501:800,kcond]-mean(y1[501:800,kcond])
        y1_tmp=y1[500:799,kcond]-mean(y1[500:799,kcond]) #because VAR() uses y_t ~ y_{t-1} + ... + y_{t-p} + u_t and not u_{t-1}
        y1_tmp=y1_tmp-mean(y1_tmp)
        if(type_inter==1){
                y_with=log(Y_with[501:800,,kcond])
                y_with[,1]=y_with[,1]-mean(y_with[,1])
                y_with[,2]=y_with[,2]-mean(y_with[,2])
                species2_species1_temp=data.frame(1:nrow(y_with),y_with,y1_tmp)
        
        }else if(type_inter==2){        
                y_without=log(Y_without[501:800,,kcond])
                y_without[,1]=y_without[,1]-mean(y_without[,1])
                y_without[,2]=y_without[,2]-mean(y_without[,2])
                species2_species1_temp=data.frame(1:nrow(y_without),y_without,y1_tmp)
                
        }else if(type_inter==3){        
                z_w=z_with[501:800,,kcond]
                z_w[,1]=z_w[,1]-mean(z_w[,1])
                z_w[,2]=z_w[,2]-mean(z_w[,2])
              species2_species1_temp=data.frame(1:nrow(z_w),z_w,y1_tmp)
        
        }else{
                
                z_wout=z_without[501:800,,kcond]
                z_wout[,1]=z_wout[,1]-mean(z_wout[,1])
                z_wout[,2]=z_wout[,2]-mean(z_wout[,2])
                species2_species1_temp=data.frame(1:nrow(z_wout),z_wout,y1_tmp)
        }
names(species2_species1_temp)=c("time","species1","species2","temp")

##Now, we perform the analysis

name_x="species1"
name_y="species2"
name_exo="temp"

ss=correction
if(correction==0){
	ss=NULL
}
varcompet<-VAR(y=data.frame(cbind(species2_species1_temp[,name_x],species2_species1_temp[,name_y])), type="none",lag.max=20,ic="SC",exogen=species2_species1_temp[,name_exo],season=ss)
summ_varcompet=summary(varcompet)
effect_temp1=summ_varcompet$varresult$X1$coefficients["exo1",1]
effect_temp2=summ_varcompet$varresult$X2$coefficients["exo1",1]
Pval_temp1=summ_varcompet$varresult$X1$coefficients["exo1","Pr(>|t|)"]
Pval_temp2=summ_varcompet$varresult$X2$coefficients["exo1","Pr(>|t|)"]
	
effect_12=summ_varcompet$varresult$X2$coefficients["X1.l1",1]
effect_21=summ_varcompet$varresult$X1$coefficients["X2.l1",1]
# in case we want to examine bias on those as well

tmp_coeff1[type_inter,kcond]=effect_temp1
tmp_coeff2[type_inter,kcond]=effect_temp2

inter_coeff12[type_inter,kcond]=effect_12
inter_coeff21[type_inter,kcond]=effect_21
}
}

pdf(paste("GR",growth_rate,"forcing",noise,"tempcoeff",coeff_temperature,"correctionforseasonality",correction,"temperatureeffect.pdf",sep="_"))
par(mfrow=c(2,2))
hist(tmp_coeff1[1,],xlab="Temp->1 WI",main=paste("GR=",growth_rate,"forcing = ",noise,"coef = ",coeff_temperature,"corr = ",correction),cex.main=0.7,xlim=range(c(c(tmp_coeff1[1,]),coeff_temperature)))
abline(v=coeff_temperature,col="red",lwd=1.2)
hist(tmp_coeff1[2,],xlab="Temp->1 WoutI",main="",xlim=range(c(c(tmp_coeff1[2,]),coeff_temperature)))
abline(v=coeff_temperature,col="red",lwd=1.2)
hist(tmp_coeff2[1,],xlab="Temp->2 WI",main="",xlim=range(c(c(tmp_coeff1[1,]),coeff_temperature)))
abline(v=coeff_temperature,col="red",lwd=1.2)
hist(tmp_coeff2[2,],xlab="Temp->2 WoutI",main="",xlim=range(c(c(tmp_coeff1[1,]),coeff_temperature)))
abline(v=coeff_temperature,col="red",lwd=1.2)
dev.off()

pdf(paste("GR",growth_rate,"forcing",noise,"tempcoeff",coeff_temperature,"correctionforseasonality",correction,"interactions_lag1.pdf",sep="_"))
par(mfrow=c(2,2))
inter = 0
hist(inter_coeff12[1,],xlab="1->2 WI",main=paste("GR=",growth_rate,noise,"coef",coeff_temperature,"Seas",correction),cex.main=0.7,xlim=range(c(c(inter_coeff12[1,]),inter)))
hist(inter_coeff12[2,],xlab="1->2 WoutI",main="",xlim=range(c(c(inter_coeff12[2,]),inter)))
hist(inter_coeff21[1,],xlab="2->1 WI",main="",xlim=range(c(c(inter_coeff21[1,]),inter)))
hist(inter_coeff21[2,],xlab="2->1 WoutI",main="",xlim=range(c(c(inter_coeff21[2,]),inter)))
dev.off()


}
}
}
