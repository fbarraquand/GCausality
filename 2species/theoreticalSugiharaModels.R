########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################

library("vars")

########################################################################################################################
########## Sugihara two species deterministic competition model -- interactions
########################################################################################################################

tmax=300
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=0.1
Y[1]=0.2
for (t in 1:tmax)
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
pdf(file="fig/CompetitionDeterministicNL2Species.pdf",width=12,height=6)
par(cex=1.5)
minz=min(z[151:201,])
maxz=max(z[151:201,])
for (i in 1:2){
  if (i==1){plot(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz),xlab="Time",ylab="ln(Density)")}
  lines(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz))
}
dev.off()

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
IC=VARselect(z,lag.max=20) ### Order selected is 7

crit=scale(t(IC$criteria))

pdf(file="fig/LagOrderSelection2Species_chaos_interaction.pdf",width=8,height=6)
par(cex=1.5,lwd=3)
plot(1:20,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-1.3,3),col='black')
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:20,crit[,i],type="o",col=col_vec[i])} 
legend(5,2.4,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18)
dev.off()

causality(varcompet,cause="x") #p-value < 2.2e-16
causality(varcompet,cause="y") #0.07526

########################################################################################################################
########## Sugihara two species deterministic competition model -- no interactions
########################################################################################################################

tmax=300
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=0.1
Y[1]=0.2
for (t in 1:tmax)
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
pdf(file="fig/CompetitionDeterministicNL2Species.pdf",width=12,height=6)
par(cex=1.5)
minz=min(z[151:201,])
maxz=max(z[151:201,])
for (i in 1:2){
  if (i==1){plot(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz),xlab="Time",ylab="ln(Density)")}
  lines(151:201,z[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(minz,maxz))
}
dev.off()

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
VARselect(z) ### Order selected is 7

causality(varcompet,cause="x") #p-value < 2.2e-16
causality(varcompet,cause="y") #0.07526


########################################################################################################################
####################### Repeating many times with many initial conditions           ####################################
########################################################################################################################

#Initializing vectors 
ncond<-500
Pval_12_inter_lag1=Pval_21_inter_lag1=Pval_12_inter=Pval_21_inter=Pval_12_noInter=Pval_21_noInter=Pval_12_noInter_lag1=Pval_21_noInter_lag1=lag_order_inter=lag_order_noInter=rep(NA,ncond)

########################################################################################################################
########## Sugihara two species deterministic competition model -- interactions
########################################################################################################################


for (kcond in 1:ncond){

tmax=800
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=runif(1,0.1,0.7)
Y[1]=runif(1,0.1,0.7)
for (t in 1:tmax)
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}
#X
#Y
x=log(X[501:800])
y=log(Y[501:800])
# centering
x=x-mean(x)
y=y-mean(y)
z=cbind(x,y)

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
lag_order_inter[kcond] <- varcompet$p
## VARselect(z) ### Order selected is 7
## causality(varcompet,cause="x") #p-value < 2.2e-16
## causality(varcompet,cause="y")

gxy = grangertest(x,y,order = lag_order_inter[kcond]) #x causes y 
gyx = grangertest(y,x,order = lag_order_inter[kcond]) #y causes x

Pval_12_inter[kcond]=gxy$`Pr(>F)`[2]
Pval_21_inter[kcond]=gyx$`Pr(>F)`[2]


gxy = grangertest(x,y,order = 1) #x causes y 
gyx = grangertest(y,x,order = 1) #y causes x

Pval_12_inter_lag1[kcond]=gxy$`Pr(>F)`[2]
Pval_21_inter_lag1[kcond]=gyx$`Pr(>F)`[2]

}

########################################################################################################################
########## Sugihara two species deterministic competition model -- no interactions
########################################################################################################################


for (kcond in 1:ncond){
  tryCatch({
  tmax=800
  X=Y=rep(1,tmax)## Problem if I start with 1
  X[1]=runif(1,0.1,0.7)
  Y[1]=runif(1,0.1,0.7)
  for (t in 1:tmax)
  {
    X[t+1]=X[t]*(3.8-3.8*X[t]-0*Y[t])
    Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0*X[t])
  }
#  X
#  Y
  x=log(X[501:800])
  y=log(Y[501:800])
  # centering
  x=x-mean(x)
  y=y-mean(y)
  z=cbind(x,y)
  
  varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
  lag_order_noInter[kcond] <- varcompet$p

  gxy = grangertest(x,y,order = lag_order_noInter[kcond]) #x causes y 
  gyx = grangertest(y,x,order = lag_order_noInter[kcond]) #y causes x
  
  Pval_12_noInter[kcond]=gxy$`Pr(>F)`[2]
  Pval_21_noInter[kcond]=gyx$`Pr(>F)`[2]

gxy = grangertest(x,y,order = 1) #x causes y 
gyx = grangertest(y,x,order = 1) #y causes x

Pval_12_noInter_lag1[kcond]=gxy$`Pr(>F)`[2]
Pval_21_noInter_lag1[kcond]=gyx$`Pr(>F)`[2]


  },
  error=function(e) {
    message(paste("Difficulty in fitting MAR(p) model", varcompet)) #Singularity in the matrix
    # Choose a return value in case of error
    return(z)
  })
}

library(sm)

par(mfrow=c(1,2))
sm.density.compare(c(Pval_12_inter,na.omit(Pval_12_noInter)), group = c(1,0), model = "none",lwd=2)
sm.density.compare(c(Pval_21_inter,na.omit(Pval_21_noInter)), group = c(1,0), model = "none",lwd=2)
### Poor representation

DataCompet_sugiharaDeterModel = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)


### Let's try the percentage of initial conditions at which the Granger tests works
sum(Pval_12_inter<0.1)/length(Pval_12_inter) #100% at 0.1 level
sum(Pval_21_inter<0.1)/length(Pval_21_inter) #52% at 0.1 level

sum(na.omit(Pval_12_noInter>0.1))/length(na.omit(Pval_12_noInter)) #90% at 0.1 level
sum(na.omit(Pval_21_noInter>0.1))/length(na.omit(Pval_21_noInter)) #87% at 0.1 level

plot(Pval_12_inter,Pval_21_inter,xlim=c(0,1),ylim=c(0,1))
plot(Pval_12_noInter,Pval_21_noInter,xlim=c(0,1),ylim=c(0,1))
### Second hypothesis consistent with the null... 

DataCompet_sugiharaDeterModel_Granger = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)

#Write down results
write.csv(DataCompet_sugiharaDeterModel_Granger,file="results/explo/DataCompet_deterModel_Granger.csv")
