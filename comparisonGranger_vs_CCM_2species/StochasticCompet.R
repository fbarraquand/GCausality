########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################
### FB 14/08/2018  - We now include a stochastic model in this code

library("vars")

########################################################################################################################
####################### Repeating many times with many initial conditions / stochastic sequences    ####################
########################################################################################################################

#Initializing vectors 
ncond<-500
Pval_12_inter=Pval_21_inter=Pval_12_noInter=Pval_21_noInter=lag_order_inter=lag_order_noInter=rep(NA,ncond)

########################################################################################################################
########## Stochastic two species competition model -- interactions present
########################################################################################################################


for (kcond in 1:ncond){

tmax=800
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

}

########################################################################################################################
########## Stochastic two species competition model -- no interactions
########################################################################################################################

for (kcond in 1:ncond){
  
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
  
}

par(mfrow=c(1,2))
sm.density.compare(c(Pval_12_inter,Pval_12_noInter), group = c(1,0), model = "none",lwd=2)
sm.density.compare(c(Pval_21_inter,Pval_21_noInter), group = c(1,0), model = "none",lwd=2)
### Poor representation

DataCompet_sugiharaDeterModel = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)


### Let's try the percentage of initial conditions at which the Granger tests works
sum(Pval_12_inter<0.1)/length(Pval_12_inter) #100% at 0.1 level
sum(Pval_21_inter<0.1)/length(Pval_21_inter) #52% at 0.1 level

sum(Pval_12_noInter>0.1)/length(Pval_12_noInter) #90% at 0.1 level
sum(Pval_21_noInter>0.1)/length(Pval_21_noInter) #87% at 0.1 level

### Plot 4 panel figure with (a) P-values 1->2, (b) P-values 2->1 interactions and no interactions
### Other plot with lag order - should I overlay lines??
### Or should I do biplots? 
plot(Pval_12_inter,Pval_21_inter,xlim=c(0,1),ylim=c(0,1))
plot(Pval_12_noInter,Pval_21_noInter,xlim=c(0,1),ylim=c(0,1))
### Second hypothesis consistent with the null... 


DataCompet_stochModel_Granger = data.frame(lag_order_inter,Pval_12_inter,Pval_21_inter,lag_order_noInter,Pval_12_noInter,Pval_21_noInter)
#Write down results
write.csv(DataCompet_stochModel_Granger,file="results/DataCompet_stochModel_Granger.csv")

