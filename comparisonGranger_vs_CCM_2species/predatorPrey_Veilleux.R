###########################################################################################################
########### Granger-style analysis of Veilleux's data -- FBarraquand. Started in 2015 #####################
########### Re-coded 05/10/2017 ###########################################################################
###########################################################################################################

rm(list=ls())
options(scipen=999)

### Loading the Veilleux data
DB=read.table("veilleux_subset.txt",sep="")
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
#Plot data 
pdf(file = "Veilleux_Paramecium_Didinium.pdf",height=8,width=10)
par(cex=1.5)
plot(time,x,ylim=c(min(y),max(x)),ylab="ln(abundance)")
lines(time,x,col="blue",type="o",lwd=5,pch=21)
lines(time,y,col="red",type="o",lwd=5,pch=21)
legend("bottomright",c("Paramecium","Didinium"),pch=21,col=c("blue","red"),pt.bg=c("blue","red"))
dev.off()
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)
length(x)

#### Using grangertest() 

library(lmtest)
grangertest(x,y,order = 1)
grangertest(y,x,order = 1)

grangertest(x,y,order = 2)
grangertest(y,x,order = 2)

### Using VARS to estimate simultaneously model order
library("vars")
data(Canada)
VAR(Canada, p = 2, type = "none") #example
# format data
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")
varpp<-VAR(y=predprey, type="none",ic="SC",lag.max=5)
summary(varpp)
coef(varpp) ##
causality(varpp,cause="x")
causality(varpp,cause="y")
### Seems to me we can reject the null that y does not cause x and vice versa, at least at the 0.05 level.
varpp<-VAR(y=predprey, type="none",ic="SC",lag.max=5)
varpp

IC=VARselect(y=predprey, type="none",lag.max=15) ## selection by AIC (not even AICc)
IC
matplot(1:15,scale(t(IC$criteria)))

crit=scale(t(IC$criteria))

pdf(file="LagOrderSelection2Species.pdf",width=8,height=6)
par(cex=1.5,lwd=3)
plot(1:15,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-2,2))
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i+1])}
legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18)
dev.off()

residuals(varpp)
acf(residuals(varpp))

varpp1<-VAR(y=predprey, type="none",p=1) #might be a better model, but let's see.
residuals(varpp1)
acf(residuals(varpp1)) ## OK, some remaining autocorr, must be the reason. 

#Let's roll with order 2 for now

#Extracting model coefficients
a_11_l1=coef(varpp)$x[1,1]
a_12_l1=coef(varpp)$x[2,1]
a_11_l2=coef(varpp)$x[3,1]
a_12_l2=coef(varpp)$x[4,1]

a_21_l1=coef(varpp)$y[1,1]
a_22_l1=coef(varpp)$y[2,1]
a_21_l2=coef(varpp)$y[3,1]
a_22_l2=coef(varpp)$y[4,1]

Sigma=summary(varpp)$covres

############### Simulations of our fitted models ##################################################
#### Full MAR(2) model
z=matrix(0,nrow=2,ncol=nrow(DB))
z[,1]=runif(2,0,1)
z[,2]=runif(2,0,1)
z
n=nrow(DB)
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + a_12_l1*z[2,t] + a_11_l2*z[1,t-1] + a_12_l2*z[2,t-1]+eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] + a_21_l2*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z))

### Let's roll with that - obviously the data are less variable, but we shall see what happens when we CCM z

