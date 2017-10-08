### FB 08/10/2017 -- adapted from Nantes et Irstea presentations though
### Uses a strongly seasonal driver. A logical case for conditional GC. 
library(vars)

### First case with interactions

set.seed(42)

tmax=300
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))
seasonality<-2*sin(2*pi*(1:tmax)/24)	# must be enough to affect the growth rates
### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-2*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0.31*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
}
y=log(Y)

varcompet2<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
causality(varcompet2,cause="y1") #p-value 0.03859
causality(varcompet2,cause="y2") #3.443e-13
### I may need to repeat this over many simulations

############
pdf(file="2competitors_wSharedSeasonalDriver.pdf",width=12,height=6)
par(cex=1.5,lwd=2)
plot(1:tmax,y1,ylab="Shared driver u_t",col="green",type="l",xlab="Time",ylim=c(-10,max(y1)))
lines(1:tmax,y[,1],col="blue")
lines(1:tmax,y[,2],col="black")
dev.off()

### Second case without interactions
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))


for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-0*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
}
y=log(Y)

pdf(file="2competitors_wSharedSeasonalDriver_noInteractions.pdf",width=12,height=6)
par(cex=1.5,lwd=2)
plot(1:tmax,y1,ylab="Shared driver u_t",col="green",type="l",xlab="Time",ylim=c(-10,max(y1)))
lines(1:tmax,y[,1],col="blue")
lines(1:tmax,y[,2],col="black")
dev.off()

varcompet2bis<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
summary(varcompet2bis)
causality(varcompet2bis,cause="y1") #p-value 1.073e-06
causality(varcompet2bis,cause="y2") #0.1953


### Let's use another simulation with just a little more noise

Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-0*Y[t,2] + rnorm(1,0,0.3))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.3))
}
y=log(Y)

varcompet2ter<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
summary(varcompet2ter)
causality(varcompet2ter,cause="y1") #p-value 6.263e-10
causality(varcompet2ter,cause="y2") #0.5799

### So it is possible to craft examples that are difficult to infer with GC 
### --- still we may need many simulations and parameter values for a good comparison of CCM and GC

### NB might well be possible that we would get better results by restricting the fit to a VAR1
### But we have no theoretical reason to do that

############## Two things to check though!!! ###################
### I haven't tested pairwise GC on these examples. Perhaps it works better than conditional GC? 
### And I haven't CCM this, on the other hand. 
### They say you don't need necessarily corrections for seasonality in the 2012 paper -- or do we? cf Deyle PNAS Influenza

