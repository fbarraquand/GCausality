### CP April 2019
### Plot time series for 2 species (Veilleux, chaotic deterministic model, 2-species-and-a-driver model)
### All simulations are from theoreticalSugiharaModels.R (chaotic model) and StochasticCompet.R, FBarraquand 2018

graphics.off()
rm(list=ls())

pdf("fig/time_series_small_dim_log.pdf",width=14,height=8)
par(mfrow=c(2,2),cex=1.25,mar=c(2.5,4,1.5,0.5))

#### Veilleux
#0.5
DB=read.table("data/veilleux_subset_CC05a.txt",sep="")
time=DB[,1]

DB[,2]=log(DB[,2])
DB[,2]=DB[,2]-mean(DB[,2])
DB[,3]=log(DB[,3])
DB[,3]=DB[,3]-mean(DB[,3])

plot(time,DB[,2],ylim=c(min(DB[,3]),max(DB[,2])),ylab="ln(abundance)",xlab="")
lines(time,DB[,2],col="blue",type="o",lwd=3,pch=16)
lines(time,DB[,3],col="red",type="o",lwd=3,pch=16)
legend("bottomleft",c("Paramecium","Didinium"),pch=16,col=c("blue","red"),pt.bg=c("blue","red"),bty="n",inset=c(0.075,0.01))
mtext("a)",side=2,line=2,at=max(DB[,2]),las=2,cex=1.2)

#0.375
DB=read.table("data/veilleux_subset.txt",sep="")
time=DB[,1]

DB[,2]=log(DB[,2])
DB[,2]=DB[,2]-mean(DB[,2])
DB[,3]=log(DB[,3])
DB[,3]=DB[,3]-mean(DB[,3])

plot(time,DB[,2],ylim=c(min(DB[,3]),max(DB[,3])),ylab="",xlab="")
lines(time,DB[,2],col="blue",type="o",lwd=3,pch=16)
lines(time,DB[,3],col="red",type="o",lwd=3,pch=16)
legend("bottomleft",c("Paramecium","Didinium"),pch=16,col=c("blue","red"),pt.bg=c("blue","red"),bty="n",inset=c(0.075,0.01))
mtext("b)",side=2,line=2,at=max(DB[,3]),cex=1.2,las=2)


par(mar=c(4,4,0,0.5))
####Chaos example
tmax=300
X=Y=rep(1,tmax)## Problem if I start with 1
X[1]=0.1
Y[1]=0.2
#X[1]=0.1
#Y[1]=0.1
for (t in 1:tmax)
{
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.02*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.1*X[t])
}
x=log(X)
y=log(Y)
# centering
x=x-mean(x)
y=y-mean(y)

z=cbind(x,y)
#z=cbind(X,Y)
vec_col=c("black","red")
minz=min(z[151:201,])
maxz=max(z[151:201,])

#plot(151:201,X[151:201],type="o",col="blue",lwd=3,xlab="Time",ylab="abundance",ylim=c(minz-0.1,maxz),pch=16)
#lines(151:201,Y[151:201],type="o",col="red",lwd=3,pch=16)
  plot(151:201,z[151:201,1],type="o",col="blue",lwd=3,ylim=c(minz,maxz),xlab="Time",ylab="ln(abundance)",pch=16)
  lines(151:201,z[151:201,2],type="o",col="red",lwd=3,pch=16)
legend("bottomleft",c("X","Y"),pch=16,col=c("blue","red"),bty="n",inset=c(0.075,0.01))
mtext("c)",side=2,line=2,at=max(maxz),cex=1.2,las=2)
mtext("Time",side=2,line=2,at=(151+201)/2,cex=1.2,las=2)


####... and a driver
set.seed(42)

tmax=300
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))
seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates
### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
y1<-seasonality+y1noise

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-2*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0.31*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
}
y=log(Y)
y[,1]=y[,1]-mean(y[,1])
y[,2]=y[,2]-mean(y[,2])

plot(100:200,y1[100:200],ylab="",col="black",type="o",xlab="Time",ylim=c(min(y[100:200,]),max(y[100:200,])),lwd=3,pch=16)
lines(100:200,y[100:200,1],col="blue",lwd=3,t="o",pch=16)
lines(100:200,y[100:200,2],col="red",lwd=3,t="o",pch=16)
mtext("d)",side=2,line=2,at=max(y[100:200]),cex=1.2,las=2)
mtext("Time",side=2,line=2,at=(100+200)/2,cex=1.2,las=2)
legend("bottomleft",c("driver","Y1","Y2"),pch=16,col=c("black","blue","red"),bty="n",inset=c(0.075,0.01))
dev.off()

