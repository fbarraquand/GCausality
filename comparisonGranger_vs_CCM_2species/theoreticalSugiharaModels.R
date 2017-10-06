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
pdf(file="CompetitionDeterministicNL2Species.pdf",width=12,height=6)
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
########## Sugihara two species deterministic competition model -- no interactions
########################################################################################################################

### -- to fill --- ###

########################################################################################################################
########## Stochastic competition model -- interactions
########################################################################################################################

### -- to fill --- ###

########################################################################################################################
########## Stochastic competition model -- no interactions
########################################################################################################################

### -- to fill --- ###