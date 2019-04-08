graphics.off()
rm(list=ls())

set.seed(42)
library(vars)
pdf("../figures/lag_order_large_communities.pdf",width=5,height=10)
par(mfrow=c(2,1),cex=1.25,lwd=2,mar=c(2.5,4,2,0.5))
#10 species, ref
tmax=300 ## We are very optimistic, Sugihara et al. also used 3000. 
nspecies=10
sigma=0.3 # Sigma^2= 0.1
Y=matrix(1,nrow=tmax,ncol=nspecies)
for (i in 1:nspecies){
  Y[1,i]=abs(rnorm(1,1,1))}

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(4-4*Y[t,1]-2*Y[t,2]-0.4*Y[t,3]+rnorm(1,0,sigma))
  Y[t+1,2] = Y[t,2]*exp(3.1-0.31*Y[t,1]-3.1*Y[t,2]-0.93*Y[t,3]+rnorm(1,0,sigma))
  Y[t+1,3] = Y[t,3]*exp(0.12 + 0.636*Y[t,1]+0.636*Y[t,2]-2.12*Y[t,3]+rnorm(1,0,sigma))
  ### Connecting the two sets by being influenced by all three and having an effect on 5
  Y[t+1,4] = Y[t,4]*exp(3-0.111*Y[t,1]-0.111*Y[t,2]+0.131*Y[t,3]-3.8*Y[t,4]+rnorm(1,0,sigma))
  ### Influenced by 4 and also 5,6,7
  Y[t+1,5] = Y[t,5]*exp(3-2*Y[t,5]-2*Y[t,6]-0.4*Y[t,7]+0.5*Y[t+1,4]+ rnorm(1,0,sigma))
  ### 
  Y[t+1,6] = Y[t,6]*exp(3.1-0.31*Y[t,5]-3.1*Y[t,6]-0.93*Y[t,7]+rnorm(1,0,sigma))
  Y[t+1,7] = Y[t,7]*exp(0.8 + 0.636*Y[t,5]+0.636*Y[t,6]-2.12*Y[t,7]+rnorm(1,0,sigma))
  Y[t+1,8] = Y[t,8]*exp(4-4*Y[t,8]-2*Y[t,9]-0.4*Y[t,10]+rnorm(1,0,sigma))
  Y[t+1,9] = Y[t,9]*exp(3.1-0.31*Y[t,8]-3.1*Y[t,9]-0.93*Y[t,10]+rnorm(1,0,sigma))
  Y[t+1,10] = Y[t,10]*exp(0.12 + 0.636*Y[t,8]+0.636*Y[t,9]-2.12*Y[t,10]+rnorm(1,0,sigma))

}

y=log(Y)
y=y-mean(y)

yd=data.frame(y)

IC=VARselect(y=yd, type="none",lag.max=15) ## selection by AIC (not even AICc)
crit=scale(t(IC$criteria))

plot(1:15,crit[,1],ylab="Information Criteria",xlab="",type="o",ylim=c(-2,2.25),col="black")
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i])}
legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18,bty="n")
mtext("a)",las=2,side=2,line=2.5,at=2.25)


par(mar=c(4,4,0.5,0.5))

#20 species, 1 community
ksite=1
pathrandomLV = "../../20species/data/Data_wTime_abs_LV.csv"
pathrandomVAR = "../../20species/data/Data_wTime_abs_VAR.csv"
#for (path in c(pathrandomLV,pathrandomVAR))
for (path in c(pathrandomLV))
    {

    ### Selects the files and then time series

    DBall=read.csv(path)
    DB=DBall[DBall$Site==ksite,] ## Select a site
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps
    abundance_mat=as.matrix(DB[,4:23]) ### Create matrix with time series of abundances
IC=VARselect(y=abundance_mat, type="none",lag.max=15) ## selection by AIC (not even AICc)
crit=scale(t(IC$criteria))

plot(1:15,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-2,2.25),col="black")
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i])}
#legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18)

mtext("b)",las=2,side=2,line=2.5,at=2.25)

}

dev.off()
