###CP April 2019, based on FB's previous work
### Plot the values of different information criteria depending on the lag chosen in VAR(p) model for Veilleux data sets

graphics.off()
rm(list=ls())

library(vars)

pdf("Lag_choice_veilleux.pdf")
par(mfrow=c(2,1),cex=1.25,mar=c(2,4,2,0.5),lwd=3)
######0.5
DB=read.table("veilleux_subset_CC05a.txt",sep="")
time=DB[,1]

x=log(DB[,2])
x=x-mean(x)
y=log(DB[,3])
y=y-mean(y)
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")

IC=VARselect(y=predprey, type="none",lag.max=15) ## selection by AIC (not even AICc)
crit=scale(t(IC$criteria))

plot(1:15,crit[,1],ylab="Information Criteria",xlab="",type="o",ylim=c(-2,2),col="black",xaxt="n")
axis(1,at=seq(2,14,2),lab=rep("",7))
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i])}
legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18,bty="n")
mtext("a)",side=2,line=2.5,at=2,las=2,cex=1.2)



par(mar=c(4,4,0,0.5))
######0.375
### Loading the Veilleux data
DB=read.table("veilleux_subset.txt",sep="")
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
#Plot data 
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")

IC=VARselect(y=predprey, type="none",lag.max=15) ## selection by AIC (not even AICc)
crit=scale(t(IC$criteria))

plot(1:15,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-2,2),col="black")
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i])}
mtext("b)",side=2,line=2.5,at=2,las=2,cex=1.2)


dev.off()

