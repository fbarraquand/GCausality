###########################################################################################################
########### Granger-style analysis of Veilleux's data -- FBarraquand. Started in 2015 #####################
########### Re-coded 05/10/2017 ###########################################################################
###########################################################################################################

rm(list=ls())
#options(scipen=999)

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
### Lag order = 2
lag_order = 2
#### Using grangertest() 

library(lmtest)
grangertest(x,y,order = 1)
grangertest(y,x,order = 1)

gxy = grangertest(x,y,order = 2) #x causes y 
gyx = grangertest(y,x,order = 2) #y causes x

Pval_preyToPred=gxy$`Pr(>F)`[2]
Pval_predToPrey=gyx$`Pr(>F)`[2]

DataSetV = "CC_0.375"
GrangerVeilleux = data.frame(DataSetV,lag_order,Pval_preyToPred,Pval_predToPrey,stringsAsFactors=FALSE)

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

library(rEDM)

### Define what we use as library (fitted indices) and predicted ones
lib <- c(1, 65)
pred <- c(201, 500) ### the code will overload this and enable leave-one-out cross-validation, according to
### https://cran.r-project.org/web/packages/rEDM/vignettes/rEDM_tutorial.html

pred_prey=data.frame(1:n,z[1,],z[2,])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.5))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### But this is on the log scale...
# Let's go back to the usual scale

pred_prey=data.frame(1:n,exp(z[1,]),exp(z[2,]))
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.5))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### OK let's try with the real data to compare

pred_prey=data.frame(1:n,x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 60, by = 10), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.5))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

# Barely better -- perhaps this is because I logged the data

pred_prey=data.frame(DB[,1],DB[,2],DB[,3])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 60, by = 1), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 60, by = 1), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.5))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Not what we find in Sugihara et al. - but Veilleux had several time series. 
time

### Let's check it is the same data
pdf(file = "Veilleux_Paramecium_Didinium_untransformed.pdf",height=8,width=10)
par(cex=1.5)
plot(time,DB[,2],ylim=c(min(DB[,3]),max(DB[,2])),ylab="abundance")
lines(time,DB[,2],col="blue",type="o",lwd=5,pch=21)
lines(time,DB[,3],col="red",type="o",lwd=5,pch=21)
legend("bottomright",c("Paramecium","Didinium"),pch=21,col=c("blue","red"),pt.bg=c("blue","red"))
dev.off()

### Not exactly the same time series -- let's use the exact same one

###########################################################################################################
#### Analysis of the second Veilleux dataset 
###########################################################################################################

DB=read.table("veilleux_subset_CC05a.txt",sep="")
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
#Plot data 
pdf(file = "Veilleux_Paramecium_Didinium_CC05a.pdf",height=8,width=10)
par(cex=1.5)
plot(time,x,ylim=c(min(y),max(x)),ylab="ln(abundance)")
lines(time,x,col="blue",type="o",lwd=5,pch=21)
lines(time,y,col="red",type="o",lwd=5,pch=21)
legend("bottomright",c("Paramecium","Didinium"),pch=21,col=c("blue","red"),pt.bg=c("blue","red"))
dev.off()
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)
length(x)

pdf(file = "Veilleux_Paramecium_Didinium_CC05a_untransformed.pdf",height=8,width=10)
par(cex=1.5)
plot(time,DB[,2],ylim=c(min(DB[,3]),max(DB[,2])),ylab="abundance")
lines(time,DB[,2],col="blue",type="o",lwd=5,pch=21)
lines(time,DB[,3],col="red",type="o",lwd=5,pch=21)
legend("bottomright",c("Paramecium","Didinium"),pch=21,col=c("blue","red"),pt.bg=c("blue","red"))
dev.off()

### Note: They conveniently removed the first 10 points - let's do the same 
### This will good from a statistical viewpoint at least because the TS would then have same length
n=nrow(DB)
n
DB=DB[10:n,]
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
n=nrow(DB)
n
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)
length(x)

### VAR analysis of the CC05a data
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")
varpp<-VAR(y=predprey, type="none",ic="SC",lag.max=5)
summary(varpp)
coef(varpp) ##
causality(varpp,cause="x")
causality(varpp,cause="y")
### Seems to me we can reject the null that y does not cause x and vice versa, at least at the 0.05 level.
lag_order = 1
gxy = grangertest(x,y,order = lag_order) #x causes y 
gyx = grangertest(y,x,order = lag_order) #y causes x

Pval_preyToPred=gxy$`Pr(>F)`[2]
Pval_predToPrey=gyx$`Pr(>F)`[2]

DataSetV = "CC_0.5a"
GrangerVeilleux = rbind(GrangerVeilleux ,c(DataSetV,lag_order,Pval_preyToPred,Pval_predToPrey))
GrangerVeilleux$DataSetV[2] = "CC_0.5a"
write.csv(GrangerVeilleux,"results/GrangerVeilleux.csv")

IC=VARselect(y=predprey, type="none",lag.max=15) ## selection by AIC (not even AICc)
IC
matplot(1:15,scale(t(IC$criteria)))

crit=scale(t(IC$criteria))

pdf(file="LagOrderSelection2Species_CC05a.pdf",width=8,height=6)
par(cex=1.5,lwd=3)
plot(1:15,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-2,2))
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i+1])}
legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18)
dev.off()

residuals(varpp)
acf(residuals(varpp)) #OK

### Simulation according to VAR(1)

#Extracting model coefficients
a_11_l1=coef(varpp)$x[1,1]
a_12_l1=coef(varpp)$x[2,1]

a_21_l1=coef(varpp)$y[1,1]
a_22_l1=coef(varpp)$y[2,1]

Sigma=summary(varpp)$covres

############### Simulations of our fitted models ##################################################
#### Full MAR(2) model
z=matrix(0,nrow=2,ncol=nrow(DB))
z[,1]=runif(2,0,1)
z
n=nrow(DB)
for (t in 1:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + a_12_l1*z[2,t] +eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] +eps[2]
}
matplot(t(z))
matlines(t(z)) ## Fairly predator-prey cycle looking

### CCM this simulation now 

pred_prey=data.frame(1:n,z[1,],z[2,])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### But this is on the log scale...
# Let's go back to the usual scale

pred_prey=data.frame(1:n,exp(z[1,]),exp(z[2,]))
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### OK let's try with the real data to compare

pred_prey=data.frame(1:n,x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 62, by = 1), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### With the original data, non-logged

pred_prey=data.frame(DB[,1],DB[,2],DB[,3])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(10, 62, by = 1), num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(10, 62, by = 1), num_samples=100)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred, 
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey, 
                                                                               tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2) 
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)

########## FB stopped there 18:40 06/10/2017 ############################################

### Produce figure -- 

#pdf(file="CCM_VeilleuxData_and_VARsim.pdf",height=10,width=10)


pdf(file="CCM_VeilleuxData.pdf",height=5,width=10)

par(mfrow=c(1,2),cex=1.25)
### Veilleux data -- logged first

pred_prey=data.frame(1:n,x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(5, 62, by = 1), num_samples=100)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred, 
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey, 
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "Real data (logged)")
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2) 
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)


### Try to output the confidence intervals
### With the original data, non-logged

pred_prey=data.frame(DB[,1],DB[,2],DB[,3])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(5, 62, by = 1), num_samples=100)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred, 
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey, 
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "Real data")
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2) 
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
dev.off()

### VAR model -- log-scale first
pdf(file="CCM_VARsim.pdf",height=5,width=10)

par(mfrow=c(1,2),cex=1.25)
pred_prey=data.frame(1:n,z[1,],z[2,])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(5, 62, by = 1), num_samples=100)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred, 
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey, 
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main ="VAR(1) - log(abundance)")
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2) 
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)


### Try to output the confidence intervals
### But this is on the log scale...
# Let's go back to the usual scale

pred_prey=data.frame(1:n,exp(z[1,]),exp(z[2,]))
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey", 
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey", 
                      lib_sizes = seq(5, 62, by = 1), num_samples=100)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred, 
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey, 
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "VAR(1) - abundance")
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2) 
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue", 
      lty = 2, lwd = 2)

dev.off()

############ Output the nonlinearity parameters -- do that later // code to modify below. 
# smap_output <- s_map(composite_ts, composite_lib, composite_pred, E = 8)
# par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
# plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", 
#      ylab = "Forecast Skill (rho)")
