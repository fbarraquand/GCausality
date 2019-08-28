###########################################################################################################
########### CCM analysis of Veilleux's data set                                                ############ 
########### From FBarraquand "predatorPrey_Veilleux.R"                                         ############
########### CP April 2019: Plot the CCM for the two Veilleux data sets (0.5 and 0.375),        ############
########### with confidence intervals, comparing with simulated values vrom VAR(p) models      ############
###########################################################################################################

graphics.off()
rm(list=ls())

library(vars)
library(rEDM)

pdf("fig/Veilleux_simulation_CCM.pdf",height=10,width=10)
par(mfrow=c(2,2),cex=1.25,mar=c(1,4,3,0.5))


#O5
DB=read.table("data/veilleux_subset_CC05a.txt",sep="")
n=nrow(DB)
DB=DB[10:n,] #Sugihara et al. removed the first 10 lines
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)


#################CCM
pred_prey=data.frame(1:length(x),x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=FALSE)
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "Real data",xlab="",xaxt="n",xlim=c(5,60))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "bottomright", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, bty="n")
mtext("a)",side=2,line=2,at=1.1,las=2,cex=1.2)


### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)

#################Simulated data
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")

varpp<-VAR(y=predprey, type="none",ic="SC",lag.max=5)
#Extracting model coefficients
a_11_l1=coef(varpp)$x[1,1]
a_12_l1=coef(varpp)$x[2,1]

a_21_l1=coef(varpp)$y[1,1]
a_22_l1=coef(varpp)$y[2,1]

Sigma=summary(varpp)$covres

############### Simulations of our fitted models ##################################################
#### Full MAR(2) model
plot(0, 0, type = "n", xlab = "", ylab = "", ylim = c(0, 1.1),xlim=c(5,60),xaxt="n",main="MAR-simulated data")
for(simu in 1:100){
print(simu)
z=matrix(0,nrow=2,ncol=nrow(DB))
z[,1]=runif(2,0,1)
#z[,2]=runif(2,0,1)
n=nrow(DB)
for (t in 1:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + a_12_l1*z[2,t] +eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] +eps[2]
}

pred_prey=data.frame(1:n,z[1,],z[2,])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), replace=FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 62, by = 1), replace=FALSE)

prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))


if(simu==1){
alpha=1
}else{
alpha=0.1
}
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
lines(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = rgb(1,0,0,alpha))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), type="l",col = rgb(0,0,1,alpha))
}
mtext("b)",side=2,line=2,at=1.1,las=2,cex=1.2)

###################0375
DB=read.table("data/veilleux_subset.txt",sep="")
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)


#################CCM
par(mar=c(4,4,0,0.5))
pred_prey=data.frame(1:length(x),x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=FALSE)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "",xlab="Library size",xlim=c(5,60))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")

### Try to output the confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
mtext("c)",side=2,line=2,at=1.1,las=2,cex=1.2)



#################Simulated data
predprey=data.frame(cbind(x,y))
names(predprey)=c("x","y")

varpp<-VAR(y=predprey, type="none",ic="SC",lag.max=5)
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
plot(0, 0, type = "n", xlab = "Library size", ylab = "", ylim = c(0, 1.1),xlim=c(5,60))
for(simu in 1:100){
print(simu)
z=matrix(0,nrow=2,ncol=nrow(DB))
z[,1]=runif(2,0,1)
z[,2]=runif(2,0,1)
n=nrow(DB)
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + a_12_l1*z[2,t] + a_11_l2*z[1,t-1] + a_12_l2*z[2,t-1]+eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] + a_21_l2*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}

pred_prey=data.frame(1:n,z[1,],z[2,])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 62, by = 1), replace=FALSE)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 62, by = 1), replace=FALSE)

prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))


if(simu==1){
alpha=1
}else{
alpha=0.1
}
lines(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = rgb(1,0,0,alpha))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), type="l",col = rgb(0,0,1,alpha))
}
mtext("d)",side=2,line=2,at=1.1,las=2,cex=1.2)
dev.off()
