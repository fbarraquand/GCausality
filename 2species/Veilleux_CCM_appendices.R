###########################################################################################################
########### CCM analysis of Veilleux's data set                                                ############ 
########### From FBarraquand "predatorPrey_Veilleux.R"                                         ############
########### CP April 2019: Plot the CCM for the two Veilleux data sets (0.5 and 0.375),        ############
########### with confidence intervals, comparing raw values and logged(values)                 ############
###########################################################################################################

graphics.off()
rm(list=ls())

b_replace=FALSE
#If replace is false, one data point cannot be used several times when building a library in CCM analyses

library(rEDM)

if(b_replace){
pdf("fig/CCM_Veilleux2DataSets_withreplace.pdf",height=10,width=10)
}else{
pdf("fig/CCM_Veilleux2DataSets_noreplace.pdf",height=10,width=10)
}
par(mfrow=c(2,2),cex=1.25,mar=c(1,4,3,0.5))



#O5
DB=read.table("data/veilleux_subset_CC05a.txt",sep="")
n=nrow(DB)
DB=DB[10:n,] #Sugihara et al. removed the first 10 lines
time=DB[,1]
x=log(DB[,2]) #prey 
y=log(DB[,3]) #predator
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)

pred_prey=data.frame(1:length(x),x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=b_replace)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=b_replace)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1), main = "Real data (logged)",xlab="",xaxt="n",xlim=c(5,60))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "bottomright", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, bty="n")
legend("topleft","a)",bty="n")

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

par(mar=c(1,1,3,3))

pred_prey=data.frame(DB[,1],DB[,2],DB[,3])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=b_replace)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), num_samples=100,replace=b_replace)
### num_samples=100 necessary to estimate sd.rho
prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", ylim = c(0, 1.1), main = "Real data",xaxt="n",yaxt="n",xlim=c(5,60),xlab="",ylab="")
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend("topleft","b)",bty="n")

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)


#0375
par(mar=c(4,4,0,0.5))
DB=read.table("data/veilleux_subset.txt",sep="")
time=DB[,1]
x=log(DB[,2]) #prey # also done with non-log transformed data. 
y=log(DB[,3]) #predator
x=x-mean(x) ### centering - very important for many tests. 
y=y-mean(y)

pred_prey=data.frame(1:length(x),x,y)
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1) ,replace=b_replace,num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), replace=b_replace,num_samples=100)

prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))



plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1),xlim=c(5,60))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend("topleft","c)",bty="n")


#confidence intervals
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)


par(mar=c(4,1,0,3))
pred_prey=data.frame(DB[,1],DB[,2],DB[,3])
names(pred_prey)=c("time","prey","pred")

prey_xmap_pred <- ccm(pred_prey, E = 3, lib_column = "prey",
                      target_column = "pred", lib_sizes = seq(5, 60, by = 1), replace=b_replace,num_samples=100)

pred_xmap_prey <- ccm(pred_prey, E = 3, lib_column = "pred", target_column = "prey",
                      lib_sizes = seq(5, 60, by = 1), replace=b_replace,num_samples=100)

prey_xmap_pred_means <- data.frame(ccm_means(prey_xmap_pred), sd.rho = with(prey_xmap_pred,
                                                                            tapply(rho, lib_size, sd)))
pred_xmap_prey_means <- data.frame(ccm_means(pred_xmap_prey), sd.rho = with(pred_xmap_prey,
                                                                            tapply(rho, lib_size, sd)))


plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "", ylim = c(0, 1.1),yaxt="n",xlim=c(5,60))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend("topleft","d)",bty="n")

#confidence intervals

lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho + 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(prey_xmap_pred_means$lib_size, prey_xmap_pred_means$rho - 2*prey_xmap_pred_means$sd.rho, col = "red",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho + 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)
lines(pred_xmap_prey_means$lib_size, pred_xmap_prey_means$rho - 2*pred_xmap_prey_means$sd.rho, col = "blue",
      lty = 2, lwd = 2)





dev.off()
