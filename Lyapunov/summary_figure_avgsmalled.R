####################### CP 13/06/2019
####################### This script uses the dynamics and lyapunov exponents designed and computed by Frederic Barraquand in order to plot the relation between lyapunov (degree of chaos) and recall+sensitivity of Granger Causality and Convergent Cross Mapping

rm(list=ls())
graphics.off()

#Small dimension
l_chaos=0.41 #0.43 without interactions
l_2species=-0.18 #0.35 without interactions
l_driver=-0.14 #0.19 without interactions
alpha_s=0.1
threshold_s=0.04

#Large dimension
l_10=0.33
l_20=-0.075

#Function to compute recall and specificity (the 2 following functions were written by FB)
ratesClassif <- function (estimated_mat,true_binary_mat)
{
  if (typeof(estimated_mat)=="double")
  { Pos = (abs(estimated_mat)>0)} ## estimated links
    else if (typeof(estimated_mat)=="logical")
    {Pos = estimated_mat}
    else {print("Unknow type for estimated interaction matrix")}

  TruePos =  (abs(true_binary_mat)>0) ## work as well with non-binary
  TP = sum((TruePos == TRUE) & (Pos == TRUE))
  FN = sum((TruePos == TRUE) & (Pos == FALSE))
  FP = sum((TruePos == FALSE) & (Pos == TRUE))
  TN = sum((TruePos == FALSE) & (Pos == FALSE))
  return(c(TP,FN,FP,TN))
}

#### Devise function to compute recall, precision, and ROC curves
diagnosticsClassif<- function (vector_classif){
  TP = vector_classif[1]
  FN = vector_classif[2]
  FP = vector_classif[3]
  TN = vector_classif[4]

  FPR = FP / (FP + TN)  # = fall-out rate  = 1 - specificity = 1 - true negative rate
  TPR = TP / (TP+FN) # = Recall = True positives over relevant elements (power to find the true positives)
  ### ROC = c(FPR,TPR)
  Precision = TP/(TP+FP)  #True positives over tested positives (percentage positives correct)
  return(c(FPR,TPR,Precision))
}


pdf("diagnostics_vs_lyapunov.pdf")
#pdf("total_recall_smalldim.pdf")
par(mfrow=c(2,2))

#Small dim
plot(0,0,xlab="Lyapunov exponent",ylab="Recall",t="n",xlim=c(0.5,3.5),ylim=c(0,1))
#Chaos
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
tp=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s),((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s),((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(1,rec,pch=16,col="blue")

tp=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1),((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1),((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(1,rec,pch=18,col="red")

#Stochastic 2 species
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
tp=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s),((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s),((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(2,rec,pch=16,col="blue")

tp=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1),((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1),((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(2,rec,pch=18,col="red")

#Stochastic and driver
tab_GC=read.csv('../twoSpecies_andDriver/DataCompet_driver_inter_factorized_GC_otf.csv')
tab_inter=tab_GC[1:500,]
tab_inter=na.exclude(tab_inter)
#tab_nointer=tab_GC[501:1000,]
#Pairwise
tp=sum(c((tab_inter$Pval_12_inter_GC_noexo<alpha_s)&(tab_inter$log_12_inter_noexo>threshold_s),((tab_inter$Pval_21_inter_GC_noexo<alpha_s)&(tab_inter$log_21_inter_noexo>threshold_s))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_GC_noexo>alpha_s)|(tab_inter$log_12_inter_noexo<threshold_s),((tab_inter$Pval_21_inter_GC_noexo>alpha_s)|(tab_inter$log_21_inter_noexo<threshold_s))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(3,rec,pch=16,col="blue")
#Conditional
tp=sum(c((tab_inter$Pval_12_inter_GC_exo<alpha_s)&(tab_inter$log_12_inter_exo>threshold_s),((tab_inter$Pval_21_inter_GC_exo<alpha_s)&(tab_inter$log_21_inter_exo>threshold_s))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_GC_exo>alpha_s)|(tab_inter$log_12_inter_exo<threshold_s),((tab_inter$Pval_21_inter_GC_exo>alpha_s)|(tab_inter$log_21_inter_exo<threshold_s))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(3,rec,pch=1,col="blue",lwd=2)

tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
tp=sum(c((tab_inter$Pval_12_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_12>0.1),((tab_inter$Pval_21_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_21>0.1))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_12<0.1),((tab_inter$Pval_21_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_21<0.1))))/(2*nrow(tab_inter))
rec=tp/(tp+fn)
points(3,rec,pch=18,col="red")

#legend("bottomright",c("GC (+pairwise)","GC conditional","CCM"),pch=c(16,1,18),col=c("black","black","red"),bty="n")


################NO INTER small dim
plot(0,0,xlab="Lyapunov exponent",ylab="Specificity",t="n",xlim=c(0.5,3.5),ylim=c(0,1))
#Chaos
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
fp=sum(c((tab_inter$Pval_12_noInter_GC<alpha_s)&(tab_inter$log_12_noInter>threshold_s),((tab_inter$Pval_21_noInter_GC<alpha_s)&(tab_inter$log_21_noInter>threshold_s))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_noInter_GC>alpha_s)|(tab_inter$log_12_noInter<threshold_s),((tab_inter$Pval_21_noInter_GC>alpha_s)|(tab_inter$log_21_noInter<threshold_s))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(1,spec,pch=16,col="blue")

fp=sum(c((tab_inter$Pval_12_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_noInter_v2>0.1),((tab_inter$Pval_21_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_noInter_v2>0.1))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_noInter_v2<0.1),((tab_inter$Pval_21_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_noInter_v2<0.1))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(1,spec,pch=18,col="red")

#Stochastic 2 species
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
fp=sum(c((tab_inter$Pval_12_noInter_GC<alpha_s)&(tab_inter$log_12_noInter>threshold_s),((tab_inter$Pval_21_noInter_GC<alpha_s)&(tab_inter$log_21_noInter>threshold_s))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_noInter_GC>alpha_s)|(tab_inter$log_12_noInter<threshold_s),((tab_inter$Pval_21_noInter_GC>alpha_s)|(tab_inter$log_21_noInter<threshold_s))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(2,spec,pch=16,col="blue")

fp=sum(c((tab_inter$Pval_12_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_noInter_v2>0.1),((tab_inter$Pval_21_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_noInter_v2>0.1))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_noInter_v2<0.1),((tab_inter$Pval_21_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_noInter_v2<0.1))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(2,spec,pch=18,col="red")

#Stochastic and driver
tab_GC=read.csv('../twoSpecies_andDriver/DataCompet_driver_inter_factorized_GC_otf.csv')
tab_inter=tab_GC[501:1000,]
tab_inter=na.exclude(tab_inter)
#Pairwise
fp=sum(c((tab_inter$Pval_12_inter_GC_noexo<alpha_s)&(tab_inter$log_12_inter_noexo>threshold_s),((tab_inter$Pval_21_inter_GC_noexo<alpha_s)&(tab_inter$log_21_inter_noexo>threshold_s))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_inter_GC_noexo>alpha_s)|(tab_inter$log_12_inter_noexo<threshold_s),((tab_inter$Pval_21_inter_GC_noexo>alpha_s)|(tab_inter$log_21_inter_noexo<threshold_s))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(3,spec,pch=16,col="blue")
#Conditional
fp=sum(c((tab_inter$Pval_12_inter_GC_exo<alpha_s)&(tab_inter$log_12_inter_exo>threshold_s),((tab_inter$Pval_21_inter_GC_exo<alpha_s)&(tab_inter$log_21_inter_exo>threshold_s))))/(2*nrow(tab_inter))
tn=sum(c((tab_inter$Pval_12_inter_GC_exo>alpha_s)|(tab_inter$log_12_inter_exo<threshold_s),((tab_inter$Pval_21_inter_GC_exo>alpha_s)|(tab_inter$log_21_inter_exo<threshold_s))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(3,spec,pch=1,col="blue",lwd=2)

tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_noIntersp1sp2factorized_CCM_otf.csv")
tp=sum(c((tab_inter$Pval_12_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_12>0.1),((tab_inter$Pval_21_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_21>0.1))))/(2*nrow(tab_inter))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_12<0.1),((tab_inter$Pval_21_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_21<0.1))))/(2*nrow(tab_inter))
spec=fp/(fp+tn)
points(3,spec,pch=18,col="red")


#######################################TO DO: LARGE DIMENSION

dev.off() 
