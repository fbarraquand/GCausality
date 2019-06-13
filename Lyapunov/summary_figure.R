####################### CP 13/06/2019
####################### This script uses the dynamics and lyapunov exponents designed and computed by Frederic Barraquand in order to plot the relation between lyapunov (degree of chaos) and recall+sensitivity of Granger Causality and Convergent Cross Mapping

rm(list=ls())
graphics.off()

#Small dimension
l_chaos=0.41
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


#pdf("diagnostics_vs_lyapunov.pdf")
pdf("total_recall_smalldim.pdf")
#par(mfrow=c(2,2))

#Small dim
plot(0,0,xlab="Lyapunov exponent",ylab="Recall",t="n",xlim=c(-0.5,0.5),ylim=c(0,1))
#Chaos
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
	tp=((tab_inter$Pval_12_inter_GC[n]<alpha_s)&(tab_inter$log_12_inter[n]>threshold_s))+((tab_inter$Pval_21_inter_GC[n]<alpha_s)&(tab_inter$log_21_inter[n]>threshold_s))
	fn=((tab_inter$Pval_12_inter_GC[n]>alpha_s)|(tab_inter$log_12_inter[n]<threshold_s))+((tab_inter$Pval_21_inter_GC[n]>alpha_s)|(tab_inter$log_21_inter[n]<threshold_s))
	rec_tmp=max(0,tp/(tp+fn))
	if(rec_tmp==0){
		rec[1]=rec[1]+1
	}else if(rec_tmp==0.5){
		rec[2]=rec[2]+1
	}else{
		rec[3]=rec[3]+1
	}
}
points(rep(l_chaos,3),c(0,0.5,1),pch=16,cex=5*rec/n)

rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_CCM_surr[n]<alpha_s)&(tab_inter$RhoLMax_12_inter_v2[n]>0.1))+((tab_inter$Pval_21_inter_CCM_surr[n]<alpha_s)&(tab_inter$RhoLMax_21_inter_v2[n]>0.1))
        fn=((tab_inter$Pval_12_inter_CCM_surr[n]>alpha_s)|(tab_inter$RhoLMax_12_inter_v2[n]<0.1))+((tab_inter$Pval_21_inter_CCM_surr[n]>alpha_s)|(tab_inter$RhoLMax_21_inter_v2[n]<0.1))
        rec_tmp=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
points(rep(l_chaos,3),c(0,0.5,1),pch=18,cex=3*rec/n,col="red")


#Stochastic 2 species
tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
#tab_nointer=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")
rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_GC[n]<alpha_s)&(tab_inter$log_12_inter[n]>threshold_s))+((tab_inter$Pval_21_inter_GC[n]<alpha_s)&(tab_inter$log_21_inter[n]>threshold_s))
        fn=((tab_inter$Pval_12_inter_GC[n]>alpha_s)|(tab_inter$log_12_inter[n]<threshold_s))+((tab_inter$Pval_21_inter_GC[n]>alpha_s)|(tab_inter$log_21_inter[n]<threshold_s))
        rec_tmp=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
	points(rep(l_2species,3),c(0,0.5,1),pch=16,cex=5*rec/n)

rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_CCM_surr[n]<alpha_s)&(tab_inter$RhoLMax_12_inter_v2[n]>0.1))+((tab_inter$Pval_21_inter_CCM_surr[n]<alpha_s)&(tab_inter$RhoLMax_21_inter_v2[n]>0.1))
        fn=((tab_inter$Pval_12_inter_CCM_surr[n]>alpha_s)|(tab_inter$RhoLMax_12_inter_v2[n]<0.1))+((tab_inter$Pval_21_inter_CCM_surr[n]>alpha_s)|(tab_inter$RhoLMax_21_inter_v2[n]<0.1))
        rec_tmp=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
points(rep(l_chaos,3),c(0,0.5,1),pch=18,cex=3*rec/n,col="red")



#Stochastic and driver
tab_GC=read.csv('../twoSpecies_andDriver/DataCompet_driver_inter_factorized_GC_otf.csv')
tab_inter=tab_GC[1:500,]
#tab_nointer=tab_GC[501:1000,]
#Pairwise
rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_GC_noexo[n]<alpha_s)&(tab_inter$log_12_inter_noexo[n]>threshold_s))+((tab_inter$Pval_21_inter_GC_noexo[n]<alpha_s)&(tab_inter$log_21_inter_noexo[n]>threshold_s))
        fn=((tab_inter$Pval_12_inter_GC_noexo[n]>alpha_s)|(tab_inter$log_12_inter_noexo[n]<threshold_s))+((tab_inter$Pval_21_inter_GC_noexo[n]>alpha_s)|(tab_inter$log_21_inter_noexo[n]<threshold_s))
        rec_=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
	points(rep(l_driver,3),c(0,0.5,1),pch=1,cex=5*rec/n)
#Conditional
rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_GC_exo[n]<alpha_s)&(tab_inter$log_12_inter_exo[n]>threshold_s))+((tab_inter$Pval_21_inter_GC_exo[n]<alpha_s)&(tab_inter$log_21_inter_exo[n]>threshold_s))
        fn=((tab_inter$Pval_12_inter_GC_exo[n]>alpha_s)|(tab_inter$log_12_inter_exo[n]<threshold_s))+((tab_inter$Pval_21_inter_GC_exo[n]>alpha_s)|(tab_inter$log_21_inter_exo[n]<threshold_s))
        rec_tmp=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
	points(rep(l_driver,3),c(0,0.5,1),pch=16,cex=5*rec/n)

tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
rec=rep(0,3)
for(n in 1:nrow(tab_inter)){
        tp=((tab_inter$Pval_12_inter_CCM_surr_season[n]<alpha_s)&(tab_inter$Rho_12[n]>0.1))+((tab_inter$Pval_21_inter_CCM_surr_season[n]<alpha_s)&(tab_inter$Rho_21[n]>0.1))
        fn=((tab_inter$Pval_12_inter_CCM_surr_season[n]>alpha_s)|(tab_inter$Rho_12[n]<0.1))+((tab_inter$Pval_21_inter_CCM_surr_season[n]>alpha_s)|(tab_inter$Rho_21[n]<0.1))
        rec_tmp=max(0,tp/(tp+fn))
        if(rec_tmp==0){
                rec[1]=rec[1]+1
        }else if(rec_tmp==0.5){
                rec[2]=rec[2]+1
        }else{
                rec[3]=rec[3]+1
        }
}
        points(rep(l_driver,3),c(0,0.5,1),pch=18,cex=3*rec/n,col="red")

legend("bottomright",c("GC (+pairwise)","GC conditional","CCM"),pch=c(16,1,18),col=c("black","black","red"),bty="n")

dev.off()

if(1==0){
plot(0,0,xlab="Lyapunov exponent",ylab="Specificity",t="n",xlim=c(-0.5,0.5),ylim=c(0,1))
tab_nointer=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")
for(n in 1:nrow(tab_nointer)){
	fp=((tab_nointer$Pval_12_noInter_GC[n]<alpha_s)&(tab_nointer$log_12_noInter[n]>threshold_s))+((tab_nointer$Pval_21_noInter_GC[n]<alpha_s)&(tab_nointer$log_21_noInter[n]>threshold_s))
	tn=((tab_nointer$Pval_12_noInter_GC[n]>alpha_s)|(tab_nointer$log_12_noInter[n]<threshold_s))+((tab_nointer$Pval_21_noInter_GC[n]>alpha_s)|(tab_nointer$log_21_noInter[n]<threshold_s))
	spec=fp/(fp+tn)
	points(l_chaos,spec)
}



tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
tab_nointer=read.csv("../twoSpecies_andDriver/DataCompet_driver_noIntersp1sp2factorized_CCM_otf.csv")

} 
