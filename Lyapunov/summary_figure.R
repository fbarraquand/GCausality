####################### CP 13/06/2019
####################### This script uses the dynamics and lyapunov exponents designed and computed by Frederic Barraquand in order to plot the relation between lyapunov (degree of chaos) and recall+sensitivity of Granger Causality and Convergent Cross Mapping

rm(list=ls())
graphics.off()

library(Hmisc)

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
{ ### modified to avoid counting intrap. interactions
  if (typeof(estimated_mat)=="double")
  { Pos = (abs(estimated_mat)>0)} ## estimated links
  else if (typeof(estimated_mat)=="logical")
  {Pos = estimated_mat}
  else {print("Unknow type for estimated interaction matrix")}

  TruePos =  (abs(true_binary_mat)>0) ## work as well with non-binary
  TP = sum((TruePos[lower.tri(TruePos) | upper.tri(TruePos)] == TRUE) & (Pos[lower.tri(Pos) | upper.tri(Pos)] == TRUE))
  FN = sum((TruePos[lower.tri(TruePos) | upper.tri(TruePos)] == TRUE) & (Pos[lower.tri(Pos) | upper.tri(Pos)] == FALSE))
  FP = sum((TruePos[lower.tri(TruePos) | upper.tri(TruePos)] == FALSE) & (Pos[lower.tri(Pos) | upper.tri(Pos)] == TRUE))
  TN = sum((TruePos[lower.tri(TruePos) | upper.tri(TruePos)] == FALSE) & (Pos[lower.tri(Pos) | upper.tri(Pos)] == FALSE))
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


pdf("diagnostics_vs_lyapunov_pvalrplus1_newcolors.pdf",width=8)
layout(matrix(c(1,2,3,4),2,2,byrow=F))
par(mar=c(2,4,0.5,0.5))

#Small dim recall
plot(0,0,xlab="",ylab="Recall/Sensitivity",t="n",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n")
mtext("a)",side=2,las=2,at=1.,cex=0.75,line=3)
axis(1,c(1,2,3),rep("",3))
#Chaos
#tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
tab_inter=read.csv("../2species/results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
tp=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s),((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s))))
fn=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s),((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s))))

tp_strong=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s)))
tp_weak=sum(c((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s)))
fn_strong=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s)))
fn_weak=sum(c((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s)))

rec=tp/(tp+fn)
rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
points(1,rec,pch=16,col="blue",cex=2)
points(1,rec_strong,pch=16,col="darkblue",cex=1.25)
points(1,rec_weak,pch=16,col="cyan",cex=1.25)

tp=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1),((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1))))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1),((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1))))

tp_strong=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1)))
fn_strong=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1)))
tp_weak=sum(c((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1)))
fn_weak=sum(c((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1)))

rec=tp/(tp+fn)
rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
points(1.1,rec,pch=18,col="grey",cex=2)
points(1.1,rec_strong,pch=18,col="black",cex=1.25)
points(1.1,rec_weak,pch=18,col="gray90",cex=1.25)

#Stochastic 2 species
#tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
tab_inter=read.csv("../2species/results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
tp=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s),((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s))))
fn=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s),((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s))))
tp_strong=sum(c((tab_inter$Pval_21_inter_GC<alpha_s)&(tab_inter$log_21_inter>threshold_s)))
tp_weak=sum(c((tab_inter$Pval_12_inter_GC<alpha_s)&(tab_inter$log_12_inter>threshold_s)))
fn_strong=sum(c((tab_inter$Pval_21_inter_GC>alpha_s)|(tab_inter$log_21_inter<threshold_s)))
fn_weak=sum(c((tab_inter$Pval_12_inter_GC>alpha_s)|(tab_inter$log_12_inter<threshold_s)))

rec=tp/(tp+fn)
rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
points(2,rec,pch=16,col="blue",cex=2)
points(2,rec_strong,pch=16,col="darkblue",cex=1.25)
points(2,rec_weak,pch=16,col="cyan",cex=1.25)


tp=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1),((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1))))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1),((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1))))

tp_strong=sum(c((tab_inter$Pval_21_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_inter_v2>0.1)))
fn_strong=sum(c((tab_inter$Pval_21_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_inter_v2<0.1)))
tp_weak=sum(c((tab_inter$Pval_12_inter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_inter_v2>0.1)))
fn_weak=sum(c((tab_inter$Pval_12_inter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_inter_v2<0.1)))


rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
rec=tp/(tp+fn)
points(2.1,rec,pch=18,col="grey",cex=2)
points(2.1,rec_strong,pch=18,col="black",cex=1.25)
points(2.1,rec_weak,pch=18,col="grey90",cex=1.25)

#Stochastic and driver
#tab_GC=read.csv('../2species_driver/DataCompet_driver_inter_factorized_GC_otf.csv')
tab_GC=read.csv('../2species_driver/results/DataCompet_driver_intersp1sp2factorized_GC_otf_with_F_Wald_test.csv')
tab_inter=tab_GC[1:500,]
tab_inter=na.exclude(tab_inter)
#tab_nointer=tab_GC[501:1000,]
#Pairwise
tp=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest<alpha_s)&(tab_inter$log_12_inter_noexo>threshold_s),((tab_inter$Pval_21_inter_GC_no_exo_Ftest<alpha_s)&(tab_inter$log_21_inter_noexo>threshold_s))))
fn=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest>alpha_s)|(tab_inter$log_12_inter_noexo<threshold_s),((tab_inter$Pval_21_inter_GC_no_exo_Ftest>alpha_s)|(tab_inter$log_21_inter_noexo<threshold_s))))

tp_strong=sum(c((tab_inter$Pval_21_inter_GC_no_exo_Ftest<alpha_s)&(tab_inter$log_21_inter_noexo>threshold_s)))
tp_weak=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest<alpha_s)&(tab_inter$log_12_inter_noexo>threshold_s)))
fn_strong=sum(c((tab_inter$Pval_21_inter_GC_no_exo_Ftest>alpha_s)|(tab_inter$log_21_inter_noexo<threshold_s)))
fn_weak=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest>alpha_s)|(tab_inter$log_12_inter_noexo<threshold_s)))

rec=tp/(tp+fn)
rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
points(3,rec,pch=16,col="blue",cex=2)
points(3,rec_strong,pch=16,col="darkblue",cex=1.25)
points(3,rec_weak,pch=16,col="cyan",cex=1.25)

#Conditional
tp=sum(c((tab_inter$Pval_12_inter_GC_exo<alpha_s)&(tab_inter$log_12_inter_exo>threshold_s),((tab_inter$Pval_21_inter_GC_exo<alpha_s)&(tab_inter$log_21_inter_exo>threshold_s))))
fn=sum(c((tab_inter$Pval_12_inter_GC_exo>alpha_s)|(tab_inter$log_12_inter_exo<threshold_s),((tab_inter$Pval_21_inter_GC_exo>alpha_s)|(tab_inter$log_21_inter_exo<threshold_s))))
tp_strong=sum(c((tab_inter$Pval_21_inter_GC_exo<alpha_s)&(tab_inter$log_21_inter_exo>threshold_s)))
tp_weak=sum(c((tab_inter$Pval_12_inter_GC_exo<alpha_s)&(tab_inter$log_12_inter_exo>threshold_s)))
fn_strong=sum(c((tab_inter$Pval_21_inter_GC_exo>alpha_s)|(tab_inter$log_21_inter_exo<threshold_s)))
fn_weak=sum(c((tab_inter$Pval_12_inter_GC_exo>alpha_s)|(tab_inter$log_12_inter_exo<threshold_s)))

rec=tp/(tp+fn)
rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
points(3.05,rec,pch=1,col="blue",lwd=2,cex=2)
points(3.05,rec_strong,pch=1,col="darkblue",cex=1.25,lwd=1.25)
points(3.05,rec_weak,pch=1,col="cyan",cex=1.25,lwd=1.25)

#tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
tab_inter=read.csv("../2species_driver/results/DataCompet_driver_intersp1sp2factorized_CCM_otf_test.csv")
tp=sum(c((tab_inter$Pval_12_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_12>0.1),((tab_inter$Pval_21_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_21>0.1))))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_12<0.1),((tab_inter$Pval_21_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_21<0.1))))

tp_strong=sum(c((tab_inter$Pval_21_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_21>0.1)))
fn_strong=sum(c((tab_inter$Pval_21_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_21<0.1)))
tp_weak=sum(c((tab_inter$Pval_12_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_12>0.1)))
fn_weak=sum(c((tab_inter$Pval_12_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_12<0.1)))

rec_strong=tp_strong/(tp_strong+fn_strong)
rec_weak=tp_weak/(tp_weak+fn_weak)
rec=tp/(tp+fn)
points(3.1,rec,pch=18,col="grey",cex=2)
points(3.1,rec_strong,pch=18,col="black",cex=1.25)
points(3.1,rec_weak,pch=18,col="grey90",cex=1.25)

legend("bottomleft",c("GC (+pairwise)","GC conditional","CCM"),pch=c(16,1,18),col=c("blue","blue","grey"),bty="n")

################NO INTER small dim
plot(0,0,xlab="Lyapunov exponent",ylab="Specificity",t="n",xlim=c(0.5,3.5),ylim=c(0,1),xaxt="n")
mtext("c)",side=2,las=2,at=1.,cex=0.75,line=3)
axis(1,c(1,2,3),c("Deterministic chaos","2 spp","2 spp and driver"))
#Chaos
#tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")
tab_inter=read.csv("../2species/results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
fp=sum(c((tab_inter$Pval_12_noInter_GC<alpha_s)&(tab_inter$log_12_noInter>threshold_s),((tab_inter$Pval_21_noInter_GC<alpha_s)&(tab_inter$log_21_noInter>threshold_s))))
tn=sum(c((tab_inter$Pval_12_noInter_GC>alpha_s)|(tab_inter$log_12_noInter<threshold_s),((tab_inter$Pval_21_noInter_GC>alpha_s)|(tab_inter$log_21_noInter<threshold_s))))
spec=tn/(fp+tn)
points(1,spec,pch=16,col="blue",cex=2)

fp=sum(c((tab_inter$Pval_12_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_noInter_v2>0.1),((tab_inter$Pval_21_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_noInter_v2>0.1))))
tn=sum(c((tab_inter$Pval_12_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_noInter_v2<0.1),((tab_inter$Pval_21_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_noInter_v2<0.1))))
spec=tn/(fp+tn)
points(1.1,spec,pch=18,col="gray",cex=2)

#Stochastic 2 species
#tab_inter=read.csv("../comparisonGranger_vs_CCM_2species/results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")
tab_inter=read.csv("../2species/results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")
tab_inter=na.exclude(tab_inter)
fp=sum(c((tab_inter$Pval_12_noInter_GC<alpha_s)&(tab_inter$log_12_noInter>threshold_s),((tab_inter$Pval_21_noInter_GC<alpha_s)&(tab_inter$log_21_noInter>threshold_s))))
tn=sum(c((tab_inter$Pval_12_noInter_GC>alpha_s)|(tab_inter$log_12_noInter<threshold_s),((tab_inter$Pval_21_noInter_GC>alpha_s)|(tab_inter$log_21_noInter<threshold_s))))
spec=tn/(fp+tn)
points(2,spec,pch=16,col="blue",cex=2)

fp=sum(c((tab_inter$Pval_12_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_12_noInter_v2>0.1),((tab_inter$Pval_21_noInter_CCM_surr<alpha_s)&(tab_inter$RhoLMax_21_noInter_v2>0.1))))
tn=sum(c((tab_inter$Pval_12_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_12_noInter_v2<0.1),((tab_inter$Pval_21_noInter_CCM_surr>alpha_s)|(tab_inter$RhoLMax_21_noInter_v2<0.1))))
spec=tn/(fp+tn)
points(2.1,spec,pch=18,col="grey",cex=2)

#Stochastic and driver
#tab_GC=read.csv('../twoSpecies_andDriver/DataCompet_driver_inter_factorized_GC_otf.csv')
tab_GC=read.csv('../2species_driver/results/DataCompet_driver_intersp1sp2factorized_GC_otf_with_F_Wald_test.csv')
tab_inter=tab_GC[501:1000,]
tab_inter=na.exclude(tab_inter)
#Pairwise
fp=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest<alpha_s)&(tab_inter$log_12_inter_noexo>threshold_s),((tab_inter$Pval_21_inter_GC_no_exo_Ftest<alpha_s)&(tab_inter$log_21_inter_noexo>threshold_s))))
tn=sum(c((tab_inter$Pval_12_inter_GC_noexo_Ftest>alpha_s)|(tab_inter$log_12_inter_noexo<threshold_s),((tab_inter$Pval_21_inter_GC_no_exo_Ftest>alpha_s)|(tab_inter$log_21_inter_noexo<threshold_s))))
spec=tn/(fp+tn)
points(3,spec,pch=16,col="blue",cex=2)
#Conditional
fp=sum(c((tab_inter$Pval_12_inter_GC_exo<alpha_s)&(tab_inter$log_12_inter_exo>threshold_s),((tab_inter$Pval_21_inter_GC_exo<alpha_s)&(tab_inter$log_21_inter_exo>threshold_s))))
tn=sum(c((tab_inter$Pval_12_inter_GC_exo>alpha_s)|(tab_inter$log_12_inter_exo<threshold_s),((tab_inter$Pval_21_inter_GC_exo>alpha_s)|(tab_inter$log_21_inter_exo<threshold_s))))
spec=tn/(fp+tn)
points(3.05,spec,pch=1,col="blue",lwd=2,cex=2)

#tab_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_noIntersp1sp2factorized_CCM_otf_test.csv")
tab_inter=read.csv("../2species_driver/results/DataCompet_driver_noIntersp1sp2factorized_CCM_otf_test.csv")
tp=sum(c((tab_inter$Pval_12_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_12>0.1),((tab_inter$Pval_21_inter_CCM_surr_season<alpha_s)&(tab_inter$Rho_21>0.1))))
fn=sum(c((tab_inter$Pval_12_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_12<0.1),((tab_inter$Pval_21_inter_CCM_surr_season>alpha_s)|(tab_inter$Rho_21<0.1))))
spec=tn/(fp+tn)
points(3.1,spec,pch=18,col="grey",cex=2)


modelType = c("refLV","refVAR","randomLV","randomVAR")
interaction_matrix = rbind(c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,0,0,0,0,0,0,0),
                           c(1,1,1,1,0,0,0,0,0,0),
                           c(0,0,0,1,1,1,1,0,0,0),
                           c(0,0,0,0,1,1,1,0,0,0),
                           c(0,0,0,0,1,1,1,0,0,0),
                           c(0,0,0,0,0,0,0,1,1,1),
                           c(0,0,0,0,0,0,0,1,1,1),
                           c(0,0,0,0,0,0,0,1,1,1))


alpha_level=0.2
nsite=25

FPR=array(NA,dim=c(2,nsite,length(modelType)))
TPR=array(NA,dim=c(2,nsite,length(modelType)))

mat_inter_per_site=array(0,dim=c(20,20,nsite,length(modelType)))
val="p_pairwise_adj"
for (m in 1:length(modelType)){
        if(m>2){ #20 species, we need to duplicate the interaction matrix
                null_mat = matrix(0,10,10)
                interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
                interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
                causality_matrix = interaction_matrix_tmp
        }else{
                causality_matrix=interaction_matrix
        }

        model=modelType[m]
                tab=read.csv(paste("../10species/results/large_GC_per_interaction_",model,".csv",sep=""))
                tab=na.exclude(tab)
                nsite=length(unique(tab$site))
                nspecies=length(unique(tab$sp1))
                mat_inter=matrix(NA,nspecies,nspecies)

                nspecies=nrow(causality_matrix)
for(i in 1:nspecies){
        for(j in 1:nspecies){
                if(i != j){
                        id=which((tab$sp1==i)&(tab$sp2==j)) #Here, we don't need to correct for i/j because it was computed the right way in large_simulation_output (tab[i,j] is effect of j on i)
                        mat_inter[i,j]=sum(tab[id,val]<alpha_level)/nsite
                        for(k in unique(tab$site)){
                                id=which((tab$sp1==i)&(tab$sp2==j)&(tab$site==k))
                                if(val=='mat_simone'){ #If we're using from the simone-package, we did not store a p-value, but a magnitude of effect. When it is different from 0, j causes i
                                        if(abs(tab[id,val])>0){
                                                mat_inter_per_site[i,j,k-min(tab$site)+1,m]=1.
                                        }

                                }else{
                                        if(tab[id,val]<alpha_level){
                                                mat_inter_per_site[i,j,k-min(tab$site)+1,m]=1.
                                        }
                                }
                        }
                }
        }
}




        for(k in unique(tab$site)){
                rates = ratesClassif(mat_inter_per_site[1:nspecies,1:nspecies,k-min(tab$site)+1,m],causality_matrix)
                resultsC = diagnosticsClassif(rates)
    FPR[1,k,m] = resultsC[1]
    TPR[1,k,m] = resultsC[2]
}
}

#### CCM

mat_inter_per_site=array(0,dim=c(20,20,nsite,length(modelType)))
val="pvalCP_adj" #We can also use the p-value of the Cobey-Baskerville method, and ignore the BH-adjustment

m=0
for (model in modelType){
m=m+1
if(m<3){
tab=read.csv(paste("../10species/results/10species_CCM_per_interaction_",model,".csv",sep=""))
causality_matrix=interaction_matrix
}else{
tab=read.csv(paste("../20species/results/20species_CCM_per_interaction",model,"tot.csv",sep="_"))
#if(m==3){
#tab1=read.csv(paste("../20species/results/20species_CCM_per_interaction_",model,"_k1_k5.csv",sep=""))
#tab1=tab1[,2:ncol(tab1)]
#tab2=read.csv(paste("../20species/results/20species_CCM_per_interaction_",model,"_k6_k10.csv",sep=""))
#tab2=tab2[,2:ncol(tab2)]
#tab3=read.csv(paste("../20species/results/20species_CCM_per_interaction_",model,"_k11_k15.csv",sep=""))
#tab4=read.csv(paste("../20species/results/20species_CCM_per_interaction_",model,"_k16_k20.csv",sep=""))
#tab5=read.csv(paste("../20species/results/20species_CCM_per_interaction_",model,"_k21_k25.csv",sep=""))
#tab=rbind(tab1,tab2,tab3,tab4,tab5)
#}else{

                null_mat = matrix(0,10,10)
                interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
                interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
                causality_matrix = interaction_matrix_tmp

}
tab=na.exclude(tab)

 nsite=length(unique(tab$site))
nspecies=length(unique(tab$sp1))

mat_inter=matrix(NA,nspecies,nspecies)

for(i in 1:nspecies){
        for(j in 1:nspecies){
                if(i != j){
                        id=which((tab$sp1==i)&(tab$sp2==j)) #The code in analysis*, for CCM and GC, is written so that sp1 causes sp2, which is matrix_interaction[sp2,sp1]
                        mat_inter[j,i]=sum(tab[id,val]<alpha_level)/nsite
                        }
                        for(k in 1:nsite){
                                id=which((tab$sp1==i)&(tab$sp2==j)&(tab$site==k))
                                if(tab[id,val]<alpha_level){
                                        mat_inter_per_site[j,i,k-min(tab$site)+1,m]=1.
                                }
                        }
                }
        }
                nspecies=nrow(causality_matrix)
        for(k in unique(tab$site)){
                rates = ratesClassif(mat_inter_per_site[1:nspecies,1:nspecies,k-min(tab$site)+1,m],causality_matrix)
                resultsC = diagnosticsClassif(rates)
    FPR[2,k,m] = resultsC[1]
    TPR[2,k,m] = resultsC[2]
}
}


##Recall
#modelType = c("refLV","refVAR","randomLV","randomVAR")
#modelType=c("refLV","randomLV","refVAR","randomVAR")
list_m=c(1,3,2,4)
tmp_m=0
list_pch=c(16,18,16,18)
plot(0,0,t="n",xlim=c(0.5,4.5),ylim=c(0,1),xlab="",ylab="",xaxt="n")
mtext("b)",side=2,las=2,at=1.,cex=0.75,line=3)
axis(1,1:4,rep("",4))
for(m in list_m){
	tmp_m=tmp_m+1
	mean_TPR=mean(TPR[1,,m])
	sd_TPR=sd(TPR[1,,m])
	errbar(tmp_m,mean_TPR,mean_TPR+sd_TPR,mean_TPR-sd_TPR,add=T,pch=16,col="blue",errbar.col="blue",cex=2)
	mean_TPR=mean(TPR[2,,m])
	sd_TPR=sd(TPR[2,,m])
	errbar(tmp_m+0.15,mean_TPR,mean_TPR+sd_TPR,mean_TPR-sd_TPR,add=T,pch=18,col="grey",errbar.col="grey",cex=2)
}

tmp_m=0
plot(0,0,t="n",xlim=c(0.5,4.5),ylim=c(0,1),xlab="",ylab="",xaxt="n")
mtext("d)",side=2,las=2,at=1.,cex=0.75,line=3)
axis(1,1:4,c("Ricker 10","Ricker 20","MAR 10","MAR 20"),las=1)
for(m in list_m){
	tmp_m=tmp_m+1
	mean_FPR=mean(1-FPR[1,,m])
	sd_FPR=sd(1-FPR[1,,m])
	errbar(tmp_m,mean_FPR,mean_FPR+sd_FPR,mean_FPR-sd_FPR,add=T,col="blue",pch=16,errbar.col="blue",cex=2)
	mean_FPR=mean(1-FPR[2,,m])
	sd_FPR=sd(1-FPR[2,,m])
	errbar(tmp_m+0.15,mean_FPR,mean_FPR+sd_FPR,mean_FPR-sd_FPR,add=T,col="grey",pch=18,errbar.col="grey",cex=2)
}


legend("bottomright",c("GC","CCM"),pch=c(16,18),col=c("blue","grey"),bty="n",cex=1)


dev.off() 
