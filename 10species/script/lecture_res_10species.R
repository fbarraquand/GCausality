########################################################################################################
### Read results from analysis_CCM_byCP.R to compute diagnostics results (false positives and so on, to draw ROC plot),as well as results for each interaction
### From FBarraquand, analysis_Simone_Clustering_woutIntraSp.R 
###CP April 2019
########################################################################################################

graphics.off()
rm(list=ls())

########################### Utilitary functions######################
###########  Written by FB
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

############################## End of utilitary functions

#Inits
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
nspecies=nrow(interaction_matrix)
nsite=25


modelType = c("refLV","refVAR","randomLV","randomVAR") #randomVAR running on Fred's
mat_inter_per_site=array(0,dim=c(20,20,nsite,length(modelType)))
val="pvalCP_adj" #We can also use the p-value of the Cobey-Baskerville method, and ignore the BH-adjustment
#val="pvalCobeyBaskerville_adj" #We can also use the p-value of the Cobey-Baskerville method, and ignore the BH-adjustment

pdf(paste("../figures/large_example_CCM_",val,"_with_legend.pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(2,2,4,1),xpd=TRUE)
m=0

for (model in modelType){
m=m+1
if(m<3){
tab=read.csv(paste("../results/10species_CCM_per_interaction_",model,".csv",sep=""))
causality_matrix=interaction_matrix
}else{
if(m==3){
tab1=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_k1_k5.csv",sep=""))
tab1=tab1[,2:ncol(tab1)]
tab2=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_k6_k10.csv",sep=""))
tab2=tab2[,2:ncol(tab2)]
tab3=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_k11_k15.csv",sep=""))
tab4=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_k16_k20.csv",sep=""))
tab5=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_k21_k25.csv",sep=""))
tab=rbind(tab1,tab2,tab3,tab4,tab5)
}else{
tab=read.csv("../../20species/results/20species_CCM_per_interaction_randomVAR.csv")
}

                null_mat = matrix(0,10,10)
                interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
                interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
                causality_matrix = interaction_matrix_tmp

}
tab=na.exclude(tab)

 nsite=length(unique(tab$site))
nspecies=length(unique(tab$sp1))

mat_inter=matrix(NA,nspecies,nspecies)

#For now, I have removed the transparency, which was nice but sometimes misleading
plot(0,0,t="n",xlim=c(0.5,nspecies+0.5),ylim=c(0.5,nspecies+0.5),ylab="",xlab="",main="")
if(m<3){
points(0.8,0.8,col="black",bg="white",cex=5,pch=21)
points(0.8,0.8,col="black",bg="white",cex=0.5,pch=21)

lines(c(0.7,0.9),c(-0.3,-0.3),lty=1)
lines(c(0.7,0.7),c(0.8,-0.3),lty=3)
lines(c(0.9,0.9),c(0.8,-0.3),lty=3)
text(0.25,-0.3,"10%",pos=2,cex=0.5)

lines(c(0.3,1.3),c(-0.75,-0.75),lty=1)
lines(c(0.3,0.3),c(0.8,-0.75),lty=3)
lines(c(1.3,1.3),c(0.8,-0.75),lty=3)
text(0.25,-0.75,"100%",pos=2,cex=0.5)

}else{
points(0.8,0.8,col="black",bg="white",cex=2.5,pch=21)
points(0.8,0.8,col="black",bg="white",cex=0.5,pch=21)

lines(c(0.7,0.9),c(-0.9,-0.9),lty=1)
lines(c(0.7,0.7),c(0.8,-0.9),lty=3)
lines(c(0.9,0.9),c(0.8,-0.9),lty=3)
text(0.95,-0.9,"10%",pos=4,cex=0.5)

lines(c(0.3,1.3),c(-1.5,-1.5),lty=1)
lines(c(0.3,0.3),c(0.8,-1.5),lty=3)
lines(c(1.3,1.3),c(0.8,-1.5),lty=3)
text(0.95,-1.5,"50%",pos=4,cex=0.5)
}

for(i in 1:nspecies){
	for(j in 1:nspecies){
		if(i != j){
			id=which((tab$sp1==i)&(tab$sp2==j)) #The code in analysis*, for CCM and GC, is written so that sp1 causes sp2, which is matrix_interaction[sp2,sp1]
			mat_inter[j,i]=sum(tab[id,val]<alpha_level)/nsite
			if(mat_inter[j,i]>0){
				if(causality_matrix[j,i]==1){
					#colo=rgb(0,0,1,mat_inter[j,i]) #Blue is right, true positives
					colo=rgb(0,0,1,1) #Blue is right, true positives
				}else{
					#colo=rgb(1,0,0,mat_inter[j,i]) #Red is wrong, false positives
					colo=rgb(1,0,0,1) #Red is wrong, false positives
				}
				points(j,i,col=colo,cex=5*mat_inter[j,i],pch=16)
			}else{
				if(causality_matrix[j,i]==1){ #false negatives
					print(mat_inter[j,i])
					points(j,i,col="black",cex=2.5,pch=16)
				}
			}
			for(k in 1:nsite){
				id=which((tab$sp1==i)&(tab$sp2==j)&(tab$site==k))
				if(tab[id,val]<alpha_level){
					mat_inter_per_site[j,i,k-min(tab$site)+1,m]=1.
				}
			}
		}
	}
}



}

dev.off()

### ROC curve
pdf(paste("../figures/ROC_pairwiseCCM_large_",val,".pdf",sep=""),width=8,height=8)
par(cex=1.5)
colo=c("black","yellow","blue","red")
plot(0,0,t="n",xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise CCM")
abline(a=0,b=1,lwd=2)
legend("right",legend=modelType,col=colo,pch=19,cex=0.8,bty="n")

for(m in 1:length(modelType)){
 if(m>2){
                null_mat = matrix(0,10,10)
                interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
                interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
                causality_matrix = interaction_matrix_tmp
        }else{
                causality_matrix=interaction_matrix
        }
                nspecies=nrow(causality_matrix)
        for(k in unique(tab$site)){
                rates = ratesClassif(mat_inter_per_site[1:nspecies,1:nspecies,k-min(tab$site)+1,m],causality_matrix)
                resultsC = diagnosticsClassif(rates)
    modelT= model
    FPR = resultsC[1]
    TPR = resultsC[2]
    Precision = resultsC[3]
    points(FPR,TPR,pch=19,col=colo[m])
}
}
dev.off()
