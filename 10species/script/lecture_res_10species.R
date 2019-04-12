graphics.off()
rm(list=ls())


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


        coln=c("site","sp1","sp2","E1","E2","pvalCobeyBaskerville","rhomax","deltarho","pvalCobeyBaskerville_adj","pvalCP","pvalCP_adj")
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
modelType = c("refLV","refVAR") #For now, we're just gonna focus on 10 species
mat_inter_per_site=array(0,dim=c(nspecies,nspecies,nsite,length(modelType)))
#val="pvalCP_adj"
val="pvalCP_adj"
pdf(paste("../figures/10species_example_CCM_",val,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(2,2,4,1))
m=0
for (model in modelType){
m=m+1
tab=read.csv(paste("../results/10species_CCM_per_interaction_",model,".csv",sep=""))
tab=na.exclude(tab)
 nsite=length(unique(tab$site))

nspecies=length(unique(tab$sp1))

mat_inter=matrix(NA,nspecies,nspecies)


plot(0,0,t="n",xlim=c(0.5,10.5),ylim=c(0.5,10.5),ylab="",xlab="",main=model)
for(i in 1:nspecies){
	for(j in 1:nspecies){
		if(i != j){
			id=which((tab$sp1==i)&(tab$sp2==j))
			mat_inter[j,i]=sum(tab[id,val]<alpha_level)/nsite
			if(interaction_matrix[j,i]==1){
				colo=rgb(0,0,1,mat_inter[j,i]) #Blue is right
			}else{
				colo=rgb(1,0,0,mat_inter[j,i]) #Red is wrong
			}
			points(j,i,col=colo,cex=5*mat_inter[j,i],pch=16)
			for(k in unique(tab$site)){
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

pdf(paste("../figures/ROC_pairwiseCCM_10species_",val,".pdf",sep=""),width=8,height=8)
colo=c("black","yellow")
plot(0,0,t="n",xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise CCM")
abline(a=0,b=1,lwd=2)
legend("right",legend=modelType,col=colo,pch=19,cex=0.8)

for(m in 1:2){
	for(k in unique(tab$site)){
		rates = ratesClassif(mat_inter_per_site[,,k-min(tab$site)+1,m],interaction_matrix)
	        resultsC = diagnosticsClassif(rates)
    modelT= model
    FPR = resultsC[1]
    TPR = resultsC[2]
    Precision = resultsC[3]
    points(FPR,TPR,pch=19,col=colo[m])
}
}
dev.off()
