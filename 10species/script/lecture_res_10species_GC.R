### CP April 2019
### Read results from large_simulation_output_per_interaction to compute diagnostics results (false positives and so on, to draw ROC plot) for Granger-causality, as well as results for each interaction

graphics.off()
rm(list=ls())

modelType = c("refLV","refVAR","randomLV","randomVAR")

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

mat_inter_per_site=array(0,dim=c(20,20,nsite,length(modelType)))
val="p_pairwise_adj"
#val="mat_simone"

pdf(paste("../figures/large_example_GC_",val,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2),mar=c(2,2,4,1))

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
		tab=read.csv(paste("../results/large_GC_per_interaction_",model,".csv",sep=""))
		tab=na.exclude(tab)
 		nsite=length(unique(tab$site))
		nspecies=length(unique(tab$sp1))
		mat_inter=matrix(NA,nspecies,nspecies)


plot(0,0,t="n",xlim=c(0.5,nspecies+.5),ylim=c(0.5,nspecies+.5),ylab="",xlab="",main=model)
for(i in 1:nspecies){
        for(j in 1:nspecies){
                if(i != j){
                        id=which((tab$sp1==i)&(tab$sp2==j)) #Here, we don't need to correct for i/j because it was computed the right way in large_simulation_output (tab[i,j] is effect of j on i)
			if(val=='mat_simone'){ #If we're using from the simone-package, we did not store a p-value, but a magnitude of effect. When it is different from 0, j causes i
				mat_inter[i,j]=sum(abs(tab[id,val])>0)/nsite
                        }else{
				mat_inter[i,j]=sum(tab[id,val]<alpha_level)/nsite
			}
			if(mat_inter[i,j]>0){
	                        if(causality_matrix[i,j]==1){
        	                        colo=rgb(0,0,1,1) #Blue is true positive
                	        }else{
                        	        colo=rgb(1,0,0,1) #Red is false positive
                        	}
                        	points(i,j,col=colo,cex=5*mat_inter[i,j],pch=16)
                        }else{
                                if(causality_matrix[i,j]==1){ #false negatives
                                        points(i,j,col="black",cex=2.5,pch=16)
				}
                        } 
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



}

dev.off()

pdf(paste("../figures/ROC_pairwiseGC_large_",val,".pdf",sep=""),width=8,height=8)
colo=c("black","yellow","blue","red")
plot(0,0,t="n",xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise GC")
abline(a=0,b=1,lwd=2)
legend("right",legend=modelType,col=colo,pch=19,cex=0.8)

for(m in 1:4){
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
    FPR = resultsC[1]
    TPR = resultsC[2]
    Precision = resultsC[3]
    points(FPR,TPR,pch=19,col=colo[m])
}
}
dev.off()

