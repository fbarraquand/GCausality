########################################################################################################
### Read results from analysis_CCM_byCP.R to compute diagnostics results (false positives and so on, to draw ROC plot)
### From FBarraquand, analysis_Simone_Clustering_woutIntraSp.R 
###CP July 21st 2019: cleaned up code to have everything on panels
###CP June 2020: corrected color code to make the figure easier to read
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
modelType = c("refLV","refVAR","randomLV","randomVAR")
nice_modelType=c("Lotka-Volterra, 10sp","MAR, 10sp","Lotka-Volterra, 20sp","MAR, 20sp")
val_list=c("p_pairwise_adj","mat_simone","pvalCobeyBaskerville_adj","pvalCP_adj")
#val_list=c("p_pairwise_adj","pvalCP_adj")
nice_val=c("GC pairwise","GC SIMoNe","CCM CB2016","CCM permutation")
type_list=c("GC","GC","CCM","CCM")
margin=c("a)","b)","c)","d)")
#colo=c("black","yellow","blue","red")
colo=c("grey","yellow","black","orange")
pch_vec=c(16,17,16,17)

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
nsite=25

pdf("../figures/ROC_bestpval_formain_newsetofcolors.pdf")
par(mfrow=c(2,2))
for(v in 1:length(val_list)){
	val=val_list[v]
                if(v==1){
                        par(mar=c(2,4.5,4,1))
                }else if(v==2){
                        par(mar=c(2,2,4,1))
                }else if(v==3){
                        par(mar=c(4,4.5,2,1))
                }else if(v==4){
                        par(mar=c(4,2,2,1))
                }

	if(v>2){
		axlab="False Positive Rate (1 - specificity)"
	}else{
		axlab=""
	}
	if(v%%2==1){
		aylab="True Positive Rate (recall/sensitivity)"
	}else{
		aylab=""
	}
	plot(0,0,t="n",xlim=c(0,1),ylim=c(0,1),xlab =axlab,ylab =aylab, main = nice_val[v])
	abline(a=0,b=1,lwd=2)
	if(v==1){
		legend("bottomright",legend=nice_modelType,col=colo,pch=pch_vec,cex=1,bty="n")
	}
	mtext(margin[v],side=2,las=2,at=1.1,line=1.5,cex=0.8)


	mat_inter_per_site=array(0,dim=c(20,20,nsite,length(modelType)))

	for(m in 1:length(modelType)){
	        if(m>2){ #20 species, we need to duplicate the interaction matrix
        	        null_mat = matrix(0,10,10)
                	interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
                	interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
               		causality_matrix = interaction_matrix_tmp
        	}else{
                	causality_matrix=interaction_matrix
        	}

        	model=modelType[m]
		if(type_list[v]=="GC"){
                	tab=read.csv(paste("../results/large_GC_per_interaction_",model,".csv",sep=""))
		}else{
			if(m<3){
				tab=read.csv(paste("../results/10species_CCM_per_interaction_",model,".csv",sep=""))
			}else{
				tab=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_tot.csv",sep=""))
			}
		}
		tab=na.exclude(tab)
		nb_sp=length(unique(tab$sp1))


		for(i in 1:nb_sp){
        		for(j in 1:nb_sp){
                		if(i != j){
		                        for(k in 1:nsite){
                                		id=which((tab$sp1==i)&(tab$sp2==j)&(tab$site==k))
        					if(type_list[v]=="GC"){
                                			if(val=='mat_simone'){ 
                                        			if(abs(tab[id,val])>0){
                                                			mat_inter_per_site[i,j,k-min(tab$site)+1,m]=1.
                                        			}

                                			}else{
                                        			if(tab[id,val]<alpha_level){
                                                			mat_inter_per_site[i,j,k-min(tab$site)+1,m]=1.
                                        			}
                                			}
						}else{
                                			if(tab[id,val]<alpha_level){
                                        			mat_inter_per_site[j,i,k-min(tab$site)+1,m]=1.
                                			}
						}
					}#end k
				}
			}#end j
		}#end i

	for(k in 1:nsite){
                rates = ratesClassif(mat_inter_per_site[1:nb_sp,1:nb_sp,k-min(tab$site)+1,m],causality_matrix)
                resultsC = diagnosticsClassif(rates)
    		FPR = resultsC[1]
    		TPR = resultsC[2]
    		Precision = resultsC[3]
    		points(FPR,TPR,col=colo[m],pch=pch_vec[m])
	}#end k
	} #end model

}#nd val

dev.off()
