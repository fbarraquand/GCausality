########################################################################################################
### Read results from analysis_CCM_10species.R to compute diagnostics results (false positives and so on, to draw ROC plot),as well as results for each interaction
### From FBarraquand, analysis_Simone_Clustering_woutIntraSp.R, and CP lecture_res_10species and lecture_res_10species_GC
### CP 21st July 2019: Plot 4 panels at each time, for the methods we want
########################################################################################################

graphics.off()
rm(list=ls())


nb_sp=20

if(nb_sp==10){
modelType = c("refLV","refVAR")#,"randomLV","randomVAR")
}else{
modelType = c("randomLV","randomVAR")
}
nice_modelType=c("Lotka-Volterra","MAR","","")
#val_list=c("p_pairwise_adj","mat_simone","pvalCobeyBaskerville_adj","pvalCP_adj")

#pvalCobeyBaskerville_adj
#mat_simone

val_list=c("p_pairwise_adj","pvalCP_adj")
nice_val=c("GC pairwise","","CCM permutation","")
#val_list=c("mat_simone","pvalCobeyBaskerville_adj")
#nice_val=c("GC SIMoNe","","CCM CobeyBaskerville","")
type_list=c("GC","CCM")
margin=c("a)","b)","c)","d)")

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

pdf(paste("../figures/",nb_sp,"sp_per_inter_GC_CCM.pdf",sep=""),width=10,height=10)
par(mfrow=c(2,2),xpd=TRUE)
v=0
lab=0
for(val in val_list){
	v=v+1
	m=0
	for (model in modelType){
		m=m+1
		lab=lab+1
		if(lab==1){
			par(mar=c(2,4.5,4,1))
		}else if(lab==2){
			par(mar=c(2,2,4,1))
		}else if(lab==3){
			par(mar=c(4,4.5,2,1))
		}else if(lab==4){
			par(mar=c(4,2,2,1))
		}

		if(nb_sp==10){
			if(type_list[v]=="GC"){
                		tab=read.csv(paste("../results/large_GC_per_interaction_",model,".csv",sep=""))
			}else{
				tab=read.csv(paste("../results/10species_CCM_per_interaction_",model,".csv",sep=""))
			}
			causality_matrix=interaction_matrix
		}else{
			if(type_list[v]=="GC"){
                		tab=read.csv(paste("../results/large_GC_per_interaction_",model,".csv",sep=""))
			}else{
				tab=read.csv(paste("../../20species/results/20species_CCM_per_interaction_",model,"_tot.csv",sep=""))

			}
                	null_mat = matrix(0,10,10)
	                interaction_matrix_tmp = rbind(cbind(interaction_matrix,null_mat),cbind(null_mat,interaction_matrix))
        	        interaction_matrix_tmp[8:13,8:13] = matrix(1,6,6) ## module filled with ones
                	causality_matrix = interaction_matrix_tmp

		}
	tab=na.exclude(tab)
	nsite=length(unique(tab$site))
	mat_inter=matrix(NA,nb_sp,nb_sp)
	
	plot(0,0,t="n",xlim=c(0.5,nb_sp+0.5),ylim=c(0.5,nb_sp+0.5),ylab=nice_val[lab],xlab="",main=nice_modelType[lab],cex.lab=1.5)
	mtext(margin[lab],side=2,las=2,at=(nb_sp+0.75),line=1.75)

if(m==1&v==2){
if(nb_sp==10){
points(0.8,0.8,col="black",bg="white",cex=5,pch=21)
points(0.8,0.8,col="black",bg="white",cex=0.5,pch=21)

lines(c(0.74,0.86),c(-0.3,-0.3),lty=1)
lines(c(0.74,0.74),c(0.8,-0.3),lty=3)
lines(c(0.86,0.86),c(0.8,-0.3),lty=3)
text(0.25,-0.3,"10%",pos=2,cex=0.75)

lines(c(0.32,1.28),c(-0.75,-0.75),lty=1)
lines(c(0.32,0.32),c(0.8,-0.75),lty=3)
lines(c(1.28,1.28),c(0.8,-0.75),lty=3)
text(0.25,-0.75,"100%",pos=2,cex=0.75)
}else{
points(0.8,0.8,col="black",bg="white",cex=2.5,pch=21)
points(0.8,0.8,col="black",bg="white",cex=0.5,pch=21)

lines(c(0.7,0.9),c(-0.9,-0.9),lty=1)
lines(c(0.7,0.7),c(0.8,-0.9),lty=3)
lines(c(0.9,0.9),c(0.8,-0.9),lty=3)
text(0.95,-0.9,"10%",pos=4,cex=0.75)

lines(c(0.3,1.3),c(-1.5,-1.5),lty=1)
lines(c(0.3,0.3),c(0.8,-1.5),lty=3)
lines(c(1.3,1.3),c(0.8,-1.5),lty=3)
text(0.95,-1.5,"50%",pos=4,cex=0.75)
}
}
for(i in 1:nb_sp){
        for(j in 1:nb_sp){
                if(i != j){
                        id=which((tab$sp1==i)&(tab$sp2==j)) #Here, we don't need to correct for i/j because it was computed the right way in large_simulation_output (tab[i,j] is effect of j on i)
			if(type_list[v]=="GC"){
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

			}else if(type_list[v]=="CCM"){
                        	mat_inter[j,i]=sum(tab[id,val]<alpha_level)/nsite
                        	if(mat_inter[j,i]>0){
                                	if(causality_matrix[j,i]==1){
                                        	#colo=rgb(0,0,1,mat_inter[j,i]) #Blue is right, true positives #I used to have some transparency here, but let's ignore it
                                        	colo=rgb(0,0,1,1) #Blue is right, true positives
                                	}else{
	                                        #colo=rgb(1,0,0,mat_inter[j,i]) #Red is wrong, false positives
        	                                colo=rgb(1,0,0,1) #Red is wrong, false positives
                	                }
                        	        points(j,i,col=colo,cex=5*mat_inter[j,i],pch=16)
                        	}else{
                                	if(causality_matrix[j,i]==1){ #false negatives
                                        	points(j,i,col="black",cex=2.5,pch=16)
                                	}
                        	}
                	}
        	}
	}
}


}
}

dev.off()

