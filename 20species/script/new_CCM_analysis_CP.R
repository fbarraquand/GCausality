rm(list=ls())
graphics.off()

set.seed(5)

### Loading libraries
library(rEDM)

### Parameters
nsites = 25 ### how many samples of time series (sites or repeats)

### Useful variables
modelType = c("randomLV","randomVAR")
nrepeat = 1:nsites

nmodels = length(modelType)### For model and parameter types

alphaLevel = 0.2 # keep in mind this is corrected afterwards

### In this analysis we use only one network geometry and dimensionality
### Quantitative parameters and initial conditions for the VAR and Lotka-Volterra models are however varied

### Path to files with the time series and other useful data
pathrandomLV = "../data/Data_wTime_abs_LV.csv"
pathrandomVAR = "../data/Data_wTime_abs_VAR.csv"


######################################## Utilitary functions #############################################

path_to_file <- function(model) {
  eval(as.name(paste("path",model,sep="")))
}### use  strsplit for more fancy stuff


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

ccm_test = function(z,lag_order_inter){

  ### CCM Analysis 
  #species12=data.frame(1:nrow(z),exp(z)) #beware z is in log-scale
  species12=data.frame(1:nrow(z),z) #keeping the log
  names(species12)=c("time","sp1","sp2")


################################ Let's try stg completely different
  libsizes = 700 #Before, we only got up to 80, which is not fair because we consider 700 points in the GC test
  lm=length(libsizes)
  numsamples = 100
#####
sp1_xmap_sp2 <- ccm(species12, E = 1, lib_column = "sp1",target_column = "sp2", lib_sizes = libsizes, replace=FALSE,num_samples = 1)
rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp2"]=sample(species12[,"sp2"])
        sp1_xmap_sp2_random <- ccm(species_random, E = 1, lib_column = "sp1",target_column = "sp2", lib_sizes = libsizes, replace=FALSE,num_samples = 1)
        rho_dist[i]=sp1_xmap_sp2_random$rho
}
  Pval_1xmap2 = sum(rho_dist>sp1_xmap_sp2$rho)/numsamples 


sp2_xmap_sp1 <- ccm(species12, E = 1, lib_column = "sp2",target_column = "sp1", lib_sizes = libsizes, replace=FALSE,num_samples = 1)
rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = 1, lib_column = "sp2",target_column = "sp1", lib_sizes = libsizes, replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_2xmap1 = sum(rho_dist>sp2_xmap_sp1$rho)/numsamples 

  RhoLMax_12=sp2_xmap_sp1$rho# 1 causes 2 if 2 xmap 1
  RhoLMax_21=sp1_xmap_sp2$rho # 2 causes 1 if 1 xmap 2

  return(c(Pval_2xmap1,Pval_1xmap2,RhoLMax_12,RhoLMax_21)) ### NB we may find a way to output rho as well in a meaningful manner

}

pairwiseCCM <-function(x,alphaLevel,lagorder){ ### returns a matrix of causal links based on pairwise CCM
  ### Adjusted p-values with Benjamini-Hochberg correction
  nspecies=ncol(x) ## entry matrix has time on rows and species in columns
  p_value=matrix(0,nrow=nspecies,ncol=nspecies)
  rho_lmax=matrix(0,nrow=nspecies,ncol=nspecies)

  for (i in 1:nspecies){
    for (j in 1:nspecies){
      # cause first and effet later in grangertest()
      if (i >j){
        z=cbind(x[,i],x[,j])
        pccm=ccm_test(z,lagorder)
        p_value[i,j] = pccm[1]
        p_value[j,i] = pccm[2]
        rho_lmax[i,j] = pccm[3]
        rho_lmax[j,i] = pccm[4]
      }
    }
  }
  p_value_adj=p.adjust(p_value,method="BH")
  p_value_adj = matrix(p_value_adj,nrow=nspecies,ncol=nspecies)
  p=(p_value_adj<alphaLevel)
  return(list(p,rho_lmax))
}


################################ end of utilitary functions ##################################

#adjacency matrix for comparison to estimated matrices
interactions10species = rbind(c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,0,0,0,0,0,0,0),
                              c(1,1,1,1,1,0,0,0,0,0),
                              c(0,0,0,1,1,1,1,0,0,0),
                              c(0,0,0,0,1,1,1,0,0,0),
                              c(0,0,0,0,1,1,1,0,0,0),
                              c(0,0,0,0,0,0,0,1,1,1),
                              c(0,0,0,0,0,0,0,1,1,1),
                              c(0,0,0,0,0,0,0,1,1,1))
null_mat = matrix(0,10,10)
interaction_matrix = rbind(cbind(interactions10species,null_mat),cbind(null_mat,interactions10species))
### Adding some links between the two main compartments by adding a big connected cluster
interaction_matrix[8:13,8:13] = matrix(1,6,6) ## module filled with ones
causality_matrix = interaction_matrix

modelType

for (ksite in 1:nsites){ ### for sites or repeats
#for (ksite in 1:1){ ### for sites or repeats
  print(ksite)
  for (model in c("randomLV","randomVAR"))
    {

    ### Selects the files and then time series

    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    head(DB)
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps
    head(DB) ## 

    abundance_mat=as.matrix(DB[,4:23]) ### Create matrix with time series of abundances

    #### Pairwise CCM code 
   lag_order = 1 #lagOrder(abundance_mat)
    pCCM = pairwiseCCM(abundance_mat,alphaLevel,lag_order) ## alpha-level global higher than 5%
    rates2 = ratesClassif(pCCM[[1]],causality_matrix)
    resultsC2 = diagnosticsClassif(rates2)

    ### Output correlation matrix
    rhom=(pCCM[[2]]>0.25)
    rates = ratesClassif(rhom,causality_matrix)
    resultsC = diagnosticsClassif(rates)
  modelT= model
    nrepeat = ksite
    FPR = resultsC[1]
    TPR = resultsC[2]
    Precision = resultsC[3]
    newscoresClassif = data.frame(FPR,TPR,Precision,nrepeat,modelT)
    if ((ksite == 1)&(model=="randomLV"))
    {
      scoresClassif = newscoresClassif
    } else
    {
      scoresClassif = rbind(scoresClassif,newscoresClassif)
    }


    modelT= model
    nrepeat = ksite
    FPR = resultsC2[1]
    TPR = resultsC2[2]
    Precision = resultsC2[3]
    newscoresClassif2 = data.frame(FPR,TPR,Precision,nrepeat,modelT)
    if ((ksite == 1)&(model=="randomLV"))
    {
      scores_pCCM= newscoresClassif2
    } else
    {
      scores_pCCM = rbind(scores_pCCM,newscoresClassif2)
    }
  }
}


pdf(file = paste("../figures/ROC_CCM_BHcorrection_alpha",alphaLevel,"_newpval_log.pdf",sep=""),width=16,height = 8)
par(pty="s",mfrow=c(1,2),cex=1.5)

plot(scores_pCCM$FPR[scores_pCCM$modelT=="randomLV"],scores_pCCM$TPR[scores_pCCM$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise CCM with MY pval")
abline(a=0,b=1,lwd=2)
lines(scores_pCCM$FPR[scores_pCCM$modelT=="randomVAR"],scores_pCCM$TPR[scores_pCCM$modelT=="randomVAR"],pch=19,type="p",col="yellow")
legend("bottomright",legend=modelType,col=c("black","yellow"),pch=19,cex=0.8,bty="n")
#dev.off()

#pdf(file = paste("../figures/ROC_CCM_BHcorrection_rhobased.pdf",sep=""),width=16,height = 8)
#par(pty="s",mfrow=c(1,1),cex=1.5)

plot(scoresClassif$FPR[scoresClassif$modelT=="randomLV"],scoresClassif$TPR[scoresClassif$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise CCM with rho")
abline(a=0,b=1,lwd=2)
lines(scoresClassif$FPR[scoresClassif$modelT=="randomVAR"],scoresClassif$TPR[scoresClassif$modelT=="randomVAR"],pch=19,type="p",col="yellow")
#legend("right",legend=modelType,col=c("black","yellow"),pch=19,cex=0.8,bty="n")
dev.off()


