##############################################################
### FB 15/06/2018 -- Analysis for reference parameter set ####
### Uses Simone and pairwise GC [VAR(p)] inference        ####
### Model selection is done by BIC                        ####
### 18/06/2018 Corrected the FDR BH correction (matrix output)
##############################################################

rm(list = ls())
set.seed(5)

### Loading libraries
library(simone)

### Parameters
nsites = 25 ### how many samples of time series (sites or repeats)

### Useful variables
modelType = c("randomLV","randomVAR")
nrepeat = 1:nsites

nmodels = length(modelType)### For model and parameter types

alphaLevel = 0.2

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

library(vars)

lagOrder <-function(y,lag.max=15){ #lag-order as determined by full VAR model
  yd=data.frame(abundance_mat)
  varresult = VAR(y=yd, type="none",lag.max=lag.max)
  return(as.numeric(varresult$p))
}

pairwiseGC <-function(x,alphaLevel,lagorder){ ### returns a matrix of causal links based on pairwise GC
  ### Adjusted p-values with Benjamini-Hochberg correction
  nspecies=ncol(x) ## entry matrix has time on rows and species in columns
  p_value=matrix(0,nrow=nspecies,ncol=nspecies)
  
  for (i in 1:nspecies){
    for (j in 1:nspecies){
      # cause first and effet later in grangertest()
      if (i !=j){
        g=grangertest(x[,j],x[,i],order=lagorder)
        p_value[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
      }
    }
  }
  p_value_adj=p.adjust(p_value,method="BH")
  p_value_adj = matrix(p_value_adj,nrow=nspecies,ncol=nspecies)
  pGC=(p_value_adj<alphaLevel)
  return(pGC) 
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
  
  for (model in c("randomLV","randomVAR"))
    {
    
    ### Selects the files and then time series
  
    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    head(DB)
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps
    head(DB) ## 
  
    abundance_mat=as.matrix(DB[,4:23]) ### Create matrix with time series of abundances
 
    ### Simone Code with clustering
    ctrl<-setOptions(clusters.crit = "BIC", clusters.qmin  = 2,clusters.qmax  = 10)
    res.clust = simone(abundance_mat,type="time-course",clustering=TRUE,control=ctrl)
    g.clust=getNetwork(res.clust,"BIC") 
    plot(g.clust,type="circles")
    
    ### Compute false positives and negatives
    Ahat =g.clust$Theta
    rates = ratesClassif(Ahat,causality_matrix)
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
    
    #### Pairwise GC code 
    lag_order = lagOrder(abundance_mat)
    pGC = pairwiseGC(abundance_mat,alphaLevel,lag_order) ## alpha-level global higher than 5%
    rates2 = ratesClassif(pGC,causality_matrix)
    resultsC2 = diagnosticsClassif(rates2)
    
    
    modelT= model
    nrepeat = ksite
    FPR = resultsC2[1]
    TPR = resultsC2[2]
    Precision = resultsC2[3]
    newscoresClassif2 = data.frame(FPR,TPR,Precision,nrepeat,modelT)
    if ((ksite == 1)&(model=="randomLV"))
    {
      scores_pGC= newscoresClassif2
    } else
    {
      scores_pGC = rbind(scores_pGC,newscoresClassif2)
    }
  }
}


pdf(file = paste("../figures/ROC_Simone_BIC_clustering_woutIntraSpInteractions_BHcorrection_alpha",alphaLevel,".pdf",sep=""),width=16,height = 8)
par(pty="s",mfrow=c(1,2),cex=1.5)
plot(scoresClassif$FPR[scoresClassif$modelT=="randomLV"],scoresClassif$TPR[scoresClassif$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC Simone")
abline(a=0,b=1,lwd=2)
lines(scoresClassif$FPR[scoresClassif$modelT=="randomVAR"],scoresClassif$TPR[scoresClassif$modelT=="randomVAR"],pch=19,type="p",col="yellow")
legend("right",legend=modelType,col=c("black","yellow"),pch=19,cex=0.8)

plot(scores_pGC$FPR[scores_pGC$modelT=="randomLV"],scores_pGC$TPR[scores_pGC$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)", main = "ROC pairwise Granger")
abline(a=0,b=1,lwd=2)
lines(scores_pGC$FPR[scores_pGC$modelT=="randomVAR"],scores_pGC$TPR[scores_pGC$modelT=="randomVAR"],pch=19,type="p",col="yellow")
legend("right",legend=modelType,col=c("black","yellow"),pch=19,cex=0.8)
dev.off()
### NB I can use various criteria -- like number of edges to make variations along the ROC curve
### In the context of pairwise GC, I can change the significance level to produce the ROC curve. 

write.csv(scoresClassif,file="../results/withoutIntraSp/clustering/scoresClassif_BHcorrection.csv")
write.csv(scores_pGC,file="../results/withoutIntraSp/clustering/scores_pGC_BHcorrection.csv")
