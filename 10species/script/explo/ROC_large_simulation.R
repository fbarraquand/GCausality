########################################################################################################
### From FBarraquand, analysis_Simone_Clustering_woutIntraSp.R 
###CP April 2019
########################################################################################################
#CP 09/04/2019, this is basically a copy-paste of two FB's codes to have only one figure for the 10- and 20-species models. This only presents results for the Granger-causality and might soon be deprecated.
#CP 16/04/19 : Removed alpha hard-coding


graphics.off()
rm(list=ls())

### Loading libraries
library(simone)

### Parameters
nsites = 25 ### how many samples of time series (sites or repeats)
nmodels = 4 ### For model and parameter types
alpha=0.2

### Useful variables
modelType = c("refLV","refVAR")
#modelType = c("refLV","refVAR","randomLV","randomVAR")
nrepeat = 1:nsites

### In this analysis we use only one network geometry and dimensionality
### Quantitative parameters and initial conditions for the VAR and Lotka-Volterra models are however varied

### Path to files with the time series and other useful data
pathrefLV = "./data/ref_param_set/Data_wTime_abs_LV.csv"
pathrefVAR = "./data/ref_param_set/Data_wTime_abs_VAR.csv"


######################################## Utilitary functions #############################################
path_to_file <- function(model) {
  eval(as.name(paste("path",model,sep="")))
}### use  strsplit for more fancy stuff


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

library(vars)

lagOrder <-function(y,lag.max=15){ #lag-order as determined by full VAR model
  yd=data.frame(abundance_mat)
  varresult = VAR(y=yd, type="none",lag.max=lag.max)
  return(as.numeric(varresult$p))
}

pairwiseGC <-function(x,alphaLevel,lagorder){ ### returns a matrix of causal links based on pairwise GC
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
  pGC=(p_value<alphaLevel)
  return(pGC)
}
################################ end of utilitary functions ##################################

#adjacency matrix for comparison to estimated matrices
causality_matrix = rbind(c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,1,0,0,0,0,0,0),
                         c(0,0,0,1,1,1,1,0,0,0),
                         c(0,0,0,0,1,1,1,0,0,0),
                         c(0,0,0,0,1,1,1,0,0,0),
                         c(0,0,0,0,0,0,0,1,1,1),
                         c(0,0,0,0,0,0,0,1,1,1),
                         c(0,0,0,0,0,0,0,1,1,1))


for (ksite in 1:nsites){ ### for sites or repeats

  for (model in c("refLV","refVAR"))
    {

    ### Selects the files and then time series

    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    DB=DB[DB$Time_index %in% 201:500,] ## Select 300 last timesteps

    abundance_mat=as.matrix(DB[,4:13]) ### Create matrix with time series of abundances

    #### Simone code - clustering
    ctrl<-setOptions(clusters.crit = "BIC")
    res.clust = simone(abundance_mat,type="time-course",clustering=TRUE,control=ctrl)
    g.clust=getNetwork(res.clust,"BIC")

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
    if ((ksite == 1)&(model=="refLV"))
    {
      scoresClassif = newscoresClassif
    } else
    {
      scoresClassif = rbind(scoresClassif,newscoresClassif)
    }

    #### Pairwise GC code 
    lag_order = lagOrder(abundance_mat)
    pGC = pairwiseGC(abundance_mat,alpha,lag_order) ## alpha-level 0.05
    rates2 = ratesClassif(pGC,causality_matrix)
    resultsC2 = diagnosticsClassif(rates2)


    modelT= model
    nrepeat = ksite
    FPR = resultsC2[1]
    TPR = resultsC2[2]
    Precision = resultsC2[3]
    newscoresClassif2 = data.frame(FPR,TPR,Precision,nrepeat,modelT)
    if ((ksite == 1)&(model=="refLV"))
    {
      scores_pGC= newscoresClassif2
    } else
    {
      scores_pGC = rbind(scores_pGC,newscoresClassif2)
    }
  }
}


pdf(file = "ROC_Simone_BIC_ClusteringMixer_10_and_20_species.pdf",width=16,height = 16)
par(pty="s",mfrow=c(2,2),cex=1.5,mar=c(0.5,4,4,0.5))
plot(scoresClassif$FPR[scoresClassif$modelT=="refLV"],scoresClassif$TPR[scoresClassif$modelT=="refLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "",ylab ="True Positive Rate (recall)", main = "ROC Simone",xaxt="n")
abline(a=0,b=1,lwd=2)
lines(scoresClassif$FPR[scoresClassif$modelT=="refVAR"],scoresClassif$TPR[scoresClassif$modelT=="refVAR"],pch=19,type="p",col="blue")
legend("right",legend=modelType,col=c("black","blue"),pch=19,cex=0.8,bty="n")
mtext("a)",line=2.5,side=2,at=1,las=2,cex=1.5)


plot(scores_pGC$FPR[scores_pGC$modelT=="refLV"],scores_pGC$TPR[scores_pGC$modelT=="refLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "",ylab ="", main = "ROC pairwise Granger",xaxt="n")
abline(a=0,b=1,lwd=2)
lines(scores_pGC$FPR[scores_pGC$modelT=="refVAR"],scores_pGC$TPR[scores_pGC$modelT=="refVAR"],pch=19,type="p",col="blue")
mtext("b)",line=2.5,side=2,at=1,las=2,cex=1.5)

############## And now, 20 species
modelType = c("randomLV","randomVAR")

### Path to files with the time series and other useful data
pathrandomLV = "../20species/data/Data_wTime_abs_LV.csv"
pathrandomVAR = "../20species/data/Data_wTime_abs_VAR.csv"



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

for (ksite in 1:nsites){ ### for sites or repeats

  for (model in c("randomLV","randomVAR"))
    {

    ### Selects the files and then time series

    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps

    abundance_mat=as.matrix(DB[,4:23]) ### Create matrix with time series of abundances

    ### Simone Code with clustering
    ctrl<-setOptions(clusters.crit = "BIC", clusters.qmin  = 2,clusters.qmax  = 10)
    res.clust = simone(abundance_mat,type="time-course",clustering=TRUE,control=ctrl)
    g.clust=getNetwork(res.clust,"BIC")

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
    pGC = pairwiseGC(abundance_mat,alpha,lag_order) ## alpha-level global higher than 5%
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

par(mar=c(4,4,0.5,0.5))
plot(scoresClassif$FPR[scoresClassif$modelT=="randomLV"],scoresClassif$TPR[scoresClassif$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)")
abline(a=0,b=1,lwd=2)
lines(scoresClassif$FPR[scoresClassif$modelT=="randomVAR"],scoresClassif$TPR[scoresClassif$modelT=="randomVAR"],pch=19,type="p",col="blue")
legend("right",legend=modelType,col=c("black","blue"),pch=19,cex=0.8,bty="n")
mtext("c)",las=2,side=2,line=2.5,at=1,cex=1.5)

plot(scores_pGC$FPR[scores_pGC$modelT=="randomLV"],scores_pGC$TPR[scores_pGC$modelT=="randomLV"],pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="")
abline(a=0,b=1,lwd=2)
lines(scores_pGC$FPR[scores_pGC$modelT=="randomVAR"],scores_pGC$TPR[scores_pGC$modelT=="randomVAR"],pch=19,type="p",col="blue")
mtext("d)",las=2,side=2,line=2.5,at=1,cex=1.5)
dev.off()
