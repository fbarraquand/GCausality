##############################################################
### FB 15/06/2018 -- Analysis for reference parameter set ####
### Uses Simone and pairwise GC [VAR(p)] inference        ####
### Model selection is done by BIC                        ####
##############################################################

rm(list = ls())
set.seed(5)
library(simone)

nsites = 25 ### how many samples of time series (sites or repeats)

### Useful variables
TPR=FPR=Recall=Precision=rep(0,nsites)
TPR2=FPR2=Recall2=Precision2=rep(0,nsites)

###

#adjacency matrix for comparison
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

#### Loads the data

### LV data
DBall=read.csv("../data/ref_param_set/Data_wTime_abs_LV.csv") 
### data/random_param_set/ next

for (ksite in 1:nsites){ ### for sites or repeats
  
  DB=DBall[DBall$Site==ksite,] ## Select a site
  head(DB)
  DB=DB[DB$Time_index %in% 201:500,]
  head(DB) ## 
  
  abundance_mat=as.matrix(DB[,4:13]) ### Create matrix with time series of abundances

  #### Simone code
  ### First no clusters
  res.noclust<-simone(abundance_mat,type="time-course") #OK
  g.no<-getNetwork(res.noclust,selection = "BIC") 
  plot(g.no,type="cluster")
  plot(g.no,type="circles")
  
  ### Compute false positives and negatives
  Pos = (abs(g.no$Theta)>0)
  #Neg = (abs(g.mod$Theta)==0)
  TruePos =  (abs(causality_matrix)>0)
  #TrueNeg =  (abs(B)==0)
  
  TP = sum((TruePos == TRUE) & (Pos == TRUE))
  FN = sum((TruePos == TRUE) & (Pos == FALSE))
  FP = sum((TruePos == FALSE) & (Pos == TRUE))
  TN = sum((TruePos == FALSE) & (Pos == FALSE))
  
  Precision[ksite] = TP/(TP+FP)  #True positives over tested positives (percentage positives correct)
  Recall[ksite] = TP/(TP+FN)     #True positives over relevant elements (power to find the true positives)
  ### ROC = c(FPR,TPR)
  
  TPR[ksite] = Recall[ksite]
  FPR[ksite] = FP / (FP + TN)  # = fall-out rate  = 1 - specificity = 1 - true negative rate
  #### Pairwise GC code 

}
plot(FPR,TPR,pch=19)
plot(FPR2,TPR2,pch=19,col="blue")

### VAR-simulated data next
DBall=read.csv("../data/ref_param_set/Data_wTime_abs_VAR.csv") 

for (ksite in 1:nsites){ ### for sites or repeats
  
  DB=DBall[DBall$Site==ksite,] ## Select a site
  head(DB)
  DB=DB[DB$Time_index %in% 201:500,]
  head(DB) ## 
  
  abundance_mat=as.matrix(DB[,4:13]) ### Create matrix with time series of abundances
  
  #### Simone code
  ### First no clusters
  res.noclust<-simone(abundance_mat,type="time-course") #OK
  g.no<-getNetwork(res.noclust,selection = "BIC") 
  plot(g.no,type="cluster")
  plot(g.no,type="circles")
  
  ### Compute false positives and negatives
  Pos = (abs(g.no$Theta)>0)
  #Neg = (abs(g.mod$Theta)==0)
  TruePos =  (abs(causality_matrix)>0)
  #TrueNeg =  (abs(B)==0)
  
  TP = sum((TruePos == TRUE) & (Pos == TRUE))
  FN = sum((TruePos == TRUE) & (Pos == FALSE))
  FP = sum((TruePos == FALSE) & (Pos == TRUE))
  TN = sum((TruePos == FALSE) & (Pos == FALSE))
  
  Precision2[ksite] = TP/(TP+FP)  #True positives over tested positives (percentage positives correct)
  Recall2[ksite] = TP/(TP+FN)     #True positives over relevant elements (power to find the true positives)
  ### ROC = c(FPR,TPR)
  
  TPR2[ksite] = Recall2[ksite]
  FPR2[ksite] = FP / (FP + TN)  # = fall-out rate  = 1 - specificity = 1 - true negative rate
  #### Pairwise GC code 
  
}

par(pty="s")
plot(FPR,TPR,pch=19,xlim=c(0,1),ylim=c(0,1),xlab = "False Positive Rate (1 - specificity)",ylab ="True Positive Rate (recall)")
abline(a=0,b=1,lwd=2)
lines(FPR2,TPR2,pch=19,type="p",col="blue")

#### Devise function to compute recall, precision, and ROC curves

### NB I can use various criteria -- like number of edges to make variations along the ROC curve

### In the context of pairwise GC, I can change the significance level to produce the ROC curve. 
