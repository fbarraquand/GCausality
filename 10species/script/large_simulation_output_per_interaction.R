rm(list=ls())
graphics.off()

library(simone)
library(vars)

nsite=25

### Useful variables
modelType = c("refLV","refVAR","randomLV","randomVAR")
nrepeat = 1:nsite

### In this analysis we use only one network geometry and dimensionality
### Quantitative parameters and initial conditions for the VAR and Lotka-Volterra models are however varied

### Path to files with the time series and other useful data
pathrefLV = "./data/ref_param_set/Data_wTime_abs_LV.csv"
pathrefVAR = "./data/ref_param_set/Data_wTime_abs_VAR.csv"
pathrandomLV = "../20species/data/Data_wTime_abs_LV.csv"
pathrandomVAR = "../20species/data/Data_wTime_abs_VAR.csv"


######################################## Utilitary functions #############################################
path_to_file <- function(model) {
  eval(as.name(paste("path",model,sep="")))
}### use  strsplit for more fancy stuff


lagOrder <-function(y,lag.max=15){ #lag-order as determined by full VAR model
  yd=data.frame(abundance_mat)
  varresult = VAR(y=yd, type="none",lag.max=lag.max)
  return(as.numeric(varresult$p))
}

pairwiseGC <-function(x,lagorder){ ### returns a matrix of causal links based on pairwise GC
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
  pGC=p_value
  return(pGC) ###Warning, I have changed the function in order to have a p-value instead of a boolean
}

###################################################################

for (m in 1:length(modelType)){
	model=modelType[m]
    DBall=read.csv(path_to_file(model))
    nspecies=ncol(DBall)-3
    mat_tmp_rw=matrix(NA,nsite*nspecies*nspecies,6)
	colnames(mat_tmp_rw)=c("site","sp1","sp2","mat_simone","p_pairwise","p_pairwise_adj")
	ijk=0
for (ksite in 1:nsite){ ### for sites or repeats
    ### Selects the files and then time series

    DB=DBall[DBall$Site==ksite,] ## Select a site
	if(m<3){
    DB=DB[DB$Time_index %in% 201:500,] ## Select 300 last timesteps
    abundance_mat=as.matrix(DB[,4:13]) ### Create matrix with time series of abundances
	}else{
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps
    abundance_mat=as.matrix(DB[,4:23]) ### Create matrix with time series of abundances
	}

    #### Simone code - clustering
    ctrl<-setOptions(clusters.crit = "BIC")
    res.clust = simone(abundance_mat,type="time-course",clustering=TRUE,control=ctrl)
    g.clust=getNetwork(res.clust,"BIC")

    ### Compute false positives and negatives
    Ahat =g.clust$Theta

    #### Pairwise GC code 
    lag_order = lagOrder(abundance_mat)
    pGC = pairwiseGC(abundance_mat,lag_order) 
    p_value_adj=p.adjust(pGC,method="BH")
    pGC_adj = matrix(p_value_adj,nrow=nspecies,ncol=nspecies)


	for(i in 1:nspecies){
		for (j in 1:nspecies){
			ijk=ijk+1
			mat_tmp_rw[ijk,1]=ksite
			mat_tmp_rw[ijk,2]=i
			mat_tmp_rw[ijk,3]=j
			mat_tmp_rw[ijk,4]=Ahat[i,j]
			mat_tmp_rw[ijk,5]=pGC[i,j]
			mat_tmp_rw[ijk,6]=pGC_adj[i,j]
		}

	}
}
	write.csv(mat_tmp_rw,paste('./results/large_GC_per_interaction_',model,".csv",sep=""))
}
