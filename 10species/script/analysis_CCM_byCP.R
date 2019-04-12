rm(list=ls())
graphics.off()

library(vars)
library(rEDM)

set.seed(666)

### Parameters
nsites = 25 ### how many samples of time series (sites or repeats)

### Useful variables
#modelType = c("refLV","refVAR","randomLV","randomVAR")
modelType = c("refLV","refVAR") #For now, we're just gonna focus on 10 species
nrepeat = 1:nsites
nmodels = length(modelType)### For model and parameter types

### In this analysis we use only one network geometry and dimensionality
### Quantitative parameters and initial conditions for the VAR and Lotka-Volterra models are however varied

### Path to files with the time series and other useful data
pathrefLV = "../data/ref_param_set/Data_wTime_abs_LV.csv" #10 species
pathrefVAR = "../data/ref_param_set/Data_wTime_abs_VAR.csv" #10 species
pathrandomLV = "../../20species/data/Data_wTime_abs_LV.csv" #20
pathrandomVAR = "../../20species/data/Data_wTime_abs_VAR.csv" #20



######################################## Utilitary functions #############################################

path_to_file <- function(model) {
  eval(as.name(paste("path",model,sep="")))
}### use  strsplit for more fancy stuff

ccm_test = function(z){

  ### CCM Analysis 
  species12=data.frame(1:nrow(z),z) #beware z is in log-scale
  names(species12)=c("time","sp1","sp2")
  libsizes = seq(10, nrow(z)-10, by = 10) #Before, we only got up to 80, which is not fair because we consider 700 points in the GC test
  lm=length(libsizes)
	numsamples = 100

#CP : FIrst, we need to choose the lag order
simplex_output_predictx = simplex(species12$sp1,E=1:10)
 lag_order_inter_CCM_predictx = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]
simplex_output_predicty = simplex(species12$sp2,E=1:10)
 lag_order_inter_CCM_predicty = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


 sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx, lib_column = "sp1",
                      target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=F)
  #can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
  sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty, lib_column = "sp2", target_column = "sp1",
                      lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=F)
  #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

  ### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
  rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
  rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==max(sp1_xmap_sp2$lib_size)]

  deltarho1xmap2= mean(rho1xmap2_Lmax_random- rho1xmap2_Lmin_random)

  # Fraction of samples for which rho(L_max)<rho(L_min)
  Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

  rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
  rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==max(sp2_xmap_sp1$lib_size)]

  deltarho2xmap1= mean(rho2xmap1_Lmax_random- rho2xmap1_Lmin_random)
  Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

  ### Let's try to see those
  sp1_xmap_sp2_means <- ccm_means(sp1_xmap_sp2)
  sp2_xmap_sp1_means <- ccm_means(sp2_xmap_sp1)

  RhoLMax_12=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==max(sp1_xmap_sp2$lib_size)] # 1 causes 2 if 2 xmap 1
  RhoLMax_21=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==max(sp2_xmap_sp1$lib_size)] # 2 causes 1 if 1 xmap 2

  ###Trying out my method, just to be sure it does not work

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp2"]=sample(species12[,"sp2"])
        sp1_xmap_sp2_random <- ccm(species_random, E = lag_order_inter_CCM_predictx, lib_column = "sp1",target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp1_xmap_sp2_random$rho
}
  Pval_1xmap2_bis = sum(rho_dist>RhoLMax_21)/numsamples

rho_dist=rep(NA,numsamples)
for (i in 1:numsamples){
        species_random=species12
        species_random[,"sp1"]=sample(species12[,"sp1"])
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_inter_CCM_predicty, lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1)
        rho_dist[i]=sp2_xmap_sp1_random$rho
}
  Pval_2xmap1_bis = sum(rho_dist>RhoLMax_12)/numsamples
	tmp_val=c(Pval_2xmap1,Pval_1xmap2,RhoLMax_12,RhoLMax_21,deltarho1xmap2,deltarho2xmap1,Pval_1xmap2_bis,Pval_2xmap1_bis,lag_order_inter_CCM_predictx,lag_order_inter_CCM_predicty)

  return(tmp_val) ### NB we may find a way to output rho as well in a meaningful manner

}

pairwiseCCM <-function(x){ ### returns a matrix of causal links based on pairwise CCM
  ### Adjusted p-values with Benjamini-Hochberg correction
  nspecies=ncol(x) ## entry matrix has time on rows and species in columns
  p_value=matrix(0,nrow=nspecies,ncol=nspecies)
  rho_lmax=matrix(0,nrow=nspecies,ncol=nspecies)
  delta_rho=matrix(0,nrow=nspecies,ncol=nspecies)
  lagorder1=rep(0,nspecies)
  lagorder2=rep(0,nspecies)
  p_value_bis=matrix(0,nrow=nspecies,ncol=nspecies)

  for (i in 1:nspecies){
    for (j in 1:nspecies){
      # cause first and effet later in grangertest()
      if (i >j){
        z=cbind(x[,i],x[,j])
        pccm=ccm_test(z)
        p_value[i,j] = pccm[1]
        p_value[j,i] = pccm[2]
        rho_lmax[i,j] = pccm[3]
        rho_lmax[j,i] = pccm[4]
        delta_rho[i,j] = pccm[5]
        delta_rho[j,i] = pccm[6]
        p_value_bis[i,j] = pccm[8] #sic
        p_value_bis[j,i] = pccm[7] #sic
	lagorder1[i]=pccm[9]
	lagorder2[j]=pccm[10]
      }
    }
  }
  lagorder1[1]=lagorder2[1]
  p_value_adj=p.adjust(p_value,method="BH")
  p_value_adj = matrix(p_value_adj,nrow=nspecies,ncol=nspecies)
  p_value_bis_adj=p.adjust(p_value_bis,method="BH")
  p_value_bis_adj = matrix(p_value_bis_adj,nrow=nspecies,ncol=nspecies)
  return(list(p_value,rho_lmax,delta_rho,p_value_adj,p_value_bis,p_value_bis_adj,lagorder1,lagorder2))
}


  for (model in modelType){
	print(model)	
    mat_tmp_rw=matrix(NA,nrow=100*nsites,ncol=11)
	colnames(mat_tmp_rw)=c("site","sp1","sp2","E1","E2","pvalCobeyBaskerville","rhomax","deltarho","pvalCobeyBaskerville_adj","pvalCP","pvalCP_adj")
	ijk=0
for (ksite in 1:25){ ### for sites or repeats
  print(ksite)

    ### Selects the files and then time series

    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    DB=DB[DB$Time_index %in% 201:500,] ## Select 300 last timesteps for 10 speices

    abundance_mat=as.matrix(DB[,4:13]) 

    resCCM = pairwiseCCM(abundance_mat) ## alpha-level global higher than 5%
    for(i in 1:10){
	for(j in 1:10){
		ijk=ijk+1
		mat_tmp_rw[ijk,1:5]=c(ksite,i,j,resCCM[[7]][i],resCCM[[8]][j])
		for(res in 1:6){
		mat_tmp_rw[ijk,5+res]=resCCM[[res]][i,j]
		}
	} #j in species
	} #i in species
} #ksite
	write.csv(mat_tmp_rw,paste('../results/10species_CCM_per_interaction_',model,".csv",sep=""))
} #model
