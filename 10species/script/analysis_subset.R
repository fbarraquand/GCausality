########################################################################################################
### Trying to compute causality on the most interacting species (for the 10-species model, it's the effect of sp2 on sp1)
### Based on FBarraquand's idea and his code in comparisonGranger_vs_CCM_2species/*
### CP April 2019
########################################################################################################

library(rEDM)
library(vars)

### Parameters
nsites = 25 ### how many samples of time series (sites or repeats)
nmodels = 4 ### For model and parameter types

### Useful variables
modelType = c("refLV","refVAR","randomLV","randomVAR")
nrepeat = 1:nsites

### Path to files with the time series and other useful data
pathrefLV = "../data/ref_param_set/Data_wTime_abs_LV.csv" #10 species
pathrefVAR = "../data/ref_param_set/Data_wTime_abs_VAR.csv" #10 species
pathrandomLV = "../../20species/data/Data_wTime_abs_LV.csv" #20
pathrandomVAR = "../../20species/data/Data_wTime_abs_VAR.csv" #20


#Initializing vectors 
ncond=nsites 
Pval_12_inter_GC=Pval_21_inter_GC=matrix(NA,ncond,2)
Pval_12_inter_CCM=Pval_21_inter_CCM=matrix(NA,ncond,2)
RhoLMax_12_inter=RhoLMax_21_inter=matrix(NA,ncond,2)
index_1cause2_inter_GC=index_2cause1_inter_GC=matrix(NA,ncond,2)
index_1cause2_inter_CCM=index_2cause1_inter_CCM=matrix(NA,ncond,2)
lag_order_inter_GC=matrix(NA,ncond,2)
lag_order_inter_CCM_predictx=lag_order_inter_CCM_predicty=matrix(NA,ncond,2)


######################################## Utilitary functions #############################################

path_to_file <- function(model) {
  eval(as.name(paste("path",model,sep="")))
}### use  strsplit for more fancy stuff

######################################

sp1=1
sp2=2

m=0
  for (model in c("randomLV","randomVAR")) #This is for 20 species, it would be refLV and refVAR for the 10 species
    {
	m=m+1
    ### Selects the files and then time series

for (ksite in 1:nsites){ ### for sites or repeats
print(ksite)
    DBall=read.csv(path_to_file(model))
    DB=DBall[DBall$Site==ksite,] ## Select a site
    head(DB)
    #DB=DB[DB$Time_index %in% 201:500,] ## Select 300 last timesteps for 10 species
    DB=DB[DB$Time_index %in% 301:1000,] ## Select 700 last timesteps for 20 species
    head(DB) ## 

    #abundance_mat=as.matrix(DB[,4:13]) ## for 10 species
    abundance_mat=as.matrix(DB[,4:23]) ## for 20 species
    abundance_mat=as.matrix(abundance_mat[,c(sp1,sp2)])
    time=1:nrow(abundance_mat)
    x=abundance_mat[,1]
    y=abundance_mat[,2]
	

#LOG ? 
#Already in the generation file

#CENTER ?   
	x=x-mean(x)
	y=y-mean(y)
	z=data.frame(time,x,y)
	names(z)=c("time","sp1","sp2")
#first, plot

	pdf(paste("../figures/subset_",model,"_ksite",ksite,"20species.pdf",sep=""))
	plot(time,x,ylim=c(min(z[,2:3]),max(z[,2:3])),t="l",col="black")
	lines(time,y,col="red")
	dev.off()

#Then, GC

varcompet<-VAR(y=data.frame(cbind(x,y)), type="none",lag.max=20,ic="SC")
lag_order_inter_GC[ksite,m] <- varcompet$p

gxy = grangertest(x,y,order = lag_order_inter_GC[ksite,m]) #x causes y 
gyx = grangertest(y,x,order = lag_order_inter_GC[ksite,m]) #y causes x

Pval_12_inter_GC[ksite,m]=gxy$`Pr(>F)`[2]
Pval_21_inter_GC[ksite,m]=gyx$`Pr(>F)`[2]

if (Pval_12_inter_GC[ksite,m]<0.1)
{index_1cause2_inter_GC[ksite,m]=1} else {index_1cause2_inter_GC[ksite,m]=0}

if (Pval_21_inter_GC[ksite,m]<0.1)
{index_2cause1_inter_GC[ksite,m]=1} else {index_2cause1_inter_GC[ksite,m]=0}


#### AND CCM

#Chose E
simplex_output_predictx = simplex(x,E=1:10)
 lag_order_inter_CCM_predictx[ksite,m] = simplex_output_predictx$E[which(simmplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_inter_CCM_predicty[ksite,m] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(time,z)
names(species12)=c("time","sp1","sp2")
libsizes = seq(10, nrow(z)-10, by = 10)
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[ksite,m] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs = TRUE,num_samples = numsamples,replace=FALSE)
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[ksite,m] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs = TRUE, num_samples = numsamples,replace=FALSE) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means <- data.frame(ccm_means(sp1_xmap_sp2), sd.rho = with(sp1_xmap_sp2,
                                                                            tapply(rho, lib_size, sd)))
sp2_xmap_sp1_means <- data.frame(ccm_means(sp2_xmap_sp1), sd.rho = with(sp2_xmap_sp1,
                                                                            tapply(rho, lib_size, sd)))


tab_simu[ksite,,1]=sp1_xmap_sp2_means$rho
tab_simu[ksite,,2]=sp2_xmap_sp1_means$rho

pdf(paste("../figures/rho_f_libsize_",model,"_site",ksite,"_20species.pdf"))
       plot(libsizes,tab_simu[ksite,,1],col=rgb(1,0,0,1),ylim=c(-1,1),t="l")
	 lines(libsizes,tab_simu[ksite,,2],col=rgb(0,0,1,1))
        lines(libsizes,sp1_xmap_sp2_means$rho+2*sp1_xmap_sp2_means$sd.rho,col=rgb(1,0,0,1),lty=2)
        lines(libsizes,sp1_xmap_sp2_means$rho-2*sp1_xmap_sp2_means$sd.rho,col=rgb(1,0,0,1),lty=2)
        lines(libsizes,sp2_xmap_sp1_means$rho+2*sp2_xmap_sp1_means$sd.rho,col=rgb(0,0,1,1),lty=2)
        lines(libsizes,sp2_xmap_sp1_means$rho-2*sp2_xmap_sp1_means$sd.rho,col=rgb(0,0,1,1),lty=2)
	abline(h=0.1,col="black",lty=2)
	abline(h=0.25,col="black",lty=2)
	legend("topleft",col=c("red","blue"),c("2 causes 1","1 causes 2"),bty="n",lty=1)
dev.off()


### Using the same method for producing P-values as Cobey and Baskerville PloS One 2016
rho1xmap2_Lmin_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[1]]
rho1xmap2_Lmax_random <- sp1_xmap_sp2$rho[sp1_xmap_sp2$lib_size ==libsizes[lm]]
# Fraction of samples for which rho(L_max)<rho(L_min)
Pval_1xmap2 = sum(rho1xmap2_Lmax_random<rho1xmap2_Lmin_random)/numsamples #2 towards 1

rho2xmap1_Lmin_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[1]]
rho2xmap1_Lmax_random <- sp2_xmap_sp1$rho[sp2_xmap_sp1$lib_size ==libsizes[lm]]

Pval_2xmap1 = sum(rho2xmap1_Lmax_random<rho2xmap1_Lmin_random)/numsamples #1 towards 2

Pval_12_inter_CCM[ksite,m]=Pval_2xmap1 # 1 causes 2 if 2 xmap 1
Pval_21_inter_CCM[ksite,m]=Pval_1xmap2 # 2 causes 1 if 1 xmap 2

RhoLMax_12_inter[ksite,m]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==libsizes[lm]] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[ksite,m]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==libsizes[lm]] # 2 causes 1 if 1 xmap 2

if ((Pval_12_inter_CCM[ksite,m]<0.1)&(RhoLMax_12_inter[ksite,m]>0.1))
{index_1cause2_inter_CCM[ksite,m]=1} else {index_1cause2_inter_CCM[ksite,m]=0}

if ((Pval_21_inter_CCM[ksite,m]<0.1)&(RhoLMax_21_inter[ksite,m]>0.1))
{index_2cause1_inter_CCM[ksite,m]=1} else {index_2cause1_inter_CCM[ksite,m]=0}


	}
DataCompet_subset_10species = data.frame(1:ncond,lag_order_inter_GC[,m],Pval_12_inter_GC[,m],Pval_21_inter_GC[,m],index_1cause2_inter_GC[,m],index_2cause1_inter_GC[,m],lag_order_inter_CCM_predictx[,m],Pval_12_inter_CCM[,m],lag_order_inter_CCM_predicty[,m],Pval_21_inter_CCM[,m],index_1cause2_inter_CCM[,m],index_2cause1_inter_CCM[,m])
names(DataCompet_subset_10species)=c('site','lag_order_inter_GC','Pval_12_inter_GC','Pval_21_inter_GC','index_1cause2_inter_GC','index_2cause1_inter_GC','lag_order_inter_CCM_predictx','Pval_12_inter_CCM','lag_order_inter_CCM_predicty','Pval_21_inter_CCM','index_1cause2_inter_CCM','index_2cause1_inter_CCM')

write.csv(DataCompet_subset_10species,file=paste("../results/DataCompet_subset",model,"_20species.csv",sep=""))


}
