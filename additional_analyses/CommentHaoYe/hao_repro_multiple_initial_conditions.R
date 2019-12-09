library(dplyr)
library(purrr)
library(ggplot2)
library(rEDM)

set.seed(42)

# functions
simulate_data <- function(burn_in = 500, n = 300, 
                          init_x = runif(1,0.1,0.7), 
                          init_y = runif(1,0.1,0.7))
{
    obs_idx <- burn_in + seq(n)
    
    x <- numeric(burn_in + n)
    y <- numeric(burn_in + n)
    x[1] <- init_x
    y[1] <- init_y
    for (t in 1:(burn_in + n - 1))
    {
        x[t + 1] <- x[t] * (3.8 - 3.8 * x[t])
        y[t + 1] <- y[t] * (3.5 - 3.5 * y[t])
    }
    
    x <- log(x[obs_idx])
    y <- log(y[obs_idx])
    x <- x - mean(x)
    y <- y - mean(y)
    data.frame(x, y)
}

do_ccm <- function(dat, lib_sizes = 100, ...)
{
    bind_rows(ccm(dat, E = 1, lib_column = 1, target_column = 2, 
                  replace = FALSE, silent = TRUE, 
                  lib_sizes = lib_sizes, ...) %>% 
                  ccm_means() %>%
                  mutate(dir = "x -> y"), 
              ccm(dat, E = 1, lib_column = 2, target_column = 1, 
                  replace = FALSE, silent = TRUE, 
                  lib_sizes = lib_sizes, ...) %>%
                  ccm_means() %>%
                  mutate(dir = "y -> x"))
}

init <- readRDS("init_values.RDS")
dat <- simulate_data(init_x = init$x, init_y = init$y)
out <- do_ccm(dat)

# random shuffle surrogates
x_surr <- make_surrogate_data(dat$x)
out_surr <- map_dfr(seq_len(NCOL(x_surr)), function(j) {
    temp_dat <- data.frame(x = x_surr[, j], y = dat$y)
    do_ccm(temp_dat)
})

# random shuffle, maintaining periodicity
make_periodic_surrogate <- function(ts, num_surr = 100, period = 2)
{
    n <- length(ts)
    matrix(unlist(
        lapply(seq(num_surr), function(i) {
            surr_ts <- numeric(n)
            for (phase in seq_len(period))
            {
                idx <- seq(from = phase, to = n, by = period)
                surr_ts[idx] <- ts[sample(idx)]
            }
            surr_ts
        })
    ), ncol = num_surr)
}
x_surr_periodic <- make_periodic_surrogate(dat$x, period = 8)
out_surr_periodic <- map_dfr(seq_len(NCOL(x_surr_periodic)), function(j) {
    temp_dat <- data.frame(x = x_surr_periodic[, j], y = dat$y)
    do_ccm(temp_dat)
})








#### Addition by CP to reproduce our code, that is: multiple initial conditions and computation of p-values based on the surrogates. We're going to use the completely random surrogates and the periodicity-driven surrogates
ncond=100
tmax=800
numsamples=100

X=rep(NA,ncond)
Y=rep(NA,ncond)
simplex_output_predictx=rep(NA,ncond)
simplex_output_predicty=rep(NA,ncond)
lag_order_inter_CCM_predictx=rep(NA,ncond)
lag_order_inter_CCM_predicty=rep(NA,ncond)
  Pval_21_inter_CCM_random=  Pval_12_inter_CCM_random=  Pval_12_inter_CCM_periodic=  Pval_21_inter_CCM_periodic=rep(NA,ncond)
RhoLMax_12_inter=RhoLMax_21_inter=rep(NA,ncond)
Pval_2xmap1=Pval_1xmap2=rep(NA,ncond)

for (kcond in 1:ncond){
print(kcond)
X[kcond]=runif(1,0.1,0.7)
Y[kcond]=runif(1,0.1,0.7)
dat <- simulate_data(init_x = X[kcond], init_y =Y[kcond])

x=dat$x
y=dat$y
z=dat

simplex_output_predictx = simplex(x,E=1:10)
 lag_order_inter_CCM_predictx[kcond] = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(y,E=1:10)
 lag_order_inter_CCM_predicty[kcond] = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

 ### CCM Analysis 
species12=data.frame(501:800,z)
names(species12)=c("time","sp1","sp2")
libsizes = c(5,8,seq(10, nrow(species12)-11, by = 10))
lm=length(libsizes)
numsamples = 100
sp1_xmap_sp2 <- ccm(species12, E = lag_order_inter_CCM_predictx[kcond] , lib_column = "sp1",
                    target_column = "sp2", lib_sizes = libsizes, random_libs=T,num_samples = numsamples,replace=T)  #We want the real bootstrap from CB
#can we reconstruct 2 from 1, i.e., does 2 CCM-cause 1?
sp2_xmap_sp1 <- ccm(species12, E = lag_order_inter_CCM_predicty[kcond] , lib_column = "sp2", target_column = "sp1",
                    lib_sizes = libsizes, random_libs =T, num_samples = numsamples,replace=T) #can we reconstruct 1 from 2, i.e., does 1 CCM-cause 2?

sp1_xmap_sp2_means=ccm_means(sp1_xmap_sp2)
sp2_xmap_sp1_means=ccm_means(sp2_xmap_sp1)

################ Surrogate

RhoLMax_12_inter[kcond]=sp2_xmap_sp1_means$rho[sp2_xmap_sp1_means$lib_size==max(sp2_xmap_sp1_means$lib_size)] # 1 causes 2 if 2 xmap 1
RhoLMax_21_inter[kcond]=sp1_xmap_sp2_means$rho[sp1_xmap_sp2_means$lib_size==max(sp1_xmap_sp2_means$lib_size)] # 2 causes 1 if 1 xmap 2

rho_dist_random=rep(NA,numsamples)
rho_dist_periodic=rep(NA,numsamples)

species_random=species12
species_random_2=make_surrogate_data(species12[,"sp2"],num_surr=numsamples)
species_random_1=make_surrogate_data(species12[,"sp1"],num_surr=numsamples)

species_periodic_2=make_periodic_surrogate(species12[,"sp2"],num_surr=numsamples)
species_periodic_1=make_periodic_surrogate(species12[,"sp1"],num_surr=numsamples)

for (i in 1:numsamples){
	species_random=cbind(species12[,"sp1"],species_random_2[,i])
	colnames(species_random)=c("sp1",'sp2')
        sp1_xmap_sp2_random <- ccm(species_random, E = lag_order_inter_CCM_predictx[kcond], lib_column = "sp1",target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), replace=FALSE,num_samples = 1,random_libs=F)
        rho_dist_random[i]=sp1_xmap_sp2_random$rho
	
	species_periodic=cbind(species12[,"sp1"],species_periodic_2[,i])
	colnames(species_periodic)=c("sp1",'sp2')
        sp1_xmap_sp2_periodic <- ccm(species_periodic, E = lag_order_inter_CCM_predictx[kcond], lib_column = "sp1",target_column = "sp2", lib_sizes = max(sp1_xmap_sp2$lib_size), replace=FALSE,num_samples = 1,random_libs=F)
        rho_dist_periodic[i]=sp1_xmap_sp2_periodic$rho
}
  Pval_21_inter_CCM_random[kcond] = (1+sum(rho_dist_random>RhoLMax_21_inter[kcond]))/(1+numsamples)
  Pval_21_inter_CCM_periodic[kcond] = (1+sum(rho_dist_periodic>RhoLMax_21_inter[kcond]))/(numsamples+1)

rho_dist_random=rep(NA,numsamples)
rho_dist_periodic=rep(NA,numsamples)
for (i in 1:numsamples){
	species_random=cbind(species_random_1[,i],species12[,"sp2"])
	colnames(species_random)=c("sp1",'sp2')
        sp2_xmap_sp1_random <- ccm(species_random, E = lag_order_inter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1,random_libs=F)
        rho_dist_random[i]=sp2_xmap_sp1_random$rho
        
	species_periodic=cbind(species_periodic_1[,i],species12[,"sp2"])
	colnames(species_periodic)=c("sp1",'sp2')
        sp2_xmap_sp1_periodic <- ccm(species_periodic, E = lag_order_inter_CCM_predicty[kcond], lib_column = "sp2",target_column = "sp1", lib_sizes = max(sp2_xmap_sp1$lib_size), replace=FALSE,num_samples = 1,random_libs=F)
        rho_dist_periodic[i]=sp2_xmap_sp1_periodic$rho
}
  Pval_12_inter_CCM_random[kcond] = (1+sum(rho_dist_random>RhoLMax_12_inter[kcond]))/(1+numsamples)
  Pval_12_inter_CCM_periodic[kcond] = (1+sum(rho_dist_periodic>RhoLMax_12_inter[kcond]))/(1+numsamples)

}

DataCompet_chaos = data.frame(1:ncond,lag_order_inter_CCM_predictx,lag_order_inter_CCM_predicty,Pval_12_inter_CCM_random,Pval_21_inter_CCM_random,Pval_12_inter_CCM_periodic,Pval_21_inter_CCM_periodic)
write.csv(DataCompet_chaos,file="simu_HaoYe_chaos_withoutinter.csv")


