### FB 08/10/2017 -- adapted from Nantes et Irstea presentations though
### Uses a strongly seasonal driver. A logical case for conditional GC. 

### Modified 18/06/2017 -- added code to test for significance using surrogates
### Install new rEDM package, last development version
# library(devtools)
# devtools::install_github("ha0ye/rEDM")
# see https://cran.r-project.org/web/packages/rEDM/vignettes/rEDM-tutorial.html

rm(list=ls())
.libPaths("/home/frederic/myRpackages/")
dev.off()

## Loading packages
library(rEDM)
library(vars)

### First case with interactions
set.seed(42)

tmax=300
sigma = 0.3 #process noise = sigma^2 = 0.1
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))
seasonality<-2*sin(2*pi*(1:tmax)/24)	# must be enough to affect the growth rates
### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-2*Y[t,2] + rnorm(1,0,sigma))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0.31*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,sigma))
}
y=log(Y)

varcompet2<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
causality(varcompet2,cause="y1") #p-value 0.004054
causality(varcompet2,cause="y2") #0.00000
### I may need to repeat this over many simulations

############
pdf(file="2competitors_wSharedSeasonalDriver.pdf",width=10,height=3)
par(cex=1.5,lwd=2,mar=c(4,4,1,2))
plot(100:200,y1[100:200],ylab="log(abundance)",col="green",type="l",xlab="Time",ylim=c(-10,max(y1)))
lines(100:200,y[100:200,1],col="blue")
lines(100:200,y[100:200,2],col="black")
dev.off()

### Second case without interactions
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))


for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-0*Y[t,2] + rnorm(1,0,sigma))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,sigma))
}
y=log(Y)

pdf(file="2competitors_wSharedSeasonalDriver_noInteractions.pdf",width=10,height=3)
par(cex=1.5,lwd=2,mar=c(4,4,1,2))
plot(100:200,y1[100:200],ylab="log(abundance)",col="green",type="l",xlab="Time",ylim=c(-10,max(y1)))
lines(100:200,y[100:200,1],col="blue")
lines(100:200,y[100:200,2],col="black")
dev.off()

varcompet2bis<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
summary(varcompet2bis)
causality(varcompet2bis,cause="y1") #p-value 1.073e-06
causality(varcompet2bis,cause="y2") #0.1953


############# Pairwise GC (without conditioning on the driver) & CCM to test absolutely first... 

varcompet_noExo=VAR(y, type="none",lag.max=10,ic="SC") ## Direct pairwise GC
summary(varcompet_noExo)
causality(varcompet_noExo,cause="y1") #p-value 0.03633
causality(varcompet_noExo,cause="y2") #0.00000

### We find interactions that are not present
 
### check with grangertest() -- Wald test for nested models
grangertest(y[,1],y[,2])
grangertest(y[,2],y[,1]) ### there we find a causality. 

########## CCM analysis of this dataset ######################

### Merging to a dataframe

S12T=data.frame(1:tmax,y,y1) #species2_species1_temp
names(S12T)=c("time","species1","species2","temp")
par(mfrow=c(3,1))
plot(S12T$time,S12T$species1,type='o')
plot(S12T$time,S12T$species2,type='o')
plot(S12T$time,S12T$temp,type='o')
spectrum(S12T$temp,method="ar") #1/24 = frequency

### Old analysis with previous package below, commented out
### New analysis
### Adapted from https://cran.r-project.org/web/packages/rEDM/vignettes/rEDM-tutorial.html

### Basic CCM analysis
varnames <- colnames(S12T[2:4])
n <- NROW(S12T)
ccm_matrix <- array(NA, dim = c(length(varnames), length(varnames)), dimnames = list(varnames, 
                                                                             varnames))

for (ccm_from in varnames) {
  for (ccm_to in varnames[varnames != ccm_from]) {
    out_temp <- ccm(S12T, E = 8, lib_column = ccm_from, target_column = ccm_to, 
                    lib_sizes = n, replace = FALSE, silent = TRUE)
    ccm_matrix[ccm_from, ccm_to] <- out_temp$rho
  }
}

### ------------- Generate surrogate time series ----------------------####
num_surr <- 1000 #number of surrogate time series 
surr_species1 <- make_surrogate_data(S12T$species1, method = "seasonal", 
                                 T_period = 24, num_surr = num_surr)
surr_species2 <- make_surrogate_data(S12T$species2, method = "seasonal", 
                                     T_period = 24, num_surr = num_surr)
surr_temp <- make_surrogate_data(S12T$temp, method = "seasonal", 
                                 T_period = 24, num_surr = num_surr) #previously T_period = 12, bug

### ------------ Compute the rho values for the surrogates -------------####
### Regarding the effect of temp and species 2 on species 1 
rho_surr1 <- data.frame(temp = numeric(num_surr), species2 = numeric(num_surr))
### Regarding the effect of temp and species 1 on species 2
rho_surr2 <- data.frame(temp = numeric(num_surr), species1 = numeric(num_surr))

for (i in 1:num_surr) {
  rho_surr1$temp[i] <- ccm(cbind(S12T$species1, surr_temp[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                           replace = FALSE)$rho
  
  rho_surr1$species2[i] <- ccm(cbind(S12T$species1, surr_species2[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                          replace = FALSE)$rho
}

### P-value for rejecting the hypothesis that seasonality both causes temp and species 1
(sum(ccm_matrix['species1','temp'] < rho_surr1$temp) + 1) / 
(length(rho_surr1$temp) + 1) ### we cannot reject the null hypothesis, unclear that temp -> species 1

# ---- Rationale ------
# This is Pr(rho(species xmap real temp)< rho(species xmap surrogate temp))
# if there is a real effect of temp, rho should be higher with real temp
# so ccm_matrix['species1','temp'] < rho_surr1$temp should not be frequent
# and the P-value should be very low for temp -> species 1

### P-value for rejecting the hypothesis that seasonality both causes species 1 and species 2
# (more interesting here)
(sum(ccm_matrix['species1','species2'] < rho_surr1$species2) + 1) / 
  (length(rho_surr1$species2) + 1) #0.12 ## weakly significant  

# here we look at sp1 xmap sp2

### Perhaps we can have a more refined surrogate? to condition directly on temp?
?make_surrogate_data
signal = rnorm(100) + sin(1:100 * pi / 6)
twin_surr = make_surrogate_twin(signal, 10)
par(mfrow=c(2,1))
plot(signal,type='o',pch=19)
matplot(twin_surr)
matlines(twin_surr)

### That's what we want -- let's redo the analysis for 1 -> 2 
surr_temp_twin <- make_surrogate_data(S12T$temp, method = "twin", num_surr = num_surr)
surr_species2_twin <- make_surrogate_data(S12T$species2, method = "twin", num_surr = num_surr)
rho_surr1_twin <- data.frame(temp = numeric(num_surr), species2 = numeric(num_surr))

for (i in 1:num_surr) {
  rho_surr1_twin$temp[i] <- ccm(cbind(S12T$species1, surr_temp_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                           replace = FALSE)$rho
  
  rho_surr1_twin$species2[i] <- ccm(cbind(S12T$species1, surr_species2_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                               replace = FALSE)$rho
}

### P-value for rejecting the hypothesis that seasonality both causes temp and species 1
(sum(ccm_matrix['species1','temp'] < rho_surr1_twin$temp) + 1) / 
  (length(rho_surr1$temp) + 1) # 0.87 
### conclusion: cannot say that species1 xmap temp (equivalent to temp -> species 1)

### P-value for rejecting the hypothesis that seasonality both causes species 1 and species 2
# (more interesting here)
(sum(ccm_matrix['species1','species2'] < rho_surr1_twin$species2) + 1) / 
  (length(rho_surr1$species2) + 1) ## 0.24 
### conclusion: cannot say that species 1 xmap species 2, i.e. species 2 causes species 1

### (NB we have to check the other way around, i.e., if there is an interaction the technique can detect it. 
### cf. Cobey and Baskerville results.)

### Other causal direction -- we try to have species2 xmap species 1 to see if 1 causes 2
surr_species1_twin <- make_surrogate_data(S12T$species1, method = "twin", num_surr = num_surr)
rho_surr2_twin <- data.frame(temp = numeric(num_surr), species1 = numeric(num_surr))

for (i in 1:num_surr) {
  rho_surr2_twin$temp[i] <- ccm(cbind(S12T$species2, surr_temp_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                replace = FALSE)$rho
  
  rho_surr2_twin$species2[i] <- ccm(cbind(S12T$species2, surr_species1_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                    replace = FALSE)$rho
}

### P-value for rejecting the hypothesis that temp does not cause species 2 fluctuations
(sum(ccm_matrix['species2','temp'] < rho_surr2_twin$temp) + 1) / 
  (length(rho_surr2$temp) + 1) ### 0.000999001
### we reject the hypothesis of joint fluctuation, species 2 xmap temp

### P-value for rejecting the hypothesis that species 1 does not cause species 2 fluctuations
(sum(ccm_matrix['species2','species1'] < rho_surr2_twin$species2) + 1) / 
  (length(rho_surr1$species2) + 1) ## 0.07692308 almost significant here 
### we find that species2 xmap species1 so that 1->2. 

######## Some more ideas ################################################
### What if we used surrogate techniques in a time-domain Granger context? 
### As we did for Spectral GC!! Would that produce a fair comparison? 


############################# Using many simulations to have distributions of P-values #############

# Use of a function to test for CCM between time series x and y given temp
ccm_surrogate_testing = function(y,y1,num_surr=1000){
  
  S12T=data.frame(1:tmax,y,y1) #species2_species1_temp
  names(S12T)=c("time","species1","species2","temp")
  ### Adapted from https://cran.r-project.org/web/packages/rEDM/vignettes/rEDM-tutorial.html
  
  ### Basic CCM analysis
  varnames <- colnames(S12T[2:4])
  n <- NROW(S12T)
  ccm_matrix <- array(NA, dim = c(length(varnames), length(varnames)), dimnames = list(varnames, 
                                                                                       varnames))
  
  for (ccm_from in varnames) {
    for (ccm_to in varnames[varnames != ccm_from]) {
      out_temp <- ccm(S12T, E = 8, lib_column = ccm_from, target_column = ccm_to, 
                      lib_sizes = n, replace = FALSE, silent = TRUE)
      ccm_matrix[ccm_from, ccm_to] <- out_temp$rho
    }
  }
  
  surr_temp_twin <- make_surrogate_data(S12T$temp, method = "twin", num_surr = num_surr)
  surr_species2_twin <- make_surrogate_data(S12T$species2, method = "twin", num_surr = num_surr)
  rho_surr1_twin <- data.frame(temp = numeric(num_surr), species2 = numeric(num_surr))
  
  for (i in 1:num_surr) {
    rho_surr1_twin$temp[i] <- ccm(cbind(S12T$species1, surr_temp_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                  replace = FALSE)$rho
    
    rho_surr1_twin$species2[i] <- ccm(cbind(S12T$species1, surr_species2_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                      replace = FALSE)$rho
  }
  
  ### P-value for rejecting the hypothesis that seasonality both causes temp and species 1
  Pval_seasonCauses1 = (sum(ccm_matrix['species1','temp'] < rho_surr1_twin$temp) + 1) / 
    (length(rho_surr1_twin$temp) + 1) # 0.87 
  ### conclusion: cannot say that species1 xmap temp (equivalent to temp -> species 1)
  
  ### P-value for rejecting the hypothesis that seasonality both causes species 1 and species 2
  # (more interesting here)
  Pval_21= (sum(ccm_matrix['species1','species2'] < rho_surr1_twin$species2) + 1) / 
    (length(rho_surr1_twin$species2) + 1) ## 0.24 
  ### conclusion: cannot say that species 1 xmap species 2, i.e. species 2 causes species 1
  
  ### (NB we have to check the other way around, i.e., if there is an interaction the technique can detect it. 
  ### cf. Cobey and Baskerville results.)
  
  ### Other causal direction -- we try to have species2 xmap species 1 to see if 1 causes 2
  surr_species1_twin <- make_surrogate_data(S12T$species1, method = "twin", num_surr = num_surr)
  rho_surr2_twin <- data.frame(temp = numeric(num_surr), species1 = numeric(num_surr))
  
  for (i in 1:num_surr) {
    rho_surr2_twin$temp[i] <- ccm(cbind(S12T$species2, surr_temp_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                  replace = FALSE)$rho
    
    rho_surr2_twin$species2[i] <- ccm(cbind(S12T$species2, surr_species1_twin[,i]), E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(S12T), 
                                      replace = FALSE)$rho
  }
  
  ### P-value for rejecting the hypothesis that temp does not cause species 2 fluctuations
  Pval_seasonCauses2 = (sum(ccm_matrix['species2','temp'] < rho_surr2_twin$temp) + 1) / 
    (length(rho_surr2_twin$temp) + 1) ### 0.000999001
  ### we reject the hypothesis of joint fluctuation, species 2 xmap temp
  
  ### P-value for rejecting the hypothesis that species 1 does not cause species 2 fluctuations
  Pval_12 = (sum(ccm_matrix['species2','species1'] < rho_surr2_twin$species2) + 1) / 
    (length(rho_surr1_twin$species2) + 1) ## 0.07692308 almost significant here 
  ### we find that species2 xmap species1 so that 1->2. 
  
  return(c(Pval_12,Pval_21))
  
}



nsims=100
p12=p12_noInter=rep(NA,nsims)#initializing
p21=p21_noInter=rep(NA,nsims)

### Loop over repeats
for (krep in 1:nsims){
  
krep
  
tmax=300

### Model with interactions
Y=matrix(1,nrow=tmax,ncol=2)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))
### Without interactions
Z=matrix(1,nrow=tmax,ncol=2)
Z[1,1]=abs(rnorm(1,1,1))
Z[1,2]=abs(rnorm(1,1,1))

seasonality<-2*sin(2*pi*(1:tmax)/24)	# must be enough to affect the growth rates
### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-2*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0.31*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
  Z[t+1,1] = Z[t,1]*exp(3+0.5*y1[t] - 4*Z[t,1]-0*Z[t,2] + rnorm(1,0,0.1))
  Z[t+1,2] = Z[t,2]*exp(2.1+0.5*y1[t] -0*Z[t,1]-3.1*Z[t,2] + rnorm(1,0,0.1))
}

y=log(Y)
z=log(Z)

pccm = ccm_surrogate_testing(y,y1)
p12[krep] = pccm[1]
p21[krep] = pccm[2]


pccm = ccm_surrogate_testing(z,y1)
p12_noInter[krep]=pccm[1]
p21_noInter[krep]=pccm[2]

}

#Fraction of 1->2 when present
sum(p12<0.1)/nsims
#Fraction of 2->1 when present
sum(p21<0.1)/nsims

#Fraction of 1->2 when absent
sum(p12_noInter<0.1)/nsims
#Fraction of 2->1 when absent
sum(p21_noInter<0.1)/nsims

library(sm)

par(mfrow=c(2,2))
## three groups, each of length: length(x), length(y), length(z)
group.index <- rep(1:2, c(length(p12), length(p21)))
## collect data together and use sm.density.compare()
den <- sm.density.compare(c(p12,p21), group = group.index, model = "none") #"equal" does bootstrap for comparison

group.index <- rep(1:2, c(length(p12_noInter), length(p21_noInter)))
## collect data together and use sm.density.compare()
den <- sm.density.compare(c(p12_noInter,p21_noInter), group = group.index, model = "none") #"equal" does bootstrap for comparison



#### ------------- Old analysis 08/10/2017 ------------------------------------------####
# pdf(file="CCM_2competitors_and_envDriver.pdf",height=10,width=10)
# par(mfrow=c(2,2),cex=1.25,lwd=2)
# species1_xmap_temp <- ccm(species2_species1_temp, E = 3, lib_column = "species1", 
#                       target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# temp_xmap_species1 <- ccm(species2_species1_temp, E = 3, lib_column = "temp", target_column = "species1", 
#                       lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# a_xmap_t_means <- ccm_means(species1_xmap_temp)
# t_xmap_a_means <- ccm_means(temp_xmap_species1)
# 
# par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
# plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
#      xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
# lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
# legend(x = "topleft", legend = c("species1 xmap temp", "temp xmap species1"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)
# 
# ### Checking between species2 and species1
# 
# species1_xmap_species2 <- ccm(species2_species1_temp, E = 3, lib_column = "species1", 
#                       target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# species2_xmap_species1 <- ccm(species2_species1_temp, E = 3, lib_column = "species2", target_column = "species1", 
#                       lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
# species2_xmap_species1_means <- ccm_means(species2_xmap_species1)
# 
# #par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
# plot(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)",ylim = c(0, 1.1))
# lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = "blue")
# legend(x = "topleft", legend = c("species1 xmap species2", "species2 xmap species1"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)
# 
# # Checking between species2 and temp
# 
# species2_xmap_temp <- ccm(species2_species1_temp, E = 3, lib_column = "species2", 
#                       target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# temp_xmap_species2 <- ccm(species2_species1_temp, E = 3, lib_column = "temp", target_column = "species2", 
#                       lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)
# 
# a_xmap_t_means <- ccm_means(species2_xmap_temp)
# t_xmap_a_means <- ccm_means(temp_xmap_species2)
# 
# par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
# plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
#      xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
# lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
# legend(x = "topleft", legend = c("species2 xmap temp", "temp xmap species2"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)
# dev.off()

### We recover the correct causal order for effects of the environment
### Though it looks like species 1 affects species 2 with CCM





