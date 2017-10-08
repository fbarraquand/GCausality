### FB 08/10/2017 -- adapted from Nantes et Irstea presentations though
### Uses a strongly seasonal driver. A logical case for conditional GC. 
library(vars)

### First case with interactions

set.seed(42)

tmax=300
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
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-2*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0.31*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
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
  Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-0*Y[t,2] + rnorm(1,0,0.1))
  Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.1))
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


### Let's use another simulation with just a little more noise

# Y=matrix(1,nrow=tmax,ncol=2)
# Y[1,1]=abs(rnorm(1,1,1))
# Y[1,2]=abs(rnorm(1,1,1))
# 
# for (t in 1:(tmax-1)){
#   Y[t+1,1] = Y[t,1]*exp(3+0.5*y1[t] - 4*Y[t,1]-0*Y[t,2] + rnorm(1,0,0.3))
#   Y[t+1,2] = Y[t,2]*exp(2.1+0.5*y1[t] -0*Y[t,1]-3.1*Y[t,2] + rnorm(1,0,0.3))
# }
# y=log(Y)
# 
# varcompet2ter<-VAR(y, type="none",exogen=y1,lag.max=10,ic="SC")
# summary(varcompet2ter)
# causality(varcompet2ter,cause="y1") #p-value 6.263e-10
# causality(varcompet2ter,cause="y2") #0.5799

### So it is possible to craft examples that are difficult to infer with GC 
### --- still we may need many simulations and parameter values for a good comparison of CCM and GC

### NB might well be possible that we would get better results by restricting the fit to a VAR1
### But we have no theoretical reason to do that

############## Two things to check though!!! ###################
### I haven't tested pairwise GC on these examples. Perhaps it works better than conditional GC? 
### And I haven't CCM this, on the other hand. 
### They say you don't need necessarily corrections for seasonality in the 2012 paper -- or do we? cf Deyle PNAS Influenza

############# Pairwise GC (without conditioning on the driver) & CCM to test absolutely first... 

varcompet_noExo=VAR(y, type="none",lag.max=10,ic="SC") ## Direct pairwise GC
summary(varcompet_noExo)
causality(varcompet_noExo,cause="y1") #p-value 0.03633
causality(varcompet_noExo,cause="y2") #0.00000

### OK we do have a problem -- but what does CCM say? 

########## CCM analysis of this dataset ######################


### Merging to a dataframe

species2_species1_temp=data.frame(1:tmax,y,y1)
names(species2_species1_temp)=c("time","species1","species2","temp")

pdf(file="CCM_2competitors_and_envDriver.pdf",height=10,width=10)
par(mfrow=c(2,2),cex=1.25)
species1_xmap_temp <- ccm(species2_species1_temp, E = 3, lib_column = "species1", 
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

temp_xmap_species1 <- ccm(species2_species1_temp, E = 3, lib_column = "temp", target_column = "species1", 
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("species1 xmap temp", "temp xmap species1"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Checking between species2 and species1

species1_xmap_species2 <- ccm(species2_species1_temp, E = 3, lib_column = "species1", 
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp, E = 3, lib_column = "species2", target_column = "species1", 
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)",ylim = c(0, 1.1))
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = "blue")
legend(x = "topleft", legend = c("species1 xmap species2", "species2 xmap species1"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

# Checking between species2 and temp

species2_xmap_temp <- ccm(species2_species1_temp, E = 3, lib_column = "species2", 
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

temp_xmap_species2 <- ccm(species2_species1_temp, E = 3, lib_column = "temp", target_column = "species2", 
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t_means <- ccm_means(species2_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species2)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("species2 xmap temp", "temp xmap species2"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)
dev.off()

### We recover the correct causal order for effects of the environment
### Though it looks like species 1 affects species 2 with CMM

