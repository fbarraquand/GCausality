### CP April 2019, based on FB's previous work
### Simulates the 2-species-and-a-(seasonal)-driver model with and without interactions between species, and computes CCM for them, as well as plot the graph using surrogates

rm(list=ls())
graphics.off()

library(rEDM)

ncond=500
num_surr <- 1000 #number of surrogate time series 

set.seed(42)
tmax=300
seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates

pdf(file="CCM_2species_and_envDriver_CP_centering_with_surrogates.pdf",height=10,width=15)
par(mfcol=c(2,3),cex=1.25,lwd=2,mar=c(4,4,0.5,0.5))

Y_with=array(1,dim=c(tmax,2,ncond))
Y_without=array(1,dim=c(tmax,2,ncond))
y1=array(1,dim=c(tmax,ncond))

for(kcond in 1:ncond){
Y_with[1,1,kcond]=abs(rnorm(1,1,1))
Y_with[1,2,kcond]=abs(rnorm(1,1,1))
Y_without[1,1,kcond]=Y_with[1,1,kcond]
Y_without[1,2,kcond]=Y_with[1,2,kcond]

### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1[,kcond]<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y_with[t+1,1,kcond] = Y_with[t,1,kcond]*exp(3+0.5*y1[t,kcond] - 4*Y_with[t,1,kcond]-2*Y_with[t,2,kcond] + rnorm(1,0,0.1))
  Y_with[t+1,2,kcond] = Y_with[t,2,kcond]*exp(2.1+0.5*y1[t,kcond] -0.31*Y_with[t,1,kcond]-3.1*Y_with[t,2,kcond] + rnorm(1,0,0.1))
}

for (t in 1:(tmax-1)){
  Y_without[t+1,1,kcond] = Y_without[t,1,kcond]*exp(3+0.5*y1[t,kcond] - 4*Y_without[t,1,kcond]-0*Y_without[t,2,kcond] + rnorm(1,0,0.1))
  Y_without[t+1,2,kcond] = Y_without[t,2,kcond]*exp(2.1+0.5*y1[t,kcond] -0*Y_without[t,1,kcond]-3.1*Y_without[t,2,kcond] + rnorm(1,0,0.1))
}

}

### The first row in the plot presents CCM and surrogates for 1 simulation, while the second row presents the CCM for the remaining 499 simulations, without surrogates. So, we need different script.
y_with=log(Y_with[,,1])
y_without=log(Y_without[,,1])

y1_tmp=y1[,1]-mean(y1[,1])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")


##########SP1 and TEMP
####### with interactions
simplex_predictx = simplex(species2_species1_temp_with[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_temp <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predictsp1_with, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) 

temp_xmap_species1 <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predicttemp_with, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

alpha=1
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha),
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(0,0,1,alpha))
legend(x = "topright", legend = c("species1 xmap temp, inter", "temp xmap species1, inter","species1 xmap temp, no inter", "temp xmap species1, no inter"), col = c("red", "blue",rgb(255/256,165/256,0,alpha),rgb(155/256,79/256,150/256,alpha)), lwd = 1, bty="n", cex = 0.8)
mtext("a)",side=2,line=2.5,las=2,at=1.1)

##### without interactions
simplex_output_predictx = simplex(species2_species1_temp_without[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_temp <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predictsp1_without, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE)

temp_xmap_species1 <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predicttemp_without, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)
lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = rgb(255/256,165/256,0,alpha))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(155/256,79/256,150/256,alpha))

###Add here the surrogates
### NOTE : maybe we should try seasonal again ??? 
surr_temp_twin_with <- make_surrogate_data(species2_species1_temp_with$temp, method = "twin", num_surr = num_surr)
surr_species1_twin_with <- make_surrogate_data(species2_species1_temp_with$species1, method = "twin", num_surr = num_surr)
rho_surr1_twin_with <- list(temp =  matrix(NA,nrow=num_surr,ncol=20), species1 =  matrix(NA,nrow=num_surr,ncol=20))

surr_temp_twin_without <- make_surrogate_data(species2_species1_temp_without$temp, method = "twin", num_surr = num_surr)
surr_species1_twin_without <- make_surrogate_data(species2_species1_temp_without$species1, method = "twin", num_surr = num_surr)
rho_surr1_twin_without <- list(temp = matrix(NA,nrow=num_surr,ncol=20), species1 =  matrix(NA,nrow=num_surr,ncol=20))

for (i in 1:num_surr) {
  rho_surr1_twin_with$temp[i,] <- ccm_means(ccm(cbind(species2_species1_temp_with$species1, surr_temp_twin_with[,i]), E =  lag_order_inter_CCM_predictsp1_with, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_with$species1[i,] <- ccm_means(ccm(cbind(surr_species1_twin_with[,i], species2_species1_temp_with$temp), E =  lag_order_inter_CCM_predicttemp_with, lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_without$temp[i,] <- ccm_means(ccm(cbind(species2_species1_temp_without$species1, surr_temp_twin_without[,i]), E =  lag_order_inter_CCM_predictsp1_without, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_without$species1[i,] <- ccm_means(ccm(cbind(surr_species1_twin_without[,i], species2_species1_temp_without$temp), E =  lag_order_inter_CCM_predicttemp_without, lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),replace = FALSE))$rho
}

surr1_max=apply(rho_surr1_twin_with$temp,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(1,0,0,alpha),lty=2)
pval_surr=sum(a_xmap_t_means$rho<rho_surr1_twin$temp) /num_surr



surr1_max=apply(rho_surr1_twin_without$temp,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col= rgb(255/256,165/256,0,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_with$species1,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(0,0,1,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_without$species1,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(155/256,79/256,150/256,alpha),lty=2)


#### Compute the CCM for the remaining 499 simulations, without surrogates
for(kcond in 2:ncond){
alpha=0.1
y_with=log(Y_with[,,kcond])
y_without=log(Y_without[,,kcond])

y1_tmp=y1[,kcond]-mean(y1[,kcond])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")


#### with interactions
simplex_output_predictx = simplex(species2_species1_temp_with[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


species1_xmap_temp <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predictsp1_with, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) 

temp_xmap_species1 <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predicttemp_with, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)


if(kcond==2){
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha),
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
mtext("d)",side=2,line=2.5,las=2,at=1.1)
}else{
lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha))
}
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(0,0,1,alpha))

#### without interactions
simplex_output_predictx = simplex(species2_species1_temp_without[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


species1_xmap_temp <- ccm(species2_species1_temp_without, E = lag_order_inter_CCM_predictsp1_without, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_without, E = lag_order_inter_CCM_predicttemp_without, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = rgb(255/256,165/256,0,alpha))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(155/256,79/256,150/256,alpha))
}


############SP1 and SP2
y_with=log(Y_with[,,1])
y_without=log(Y_without[,,1])

y1_tmp=y1[,1]-mean(y1[,1])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")

#### with interactions
simplex_output_predictx = simplex(species2_species1_temp_with[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_species2 <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predictsp1_with, lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predictsp2_with, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

alpha=1
plot(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = rgb(1,0,0,alpha), xlab = "Library Size", ylab = "",ylim = c(0, 1.1),yaxt="n")
legend(x = "topright", legend = c("species1 xmap species2,inter", "species2 xmap species1,inter","species1 xmap species2,no inter", "species2 xmap species1,no inter"), col = c("red", "blue",rgb(255/256,165/256,0,alpha),rgb(155/256,79/256,150/256,alpha)), lwd = 1,bty="n", cex = 0.8)
mtext("b)",side=2,line=2.5,las=2,at=1.1)
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = rgb(0,0,1,alpha))

#### without interactions
simplex_output_predictx = simplex(species2_species1_temp_without[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_species2 <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predictsp1_without , lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predictsp2_without, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

lines(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), col =  rgb(255/256,165/256,0,alpha))
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = rgb(155/256,79/256,150/256,alpha))


###Add here the surrogates 
surr_species2_twin_with <- make_surrogate_data(species2_species1_temp_with$species2, method = "twin", num_surr = num_surr)
surr_species1_twin_with <- make_surrogate_data(species2_species1_temp_with$species1, method = "twin", num_surr = num_surr)
rho_surr1_twin_with <- list(species1 =  matrix(NA,nrow=num_surr,ncol=20), species2 =  matrix(NA,nrow=num_surr,ncol=20))

surr_species2_twin_without <- make_surrogate_data(species2_species1_temp_without$species2, method = "twin", num_surr = num_surr)
surr_species1_twin_without <- make_surrogate_data(species2_species1_temp_without$species1, method = "twin", num_surr = num_surr)
rho_surr1_twin_without <- list(species1 = matrix(NA,nrow=num_surr,ncol=20), species2 =  matrix(NA,nrow=num_surr,ncol=20))

for (i in 1:num_surr) {
  rho_surr1_twin_with$species2[i,] <- ccm_means(ccm(cbind(species2_species1_temp_with$species1, surr_species2_twin_with[,i]), E =  lag_order_inter_CCM_predictsp1_with, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_with$species1[i,] <- ccm_means(ccm(cbind(surr_species1_twin_with[,i], species2_species1_temp_with$species2), E =  lag_order_inter_CCM_predictsp2_with, lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_without$species2[i,] <- ccm_means(ccm(cbind(species2_species1_temp_without$species1, surr_species2_twin_without[,i]), E =  lag_order_inter_CCM_predictsp1_without, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),replace = FALSE))$rho
  rho_surr1_twin_without$species1[i,] <- ccm_means(ccm(cbind(surr_species1_twin_without[,i], species2_species1_temp_without$species2), E =  lag_order_inter_CCM_predictsp2_without , lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),replace = FALSE))$rho
}

surr1_max=apply(rho_surr1_twin_with$species2,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(1,0,0,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_without$species2,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col= rgb(255/256,165/256,0,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_with$species1,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(0,0,1,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_without$species1,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(155/256,79/256,150/256,alpha),lty=2)

for(kcond in 2:ncond){
alpha=0.1
y_with=log(Y_with[,,kcond])
y_without=log(Y_without[,,kcond])

y1_tmp=y1[,kcond]-mean(y1[,kcond])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")

simplex_output_predictx = simplex(species2_species1_temp_with[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_species2 <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predictsp1_with, lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predictsp2_with, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)


if(kcond==2){
plot(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = rgb(1,0,0,alpha),
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
mtext("e)",side=2,line=2.5,las=2,at=1.1)
}else{
lines(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = rgb(1,0,0,alpha))
}

lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = rgb(0,0,1,alpha))

simplex_output_predictx = simplex(species2_species1_temp_without[,"species1"],E=1:10)
 lag_order_inter_CCM_predictsp1_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


species1_xmap_species2 <- ccm(species2_species1_temp_without, E = lag_order_inter_CCM_predictsp1_without, lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_without, E = lag_order_inter_CCM_predictsp2_without, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

lines(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), col =  rgb(255/256,165/256,0,alpha))
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = rgb(155/256,79/256,150/256,alpha))
}
################SP2 and TEMP
y_with=log(Y_with[,,1])
y_without=log(Y_without[,,1])

y1_tmp=y1[,1]-mean(y1[,1])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")

#### with interactions
simplex_output_predictx = simplex(species2_species1_temp_with[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species2_xmap_temp <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predictsp2_with, lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

temp_xmap_species2 <- ccm(species2_species1_temp_with, E = lag_order_inter_CCM_predicttemp_with, lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t_means <- ccm_means(species2_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species2)

alpha=1
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha),
     xlab = "Library Size", ylab = "",yaxt="n", ylim = c(0, 1.1))
mtext("c)",side=2,line=2.5,las=2,at=1.1)
legend(x = "topright", legend = c("species2 xmap temp, inter", "temp xmap species2, inter","species2 xmap temp, no inter", "temp xmap species2, no inter"), col = c("red", "blue",rgb(255/256,165/256,0,alpha),rgb(155/256,79/256,150/256,alpha)), lwd = 1, bty="n", cex = 0.8)

lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(0,0,1,alpha))


#### without interactions
simplex_output_predictx = simplex(species2_species1_temp_without[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]

species1_xmap_temp <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predictsp2_without, lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predicttemp_without, lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = rgb(255/256,165/256,0,alpha))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(155/256,79/256,150/256,alpha))


###Add here the surrogates 
surr_temp_twin_with <- make_surrogate_data(species2_species1_temp_with$temp, method = "twin", num_surr = num_surr)
surr_species2_twin_with <- make_surrogate_data(species2_species1_temp_with$species2, method = "twin", num_surr = num_surr)
rho_surr1_twin_with <- list(temp =  matrix(NA,nrow=num_surr,ncol=20), species2 =  matrix(NA,nrow=num_surr,ncol=20))

surr_temp_twin_without <- make_surrogate_data(species2_species1_temp_without$temp, method = "twin", num_surr = num_surr)
surr_species2_twin_without <- make_surrogate_data(species2_species1_temp_without$species2, method = "twin", num_surr = num_surr)
rho_surr1_twin_without <- list(temp = matrix(NA,nrow=num_surr,ncol=20), species2 =  matrix(NA,nrow=num_surr,ncol=20))

for (i in 1:num_surr) {
  rho_surr1_twin_with$temp[i,] <- ccm_means(ccm(cbind(species2_species1_temp_with$species2, surr_temp_twin_with[,i]), E =  lag_order_inter_CCM_predictsp2_with, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_with$species2[i,] <- ccm_means(ccm(cbind(surr_species2_twin_with[,i], species2_species1_temp_with$temp), E =  lag_order_inter_CCM_predicttemp_with, lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_without$temp[i,] <- ccm_means(ccm(cbind(species2_species1_temp_without$species2, surr_temp_twin_without[,i]), E =  lag_order_inter_CCM_predictsp2_without, lib_column = 1, target_column = 2, lib_sizes = seq(5, 100, by = 5),
                           replace = FALSE))$rho
  rho_surr1_twin_without$species2[i,] <- ccm_means(ccm(cbind(surr_species2_twin_without[,i], species2_species1_temp_without$temp), E =  lag_order_inter_CCM_predicttemp_without, lib_column = 2, target_column = 1, lib_sizes = seq(5, 100, by = 5),replace = FALSE))$rho
}

surr1_max=apply(rho_surr1_twin_with$temp,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(1,0,0,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_without$temp,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col= rgb(255/256,165/256,0,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_with$species2,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(0,0,1,alpha),lty=2)

surr1_max=apply(rho_surr1_twin_without$species2,2,function(x) quantile(x,0.95))
lines(a_xmap_t_means$lib_size,surr1_max,col=rgb(155/256,79/256,150/256,alpha),lty=2)

for(kcond in 2:ncond){
alpha=0.1
alpha=0.1
y_with=log(Y_with[,,kcond])
y_without=log(Y_without[,,kcond])

y1_tmp=y1[,kcond]-mean(y1[,kcond])

y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])

species2_species1_temp_with=data.frame(1:tmax,y_with,y1_tmp)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1_tmp)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")


simplex_output_predictx = simplex(species2_species1_temp_with[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_with = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_with[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_with = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


species1_xmap_temp <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predictsp2_with , lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE)

temp_xmap_species1 <- ccm(species2_species1_temp_with, E =  lag_order_inter_CCM_predicttemp_with , lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

if(kcond==2){
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha),
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
mtext("f)",side=2,line=2.5,las=2,at=1.1)
}else{
lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = rgb(1,0,0,alpha))
}
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(0,0,1,alpha))

#### without interactions
simplex_output_predictx = simplex(species2_species1_temp_without[,"species2"],E=1:10)
 lag_order_inter_CCM_predictsp2_without = simplex_output_predictx$E[which(simplex_output_predictx$rho==max(simplex_output_predictx$rho))]

simplex_output_predicty = simplex(species2_species1_temp_without[,"temp"],E=1:10)
 lag_order_inter_CCM_predicttemp_without = simplex_output_predicty$E[which(simplex_output_predicty$rho==max(simplex_output_predicty$rho))]


species1_xmap_temp <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predictsp2_without , lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_without, E =  lag_order_inter_CCM_predicttemp_without , lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = rgb(255/256,165/256,0,alpha))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = rgb(155/256,79/256,150/256,alpha))
}

dev.off()
