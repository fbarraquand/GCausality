rm(list=ls())
graphics.off()

library(rEDM)


set.seed(42)

tmax=300
Y_with=matrix(1,nrow=tmax,ncol=2)
Y_without=matrix(1,nrow=tmax,ncol=2)
Y_with[1,1]=abs(rnorm(1,1,1))
Y_with[1,2]=abs(rnorm(1,1,1))
Y_without[1,1]=Y_with[1,1]
Y_without[1,2]=Y_with[1,2]

seasonality<-2*sin(2*pi*(1:tmax)/24)    # must be enough to affect the growth rates
### Environmental variables
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.5) )
###y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1<-seasonality+y1noise
##y2<-seasonality+y2noise

for (t in 1:(tmax-1)){
  Y_with[t+1,1] = Y_with[t,1]*exp(3+0.5*y1[t] - 4*Y_with[t,1]-2*Y_with[t,2] + rnorm(1,0,0.1))
  Y_with[t+1,2] = Y_with[t,2]*exp(2.1+0.5*y1[t] -0.31*Y_with[t,1]-3.1*Y_with[t,2] + rnorm(1,0,0.1))
}
y_with=log(Y_with)

for (t in 1:(tmax-1)){
  Y_without[t+1,1] = Y_without[t,1]*exp(3+0.5*y1[t] - 4*Y_without[t,1]-0*Y_without[t,2] + rnorm(1,0,0.1))
  Y_without[t+1,2] = Y_without[t,2]*exp(2.1+0.5*y1[t] -0*Y_without[t,1]-3.1*Y_without[t,2] + rnorm(1,0,0.1))
}
y_without=log(Y_without)


###########Should we center ?
if(1==0){
y1=y1-mean(y1)
y_with[,1]=y_with[,1]-mean(y_with[,1])
y_with[,2]=y_with[,2]-mean(y_with[,2])
y_without[,1]=y_without[,1]-mean(y_without[,1])
y_without[,2]=y_without[,2]-mean(y_without[,2])
}
#################just for now, we can still remove this
#It was expected, but centering does not change anything

species2_species1_temp_with=data.frame(1:tmax,y_with,y1)
species2_species1_temp_without=data.frame(1:tmax,y_without,y1)

names(species2_species1_temp_with)=c("time","species1","species2","temp")
names(species2_species1_temp_without)=c("time","species1","species2","temp")

pdf(file="CCM_2species_and_envDriver_CP_nocentering.pdf",height=5,width=15)
par(mfrow=c(1,3),cex=1.25,lwd=2,mar=c(4,4,0.5,0.5))

##########SP1 and TEMP
species1_xmap_temp <- ccm(species2_species1_temp_with, E = 3, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_with, E = 3, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red",
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.1))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topright", legend = c("species1 xmap temp", "temp xmap species1"), col = c("red", "blue"), lwd = 1, bty="n", cex = 0.8)
legend("topleft","a)",cex=0.8,bty="n")


species1_xmap_temp <- ccm(species2_species1_temp_without, E = 3, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_without, E = 3, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = "red",lty=2)
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue",lty=2)



############SP1 and SP2
species1_xmap_species2 <- ccm(species2_species1_temp_with, E = 3, lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_with, E = 3, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

plot(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "",ylim = c(0, 1.1),yaxt="n")
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = "blue")
legend(x = "topright", legend = c("species1 xmap species2", "species2 xmap species1"), col = c("red", "blue"), lwd = 1,bty="n", cex = 0.8)
legend("topleft","b)",cex=0.8,bty="n")

species1_xmap_species2 <- ccm(species2_species1_temp_without, E = 3, lib_column = "species1",
                      target_column = "species2", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species2_xmap_species1 <- ccm(species2_species1_temp_without, E = 3, lib_column = "species2", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

species1_xmap_species2_means <- ccm_means(species1_xmap_species2)
species2_xmap_species1_means <- ccm_means(species2_xmap_species1)

lines(species1_xmap_species2_means$lib_size, pmax(0, species1_xmap_species2_means$rho), col = "red",lty=2)
lines(species2_xmap_species1_means$lib_size, pmax(0, species2_xmap_species1_means$rho), col = "blue",lty=2)

################SP2 and TEMP
species2_xmap_temp <- ccm(species2_species1_temp_with, E = 3, lib_column = "species2",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

temp_xmap_species2 <- ccm(species2_species1_temp_with, E = 3, lib_column = "temp", target_column = "species2",
                      lib_sizes = seq(5, 100, by = 5), random_libs = FALSE)

a_xmap_t_means <- ccm_means(species2_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species2)

plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red",
     xlab = "Library Size", ylab = "",yaxt="n", ylim = c(0, 1.1))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topright", legend = c("species2 xmap temp", "temp xmap species2"), col = c("red", "blue"), lwd = 1,  cex = 0.8,bty="n")
legend("topleft","c)",cex=0.8,bty="n")


species1_xmap_temp <- ccm(species2_species1_temp_without, E = 3, lib_column = "species1",
                      target_column = "temp", lib_sizes = seq(5, 100, by = 5), replace=FALSE) #100 samples is the default, I'm too lazy to write it

temp_xmap_species1 <- ccm(species2_species1_temp_without, E = 3, lib_column = "temp", target_column = "species1",
                      lib_sizes = seq(5, 100, by = 5), replace=FALSE)

a_xmap_t_means <- ccm_means(species1_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_species1)

lines(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), col = "red",lty=2)
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue",lty=2)

dev.off()
