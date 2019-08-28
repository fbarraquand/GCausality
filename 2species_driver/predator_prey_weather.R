

library(rEDM)

### We have predator-prey-weather

### Environmental weather variable
u<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=100,sd=sqrt(0.5) )
plot(u)
spectrum(u,method="ar")

### Predator-prey system
n<-100 # sample size
# Parameters
rmax<-1
alpha<-0.1 # intraspe. compet
epsilon<-0.1 # conversion efficiency
attack_rate<-0.25 #discovery rate of prey (LV modelling)
mu<-0.05 #predator death rate
sigma=0.25 # noise term
#
Nt=matrix(1,n,2) #initializing and defining abundance vector
for (t in 1:(n-1)){
	Nt[t+1,1]=Nt[t,1]*exp(rmax+0.1*u[t]+rnorm(1,0,sigma) -alpha*Nt[t,1]-attack_rate*Nt[t,2])
	Nt[t+1,2]=Nt[t,2]*exp(epsilon*attack_rate*Nt[t,1]-mu+rnorm(1,0,sigma))
}


### Merging to a dataframe

pred_prey_temp=data.frame(1:n,Nt,u)
names(pred_prey_temp)=c("time","prey","pred","temp")


prey_xmap_temp <- ccm(pred_prey_temp, E = 3, lib_column = "prey", 
                        target_column = "temp", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

temp_xmap_prey <- ccm(pred_prey_temp, E = 3, lib_column = "temp", target_column = "prey", 
                        lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

a_xmap_t_means <- ccm_means(prey_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_prey)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap temp", "temp xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Checking between predator and prey

prey_xmap_pred <- ccm(pred_prey_temp, E = 3, lib_column = "prey", 
                        target_column = "pred", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

pred_xmap_prey <- ccm(pred_prey_temp, E = 3, lib_column = "pred", target_column = "prey", 
                        lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

prey_xmap_pred_means <- ccm_means(prey_xmap_pred)
pred_xmap_prey_means <- ccm_means(pred_xmap_prey)

#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(prey_xmap_pred_means$lib_size, pmax(0, prey_xmap_pred_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1.5))
lines(pred_xmap_prey_means$lib_size, pmax(0, pred_xmap_prey_means$rho), col = "blue")
legend(x = "topleft", legend = c("prey xmap pred", "pred xmap prey"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

# Checking between predator and temp

pred_xmap_temp <- ccm(pred_prey_temp, E = 3, lib_column = "pred", 
                      target_column = "temp", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

temp_xmap_pred <- ccm(pred_prey_temp, E = 3, lib_column = "temp", target_column = "pred", 
                      lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

a_xmap_t_means <- ccm_means(pred_xmap_temp)
t_xmap_a_means <- ccm_means(temp_xmap_pred)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("pred xmap temp", "temp xmap pred"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

# relatively low - still one might conclude that predators are caused by temp fluctuations, a little a least? 

# Should compute Granger causality on the same sample to compare. 

