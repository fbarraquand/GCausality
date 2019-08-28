###### Analyses of the sardine - anchovy system forced by SST



library(rEDM)

### Define what we use as library (fitted indices) and predicted ones
lib <- c(1, 100)
pred <- c(201, 500)


data(sardine_anchovy_sst)
head(sardine_anchovy_sst)

anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy", 
                        target_column = "np_sst", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "anchovy", 
                        lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

### Do the same with sardines 

sardine_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "sardine", 
                        target_column = "np_sst", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

sst_xmap_sardine <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "sardine", 
                        lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

a_xmap_t_means <- ccm_means(sardine_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_sardine)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.06))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("sardine xmap SST", "SST xmap sardine"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

## does not work so well... 

### And then, there's the questions whether anchovy sardine could be mistaken, using CCM, for an interacting system without SST

### Let's analyse this differently wih VAR?

plot(sardine_anchovy_sst$year,sardine_anchovy_sst$anchovy,type="o")
lines(sardine_anchovy_sst$year,sardine_anchovy_sst$sardine,col="red")
lines(sardine_anchovy_sst$year,sardine_anchovy_sst$np_sst,col="green")

### Perhaps this is going to be difficult, but let's try

SAS<-sardine_anchovy_sst[,c(2:3)] ## NB has this data been logged before centering - this is unclear
SST = sardine_anchovy_sst[,5]

var_forced=VAR(SAS,type="none",lag.max = 15,ic="SC",exogen=SST)
summary(var_forced)
### I could also fit the full MAR

### Does not seem to be strong effects of sardine -> anchovy and anchovy -> sardine, but clearly this could tested by
### comparing nested models, either through IC or LRS
### There's a significant effect of SST only on the anchovy

### Anchovy is modelled OK but not sardine // but this seems to be the case also for CCM. 

### They made a whole in the time series -- their data in the 2012 paper looks very different as well. 
