###############################################################################################
### FB 11/06/2018 -- Using Simone by Chiquet et al. on simulated ecological-style VAR1 data
### Interaction network inference using sparse method on time series ##########################
###############################################################################################

rm(list = ls())
library(simone)
set.seed(5)

################## Simulating a VAR(1) model for now ###########################
#adjacency matrix
causality_matrix = rbind(c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,0,0,0,0,0,0,0),
                         c(1,1,1,1,0,0,0,0,0,0),
                         c(0,0,0,1,1,1,1,0,0,0),
                         c(0,0,0,0,1,1,1,0,0,0),
                         c(0,0,0,0,1,1,1,0,0,0),
                         c(0,0,0,0,0,0,0,1,1,1),
                         c(0,0,0,0,0,0,0,1,1,1),
                         c(0,0,0,0,0,0,0,1,1,1))

spectral_radius=1.5 ### what would be bad for the spectral radius
while(spectral_radius>1){
B=causality_matrix*rnorm(100,-0.01,0.1)+0.8*diag(10)
## Your typical MAR(1) matrix
spectral_radius = max(abs(eigen(B)$values)) #stable
}

### Simulation for 300 timesteps
ysim=matrix(0,nrow = 300,ncol=10)
ysim[1,] = rnorm(10,0.5,0.1)
for (time_index in 1:299){
  ysim[time_index+1,] = B %*% ysim[time_index,] + rnorm(10,0,0.1)
}
matplot(ysim)
matlines(ysim) #perfect!

### Let's try SIMONE //depends on mixer
#demo(simone_timeCourse)

### First no clusters
res.noclust<-simone(ysim,type="time-course") #OK
plot(res.noclust) #indeed
res.noclust$clusters

### Obviously we should not go above about 50 edges, more like 25
g.no<-getNetwork(res.noclust,30) 
### Second argument controls the max number of edges 
### ?getNetwork
plot(g.no)
plot(g.no,type="cluster")
Bhat=g.no$Theta
B
abs(B)>0
abs(Bhat)>0 #not so bad
##other diagnostics
cor(as.vector(B),as.vector(Bhat)) #not very meaningful
## only diag
cor(diag(B),diag(Bhat))
## only off-diag?

### What if we allow for more edges? 
g.no<-getNetwork(res.noclust,50) 
### Second argument controls the max number of edges 
### ?getNetwork
plot(g.no)
plot(g.no,type="cluster")
plot(g.no,type="circles")
### OK, that's overpopulated in edges -- we need to avoid this
Bhat=g.no$Theta
B
abs(B)>0
abs(Bhat)>0 #not so bad but quite a few false positives

### Let's try with AIC/BIC 
#! note that the default "selection" in getNetwork() should not be used because it is set to length(object$clusters) = number of nodes
g.no<-getNetwork(res.noclust,selection = "BIC") 
plot(g.no,type="cluster")
plot(g.no,type="circles")
# Not bad!!
g.no$weights

g.no.aic<-getNetwork(res.noclust,selection = "AIC") 
plot(g.no.aic,type="cluster")
plot(g.no.aic,type="circles")
# overparameterized, could have been expected


### ## Let's try with clustering then #####
ctrl<-setOptions(clusters.crit = "BIC")
res.cl.new = simone(ysim,type="time-course",clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl.new,30)
g.cl
plot(g.cl)
plot(g.cl,type="circles")
### Is this the result of a hub-leaf structure or is it a modular constraint on the nodes?
g.cl$weights

### Now use BIC for full model selection 
g.cl=getNetwork(res.cl.new,selection="BIC")
g.cl
plot(g.cl)
plot(g.cl,type="circles")

### Compare the two networks, with and without clusters
plot(g.no,g.cl)
plot(g.no,g.cl,type="overlap")

# let's constrain a bit the cluster numbers to reasonable bounds and increase max number of edges
ctrl<-setOptions(clusters.qmin = 2, clusters.qmax = 5)
res.cl.new = simone(ysim,type="time-course",clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl.new,"BIC") 
g.cl
plot(g.cl,type="circles")
#still too many edges

### Compare the two networks
plot(g.no,g.cl)
plot(g.no,g.cl,type="overlap")
### Makes little difference


############## Stuff that did not work // bad syntax -- here for memory #############
#### If this is a hub/non-hub structure by default, try to condition weights
### weights = causality_matrix*0.5+0.25
### mod.net = simone(ysim,type="time-course",weights = weights, control=ctrl)
# does not work, weights is only an output
################ Try with InferEdges? ###################
# InferEdges(ysim,penalty = causality_matrix)
#InferEdges <- function(X, penalty, ctrl.inf = OptInference())

#### I have to look up that one directly in the tar.gz
### After modified simone() to include initial guess
#install.packages("/home/frederic/Documents/GC_nonlinear_communityDynamics/SimoneFiles/simone.tar.gz",repos=NULL,type="source")
apriori.matrix = causality_matrix*0.5+0.25
ctrl<-setOptions(clusters.crit = "BIC",initial.guess=apriori.matrix)
mod.net = simone(ysim,type="time-course",control=ctrl)
g.mod<- getNetwork(mod.net,"BIC") 
plot(g.mod)
g.mod$Theta

plot(g.mod,g.cl) # pareil, pas d'effet

### stronger a priori!!
apriori.matrix = causality_matrix ## zeroes and ones
ctrl<-setOptions(clusters.crit = "BIC",initial.guess=apriori.matrix)
mod.net = simone(ysim,type="time-course",control=ctrl)
g.mod<- getNetwork(mod.net,"BIC") 
plot(g.mod) # always the same

### Compute false positives and negatives
Pos = (abs(g.mod$Theta)>0)
#Neg = (abs(g.mod$Theta)==0)
TruePos =  (abs(B)>0)
#TrueNeg =  (abs(B)==0)

TP = sum((TruePos == TRUE) & (Pos == TRUE))
FN = sum((TruePos == TRUE) & (Pos == FALSE))
FP = sum((TruePos == FALSE) & (Pos == TRUE))
TN = sum((TruePos == FALSE) & (Pos == FALSE))

Precision = TP/(TP+FP)  #True positives over tested positives (percentage positives correct)
Recall = TP/(TP+FN)     #True positives over relevant elements (power to find the true positives)

Precision
Recall
### I can make distributions for these ones

#### Let's try to do that with "weights" -- I needed to change the functions for that
#apriori.matrix = (causality_matrix*0.5+0.25)
apriori.matrix = (causality_matrix*0.8+0.01)
### Looks like if we interpret these as probabilities this has the opposite effects
### Let's consider that high weight means low connection probability
apriori.matrix = (1-causality_matrix*0.9)

ctrl<-setOptions(clusters.crit = "BIC",weights=apriori.matrix)
mod.net = simone(ysim,type="time-course",control=ctrl,clustering=FALSE)
g.mod<- getNetwork(mod.net,"BIC") 
plot(g.mod)
g.mod$Theta
plot(g.mod,g.cl) # now the last network is best

### Compute false positives and negatives
Pos = (abs(g.mod$Theta)>0)
TruePos =  (abs(B)>0)
TP = sum((TruePos == TRUE) & (Pos == TRUE))
FN = sum((TruePos == TRUE) & (Pos == FALSE))
FP = sum((TruePos == FALSE) & (Pos == TRUE))
TN = sum((TruePos == FALSE) & (Pos == FALSE))
Precision = TP/(TP+FP)  
Recall = TP/(TP+FN)     
Precision
Recall


############# Try that with a nonlinear data-generating model #############


####################################################################################
########### 10 species stochastic and (very) nonlinear competition model
####################################################################################

### Two sets of 3 species connected by 1, plus 1 set of 3 unconnected species
tmax=300 ## We are very optimistic, Sugihara et al. also used 3000. 
nspecies=10
sigma=0.3 # Sigma^2= 0.1
Y=matrix(1,nrow=tmax,ncol=nspecies)
for (i in 1:nspecies){
  Y[1,i]=abs(rnorm(1,1,1))}

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(4-4*Y[t,1]-2*Y[t,2]-0.4*Y[t,3]+rnorm(1,0,sigma))
  Y[t+1,2] = Y[t,2]*exp(3.1-0.31*Y[t,1]-3.1*Y[t,2]-0.93*Y[t,3]+rnorm(1,0,sigma))
  Y[t+1,3] = Y[t,3]*exp(0.12 + 0.636*Y[t,1]+0.636*Y[t,2]-2.12*Y[t,3]+rnorm(1,0,sigma))
  ### Connecting the two sets by being influenced by all three and having an effect on 5
  Y[t+1,4] = Y[t,4]*exp(3-0.111*Y[t,1]-0.111*Y[t,2]+0.131*Y[t,3]-3.8*Y[t,4]+rnorm(1,0,sigma))
  ### Influenced by 4 and also 5,6,7
  Y[t+1,5] = Y[t,5]*exp(3-2*Y[t,5]-2*Y[t,6]-0.4*Y[t,7]+0.5*Y[t+1,4]+ rnorm(1,0,sigma))
  ### 
  Y[t+1,6] = Y[t,6]*exp(3.1-0.31*Y[t,5]-3.1*Y[t,6]-0.93*Y[t,7]+rnorm(1,0,sigma))
  Y[t+1,7] = Y[t,7]*exp(0.8 + 0.636*Y[t,5]+0.636*Y[t,6]-2.12*Y[t,7]+rnorm(1,0,sigma))
  Y[t+1,8] = Y[t,8]*exp(4-4*Y[t,8]-2*Y[t,9]-0.4*Y[t,10]+rnorm(1,0,sigma))
  Y[t+1,9] = Y[t,9]*exp(3.1-0.31*Y[t,8]-3.1*Y[t,9]-0.93*Y[t,10]+rnorm(1,0,sigma))
  Y[t+1,10] = Y[t,10]*exp(0.12 + 0.636*Y[t,8]+0.636*Y[t,9]-2.12*Y[t,10]+rnorm(1,0,sigma))
  
}

y=log(Y)
y=y-mean(y)

matplot(y)
matlines(y)
### Number 7 is annoying - but let's analyze that anyway

### First no clusters
res.noclust<-simone(y,type="time-course") #OK
plot(res.noclust) #indeed
res.noclust$clusters

### Try with BIC
g.no<-getNetwork(res.noclust,selection="BIC") 
### ?getNetwork
plot(g.no)
plot(g.no,type="cluster")
Bhat=g.no$Theta
B
abs(B)>0
abs(Bhat)>0 #not so bad

#is there any improvement if I manage to add clusters? 
res.cl.new = simone(y,type="time-course",clustering=TRUE)
g.cl=getNetwork(res.cl.new,"BIC") 
g.cl
plot(g.cl,type="circles")
#too many edges?
Bhat=g.cl$Theta
B
abs(B)>0
abs(Bhat)>0 #not so bad


### Compare the two networks
plot(g.no,g.cl)
plot(g.no,g.cl,type="overlap")
### 

### Now with weights
apriori.matrix = (1-causality_matrix*0.9)
ctrl<-setOptions(clusters.crit = "BIC",weights=apriori.matrix)
mod.net = simone(y,type="time-course",control=ctrl,clustering=FALSE)
g.mod<- getNetwork(mod.net,"BIC") 
plot(g.mod)
g.mod$Theta
plot(g.mod,g.cl) # now the last network is best

### Compute false positives and negatives
Pos = (abs(g.mod$Theta)>0)
TruePos =  (abs(B)>0)
TP = sum((TruePos == TRUE) & (Pos == TRUE))
FN = sum((TruePos == TRUE) & (Pos == FALSE))
FP = sum((TruePos == FALSE) & (Pos == TRUE))
TN = sum((TruePos == FALSE) & (Pos == FALSE))
Precision = TP/(TP+FP)  
Recall = TP/(TP+FN)     
Precision
Recall

