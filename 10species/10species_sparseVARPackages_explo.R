########################################################################################################################
########### FBarraquand 18/05/2017 - GC analysis with proper model selection on nonlinear community dynamics ###########
########################################################################################################################

library("vars")

set.seed(42)

####################################################################################
########### 10 species stochastic and nonlinear competition model
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
causality_matrix
yd=data.frame(y)

VAR(y=yd, type="none",lag.max=15)
VAR(y=yd, type="none",lag.max=15,ic="SC")

IC=VARselect(y=yd, type="none",lag.max=15) ## selection by AIC (not even AICc)
IC
matplot(1:15,scale(t(IC$criteria)))

crit=scale(t(IC$criteria))
           
pdf(file="LagOrderSelection10Species.pdf",width=8,height=6)
par(cex=1.5,lwd=3)
plot(1:15,crit[,1],ylab="Information Criteria",xlab="Number of lags",type="o",ylim=c(-2,2))
### Add other
col_vec=c("black","red","green","blue")
for (i in 2:4){lines(1:15,crit[,i],type="o",col=col_vec[i+1])}
legend(2,2,legend=c("AIC","HQ","BIC","FPE"),col=col_vec,pch=18)
dev.off()

var1=VAR(y=yd, p=1,type="none")
var2=VAR(y=yd, p=2,type="none")

causality(var1,cause="X1") #aha, not very useful. 
coef(var1)$X1
  
# Better to either compute it pairwise GC or conditional GC (based on VAR coefs)
B=matrix(0,nspecies,nspecies)
p_value=matrix(0,nspecies,nspecies) #Initializing

options(digits=2)

# Loop over affected species
names_var1=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
names_var1=names(coef(var1))

i=0
for (x in names_var1)
{
  i=i+1
  y=coef(var1)[[i]]
  for (j in 1:10){
    p_value[i,j]=y[j,4]
    B[i,j]=y[j,1]
  }
}
B
round(p_value, digits=2)
B
### Printing B
library(xtable)
Btable <- xtable(B)
print(Btable, floating=FALSE, tabular.environment="pmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

### Try pairwise GC now
p_value1=matrix(0,nrow=nspecies,ncol=nspecies)
for (i in 1:nspecies){
  for (j in 1:nspecies){
    # cause first and effet later in grangertest()
    if (i !=j){
      g=grangertest(Y[,j],Y[,i],order=1)
      p_value1[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
    }
  }
}
G_1=(p_value1<0.05)
G_1
sum(G_1==causality_matrix)/length(as.vector(G_1))
G_1=(p_value1<0.01)
G_1
sum(G_1==causality_matrix)/length(as.vector(G_1))

### Pairwise GC is better than conditional because I used the full VAR
### Very interesting stuff on this here 
### https://stats.stackexchange.com/questions/208983/estimating-a-restricted-sparse-var-in-r
restrict(var1, method = "ser")
### Not so good
### What is the AIC/BIC of the real model?
### restrict(var1, method = "man", resmat = causality_matrix)
### Does not fit what I want, i.e. the real model with 0 where it should have them. 

### Try new sparseVAR package
### install.packages("sparsevar", repos = "http://cran.us.r-project.org")
library(sparsevar)
set.seed(1)
sim <- simulateVAR(N = 20, p = 2)

results <- fitVAR(as.matrix(yd),p=1)
str(results)
##plotVAR(as.matrix(yd), results)
library(xtable)
results$A ### Seems better

Around=round(as.matrix(results$A[[1]]),digits=2)
Atable <- xtable(Around)
print(Atable, floating=FALSE, tabular.environment="pmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

### Other estimates
results2 <- fitVAR(as.matrix(yd),p=1,penalty = "ENET", alpha = 1, type.measure = "mae",
                   lambda = "lambda.1se") ### Pas terrible
results2$A

results2 <- fitVAR(as.matrix(yd),p=1,penalty = "ENET", parallel = TRUE,
                   ncores = 5, alpha = 0.95, type.measure = "mae",
                   lambda = "lambda.1se") ### Pas terrible

results2$A

results2 <- fitVAR(as.matrix(yd),p=1,penalty = "ENET", parallel = TRUE,
                   ncores = 5, alpha = 0.95, type.measure = "mse",
                   lambda = "lambda.1se") ### Pire

results3 <- fitVAR(as.matrix(yd),p=1,penalty = "ENET",alpha = 0.95, type.measure = "mse",
                   lambda = "lambda.1se")
results3$A ### Same as before at least

results3 <- fitVAR(as.matrix(yd),p=1,penalty = "ENET",alpha = 0.95, type.measure = "mse")
results3$A

results4<- fitVAR(as.matrix(yd),p=1,penalty = "ENET",alpha = 0.1, nlambda=100)
results4$A ### Closer to Ridge

results4<- fitVAR(as.matrix(yd),p=1,penalty = "ENET",alpha = 1, nlambda=100)
results4$A ### LASSO

results5<- fitVAR(as.matrix(yd),p=1,penalty = "ENET",alpha = 0.95, nlambda=100)
results5$A ### Seems like the magical 0.05 value 

### Essayons le LASSO avec p=2
results6<- fitVAR(as.matrix(yd),p=2,penalty = "ENET",alpha = 1, nlambda=10)
results6$A ### 

### Let's try to use many LASSO structure? ###
### should we jackniffe or bootstrap somehow? 
yboot<-bootstrappedVAR(results)
matplot(yboot)
matlines(yboot)
matplot(yd)
matlines(yd)
### Not clear what this does, probably not what I want. 

### Let's test a bit more sparsevar
sim <- simulateVAR(N = 10, p = 1,sparsity=0.5,rho=0.1)
sim$A ### That's not a very good model

### Simulating Something better 
B=causality_matrix*rnorm(100,-0.01,0.1)+0.8*diag(10)
## Your typical MAR(1) matrix
max(abs(eigen(B)$values)) #stable
ysim=matrix(0,nrow = 300,ncol=10)
ysim[1,] = rnorm(10,0.5,0.1)
for (time_index in 1:299){
  ysim[time_index+1,] = B %*% ysim[time_index,] + rnorm(10,0,0.1)
}
matplot(ysim)
matlines(ysim) #perfect!

res_gompertz=fitVAR(ysim,p=1)
B_hat=round(as.matrix(res_gompertz$A[[1]]),digits=2)
B_1=abs(B_hat)>0.01
sum(B_1==causality_matrix)/length(as.vector(B_1))
#0.72
########

G_1=abs(Around)>0.01
sum(G_1==causality_matrix)/length(as.vector(G_1))
#0.74

### Thus performance is not just a question of nonlinearities in the underlying dynamics
### since the linear MAR(1) is not very well estimated by LASSO either

########## I may need to use MARSS to estimate this properly with blocks of zeroes ##########
########## ?

### Let's try SIMONE by Charbonnier et al. // depends on mixer
library(simone)
data(cancer)
str(cancer,max.level = 1)

attach(cancer)
boxplot(expr,las=3,cex.axis=0.6)
matplot(expr)
matlines(expr) # are these time series?
#no clustering
res.no<-simone(expr) #OK
plot(res.no) #indeed
str(expr)

g.no<-getNetwork(res.no,30)
plot(g.no)
plot(g.no,type="cluster")

## try with clustering
ctrl<-setOptions(clusters.crit = 30)
res.cl = simone(expr,clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl,30)
plot(g.cl)
plot(g.cl,type="circles")

### Compare the two networks
plot(g.no,g.cl)
plot(g.no,g.cl,type="overlap")


####### !! OK let's try this with our data!! #########

### First no clusters
res.noclust<-simone(ysim,type="time-course") #OK
plot(res.noclust) #indeed

### Obviously we should not go above about 50 edges
### There 32 true edges but above 20 this seems a bit difficult

g.no<-getNetwork(res.noclust,30)
plot(g.no)
plot(g.no,type="cluster")

## try with clustering
ctrl<-setOptions(clusters.crit = 30)
res.cl.new = simone(ysim,type="time-course",clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl.new,30)
g.cl
plot(g.cl)
plot(g.cl,type="circles")

#Let's try with less clusters
## try with clustering
ctrl<-setOptions(clusters.crit = 20)
res.cl.new = simone(ysim,type="time-course",clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl.new,20)
g.cl
plot(g.cl,type="circles")

### Compare the two networks
plot(g.no,g.cl)
plot(g.no,g.cl,type="overlap")

### Perhaps I need more clusters
ctrl<-setOptions(clusters.crit = 40)
res.cl.new = simone(ysim,type="time-course",clustering=TRUE,control=ctrl)
g.cl=getNetwork(res.cl.new,40)
g.cl

demo(simone_timeCourse)
# 

### Try https://academic.oup.com/bioinformatics/article/26/18/i517/205683
### Can't download it?

### Try stuff from Bob O'Hara?

### One idea would really be to use a bootstrapped structure. 


######### Trying out? BigVAR http://www.wbnicholson.com/BigVAR.html
######### VARX-L: Structured Regularization for Large Vector Autoregressions with Exogenous Variables
######### 

library(BigVAR)
mod1<-constructModel(as.matrix(yd),p=1,"Basic",gran=c(150,10),RVAR=FALSE,h=1,cv="Rolling",MN=FALSE,verbose=FALSE,IC=TRUE)
### Now plotting the stuff
SparsityPlot.BigVAR.results(cv.BigVAR(mod1))
### Hmm - not very convincing

mod1<-constructModel(as.matrix(yd),p=1,"Basic",gran=c(50,10),RVAR=FALSE,h=1,cv="LOO",MN=FALSE,verbose=FALSE,IC=TRUE)
### Now plotting the stuff
SparsityPlot.BigVAR.results(cv.BigVAR(mod1))
### Looks the penalty is somehow stronger? 

mod1<-constructModel(as.matrix(yd),"SparseLag",gran=c(150,1),RVAR=TRUE,h=1,cv="LOO",verbose=TRUE,IC=TRUE)
### Now plotting the stuff
cvBV=cv.BigVAR(mod1)
SparsityPlot.BigVAR.results(cvBV)
### Hmm - not very convincing

mod2<-constructModel(as.matrix(yd),p=3,struct="HVARELEM",gran=c(50,1),RVAR=FALSE,h=1,cv="Rolling",MN=TRUE,verbose=TRUE,VARX=list(),IC=TRUE)
### Now plotting the stuff
cvBV2=cv.BigVAR(mod2)
SparsityPlot.BigVAR.results(cvBV2)
#### MN = True seems to give better results

mod2<-constructModel(ysim,p=2,struct="HVARC",gran=c(50,1),RVAR=FALSE,h=1,cv="Rolling",MN=TRUE,verbose=TRUE,VARX=list(),IC=TRUE)
### Now plotting the stuff
cvBV2=cv.BigVAR(mod2)
cvBV2
SparsityPlot.BigVAR.results(cvBV2)
### Better?

mod2<-constructModel(ysim,p=2,struct="HVARELEM",gran=c(50,1),RVAR=FALSE,h=1,cv="Rolling",MN=TRUE,verbose=TRUE,VARX=list(),IC=TRUE)
### Now plotting the stuff
cvBV2=cv.BigVAR(mod2)
cvBV2
SparsityPlot.BigVAR.results(cvBV2)

mod2<-constructModel(ysim,p=2,struct="SparseOO",gran=c(50,1),RVAR=FALSE,h=1,cv="Rolling",MN=TRUE,verbose=TRUE,VARX=list(),IC=TRUE)
### Now plotting the stuff
cvBV2=cv.BigVAR(mod2)
cvBV2
SparsityPlot.BigVAR.results(cvBV2)
### Not good

LSAIC <- VARXFit(ysim, 10, "AIC", NULL)
LSAIC$Bhat

############ Now try https://github.com/ineswilms/bigtime
#install.packages("devtools")
library("devtools")
#devtools::install_github("ineswilms/bigtime")
library('bigtime')

VARL1 <- sparseVAR(Y=ysim, VARpen="L1") # default forecast horizon is h=1
par(mfrow=c(1, 1))
LhatL1 <- lagmatrix(fit=VARL1, model="VAR", returnplot=T)
### A big ball of nonsense!

VARHLag <- sparseVAR(Y=ysim) # VARpen="HLag" is the default
par(mfrow=c(1, 1))
LhatHLag <- lagmatrix(fit=VARHLag, model="VAR", returnplot=T)

### Not too bad!!
causality_matrix

#### Let's add a maximum order of 2 ###

VARHLag <- sparseVAR(Y=ysim,p=2) # VARpen="HLag" is the default
par(mfrow=c(1, 1))
LhatHLag <- lagmatrix(fit=VARHLag, model="VAR", returnplot=T)

### And of 1 (not possible to use Hlag here) ####

VARL1 <- sparseVAR(Y=ysim, VARpen="L1",p=1) # default forecast horizon is h=1
par(mfrow=c(1, 1))
LhatL1 <- lagmatrix(fit=VARL1, model="VAR", returnplot=T)
LhatL1
Phi1=as.matrix(LhatL1$LPhi)
sum(Phi1==causality_matrix)/length(as.vector(Phi1))
#0.84 // one of the best scores but not necessarily better than the pairwise

### Use the other data
VARHLag2 <- sparseVAR(Y=as.matrix(yd)) # VARpen="HLag" is the default
par(mfrow=c(1, 1))
LhatHLag <- lagmatrix(fit=VARHLag2, model="VAR", returnplot=T)


VARL1 <- sparseVAR(Y=as.matrix(yd), VARpen="L1",VARgran = c(100,100),p=1) # default forecast horizon is h=1
par(mfrow=c(1, 1))
LhatL1 <- lagmatrix(fit=VARL1, model="VAR", returnplot=T)
#### Yiiik - ugly // absolute nonsense, either raw LASSO or absolute values. 

### Fairly bad for all of these last ones on the conditional GC, whilst working relatively well pairwise with FDR correction. 

### Perhaps only the pairwise GC is robust to changes in nonlinearities? 
### That would make senses to some extent. 

### Damn, looks like the nesting selects the highest lag on the diagonal? 
### Because the model with the nonzero 3rd element is nested within 1,2,3 lags nonzero, 
### but the model with only the first diag element is not compared. 
### Is this correct? If so I should ask, perhaps I need my own stuff... 

################### Some more ideas ##########################

# I could use the pairwise GC to define modules 
# Then fit a VAR on modules

# Or fit a VAR with the support restricted by the pairwise GC. 
# Or fit a VAR with a support defined by some grouping of the pairwise GC structure. 

# Or use the Fisher and Mehta algorithm. 
# Or something fairly similar, e.g. I start with 2x2 models. 
# Then among the ones that are found significant, I go for 3x3
# Then 4x4 [allowing some links to be zero]

# In any case it looks like we need, for any non-pairwise-method, 
# a number m of training sets (to construct LASSO or F&M)

#################### Use LIMITS from pakage seqtime ##############################

#library(devtools)  
#install_github("hallucigenia-sparsa/seqtime")  
library(seqtime)  

# Plotting the B (or equivalent) matrix
network=plotA(B-diag(10),method="network")
taxa_timeseries = t(as.matrix(yd)) # need to tranpose the matrix

Aest=limits(taxa_timeseries)$Aest ###Tryout the LIMITS algo on the Lotka-Volterra model

# True Lotka-Volterra interaction matrix
A= rbind(c(-4,-2,-0.4,0,0,0,0,0,0,0),
         c(-0.31,-3.1,-0.93,0,0,0,0,0,0,0),
         c(0.636,0.636,-2.12,0,0,0,0,0,0,0),
         c(0.111,-0.111,0.131,-3.8,0,0,0,0,0,0),
         c(0,0,0,0.5,-2,-2,-0.4,0,0,0),
         c(0,0,0,0,-0.31,-3.1,-0.93,0,0,0),
         c(0,0,0,0,0.636,0.636,-2.12,0,0,0),
         c(0,0,0,0,0,0,0,-4,-2,-0.4),
         c(0,0,0,0,0,0,0,-0.31,-3.1,-0.93),
        c(0,0,0,0,0,0,0,0.636,0.636,-2.12))
A

par(mfrow=c(1,2))
plotA(A,header="known")
plotA(Aest,header="inferred")

### From seqtime reported statistics (not necessarily very adequate)
crossCor=cor(A,Aest)
mean(diag(crossCor), na.rm=TRUE)





