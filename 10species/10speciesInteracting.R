########################################################################################################
### FBarraquand 06/10/2017 -- from previous code in 2017 // adding corrections for multiple testing
########################################################################################################

library("vars")

set.seed(40) ### results may differ based on the simulation, because some species get really variable dynamics

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

## VAR(y=yd, type="none",lag.max=15) ## stupid choice for a full VAR
VAR(y=yd, type="none",lag.max=5,ic="SC")

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
  
# Better to either compute directly pairwise GC or conditional GC (based on VAR coefs)
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
  coefi=coef(var1)[[i]]
  for (j in 1:10){
    p_value[i,j]=coefi[j,4]
    B[i,j]=coefi[j,1]
  }
}
B
round(p_value, digits=2)
B
### Printing B
library(xtable)
Btable <- xtable(B)
print(Btable, floating=FALSE, tabular.environment="pmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

### Using the p-values on the coefficients. 
p_fullMat_adjusted=p.adjust(p_value,method="BH")
p_rounded=round(matrix(p_fullMat_adjusted,nrow=10,ncol=10), digits=4)
p_rounded<0.05
causality_matrix ## quite a number of errors... 

G=(p_rounded<0.05)
G
### Percentage of correct interactions
(sum(G==causality_matrix)-nspecies)/(length(as.vector(G))-nspecies)
# 69%, 77% previously - not so bad

###############################################################################
######################## Try pairwise GC now ##################################
###############################################################################

p_value1=matrix(0,nrow=nspecies,ncol=nspecies)
for (i in 1:nspecies){
  for (j in 1:nspecies){
    # cause first and effet later in grangertest()
    if (i !=j){
      g=grangertest(y[,j],y[,i],order=1)
      p_value1[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
    }
  }
}

G_1=(p_value1<0.05)
G_1
### Percentage of correct interactions
(sum(G_1==causality_matrix)-nspecies)/(length(as.vector(G_1))-nspecies) #0.76
## Surprisingly this was 0.89 with non log-transformed Y
G_1=(p_value1<0.01)
G_1
(sum(G_1==causality_matrix)-nspecies)/(length(as.vector(G_1))-nspecies) #0.78 #0.81 previously - misleading because we find mostly no interactions...
## I should use a ROC curve? Or simply differentiate between interactions that are there and missed and vice versa.

### Better adjustment
p_ajusted_pairwise=matrix(p.adjust(p_value1,method="BH"),nrow=nspecies,ncol=nspecies)
p_adj_round=round(p_ajusted_pairwise,digits=4)

G_adj=(p_adj_round<0.05)
G_adj
### Percentage of correct interactions
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.76 #0.81 previously
G_adj=(p_value1<0.01)
G_adj
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.78 #0.81 previously
### 80% of correct interactions (compared to 90% earlier with Y on some other sims)
### I would need to do this over many repeats

### That's good / though we miss interactions rather than we infer some 
### We should really CCM this to see what we get!
### (I doubt that we get that high...)

### Limitation though: here we assumed the MAR1 model fits for each pairs... 
### What may be optimal for *all* species may not be optimal for *each* species

#####################################################################################################
####### Let's try a more complex model with potentially variable lags, with BIC selection. Lagmax=5
#####################################################################################################

p_value_wald=p_value_F=estimated_modelOrder=matrix(0,nrow=nspecies,ncol=nspecies) ### Initializing arrays for P-values
for (i in 1:nspecies){
  for (j in 1:nspecies){
    # cause first and effet later in grangertest()
    if (i !=j){
      datapair=data.frame(cbind(y[,i],y[,j]))
      names(datapair)=c("yi","yj")
      varmodel=VAR(y=datapair, type="none",ic="SC",lag.max=5)
      p_order=varmodel$p
      estimated_modelOrder[i,j]=p_order
      cause=causality(varmodel,cause="yj")
      p_value_F[i,j]=cause$Granger$p.value 
      #p-value F-test now
      ## Now use a Wald test
      g=grangertest(y[,j],y[,i],order=p_order) ## y[,i] is the explained variable
      p_value_wald[i,j] =  as.numeric(sprintf("%.5f",g$`Pr(>F)`[2]))
      
    }
  }
}

estimated_modelOrder ### To check that it does indeed change smthg

G_wald=(p_value_wald<0.05)
G_wald

G_F=(p_value_F<0.05)
G_F

### Percentage of correct interactions
(sum(G_wald==causality_matrix)-nspecies)/(length(as.vector(G_wald))-nspecies) #0.78, no change compared to order=1

G_wald=(p_value_wald<0.01) ## What about a 0.01 threshold
G_wald
(sum(G_wald==causality_matrix)-nspecies)/(length(as.vector(G_wald))-nspecies) #0.78

### Percentage of correct interactions
(sum(G_F==causality_matrix)-nspecies)/(length(as.vector(G_F))-nspecies) #0.71

G_F=(p_value_F<0.01) ## What about a 0.01 threshold
G_F
(sum(G_F==causality_matrix)-nspecies)/(length(as.vector(G_F))-nspecies) #0.76

### check both p-values
cor(as.vector(p_value_F),as.numeric(as.vector(p_value_wald)))
plot(as.vector(p_value_F),as.numeric(as.vector(p_value_wald))) ## OK, this is a bit worrying... 
### The F-test finds more interactions. While BH rarely changes the Wald p-values, perhaps it will change those


### Better adjustment
p_ajusted_pairwise=matrix(p.adjust(p_value_F,method="BH"),nrow=nspecies,ncol=nspecies)
p_adj_round=round(p_ajusted_pairwise,digits=4)

G_adj=(p_adj_round<0.05)
G_adj
### Percentage of correct interactions
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.74 not really
G_adj=(p_value1<0.01)
G_adj
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.78 

### NB what about 
G_adj=(p_adj_round<0.1)
G_adj
### Percentage of correct interactions
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.71

### Do the same with Wald p-value
### Better adjustment
p_ajusted_pairwise=matrix(p.adjust(p_value_wald,method="BH"),nrow=nspecies,ncol=nspecies)
p_adj_round=round(p_ajusted_pairwise,digits=4)

G_adj=(p_adj_round<0.05)
G_adj
### Percentage of correct interactions
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.77 
G_adj=(p_value1<0.01)
G_adj
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.78 

### NB what about 
G_adj=(p_adj_round<0.1)
G_adj
### Percentage of correct interactions
(sum(G_adj==causality_matrix)-nspecies)/(length(as.vector(G_adj))-nspecies) #0.78

### The Wald test seems in fact more reliable... 