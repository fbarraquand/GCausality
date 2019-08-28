#### FB 25/05/2018, from code Causality_Marseille.R dated oct. 24 2016
#### 5-species competition, almost the model of Sugihara (Ricker used instead of logistic, otherwise similar)
#### Useful to compare pairwise GC (+FDR correction?) to CCM performance with the same models as Sugihara et al. 2012
#### Unclear how to test properly for the statistical significance of the causal relationship with CCM in a 5x5 species context

tmax=300
Y=matrix(1,nrow=tmax,ncol=5)
Y[1,1]=abs(rnorm(1,1,1))
Y[1,2]=abs(rnorm(1,1,1))
Y[1,3]=abs(rnorm(1,1,1))
Y[1,4]=abs(rnorm(1,1,1))
Y[1,5]=abs(rnorm(1,1,1))

for (t in 1:(tmax-1)){
  Y[t+1,1] = Y[t,1]*exp(4-4*Y[t,1]-2*Y[t,2]-0.4*Y[t,3])
  Y[t+1,2] = Y[t,2]*exp(3.1-0.31*Y[t,1]-3.1*Y[t,2]-0.93*Y[t,3])
  Y[t+1,3] = Y[t,3]*exp(0.12 + 0.636*Y[t,1]+0.636*Y[t,2]-2.12*Y[t,3])
  Y[t+1,4] = Y[t,4]*exp(3.8-0.111*Y[t,1]-0.111*Y[t,2]+0.131*Y[t,3]-3.8*Y[t,4])#2.12 ->0.12
  Y[t+1,5] = Y[t,5]*exp(4.1-0.082*Y[t,1]-0.111*Y[t,2]-0.125*Y[t,3]-4.1*Y[t,5])
  
}

y=log(Y)
y=y-mean(y)

vec_col=c("black","blue","green","red","pink")
pdf(file="CompetitionDeterNL5Species.pdf",width=12,height=6)
par(cex=1.5)
miny=min(y[151:201,])
maxy=max(y[151:201,])
for (i in 1:5){
  if (i==1){plot(151:201,y[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(miny,maxy),xlab="Time",ylab="ln(Density)")}
  lines(151:201,y[151:201,i],type="o",col=vec_col[i],lwd=3,ylim=c(miny,maxy))
}
dev.off()

### Now we should recover relations between 
# and Y1  Y2  Y3  Y4  Y5
# Y1  x   x   x
# Y2  x   x   x
# Y3  x   x   x
# Y4  x   x   x   x 
# Y5  x   x   x       x
##### I think we need conditional-G for this, but let's see what we find with pairwise G. 
##### After all we can put ourselves in the nice conditions suggested by Sugihara where we don't know anything
##### Let's consider 1 then 5 lags, then 10, finally 15

p_value1=matrix(0,nrow=5,ncol=5)
p_value5=matrix(0,nrow=5,ncol=5)
p_value10=matrix(0,nrow=5,ncol=5)
p_value15=matrix(0,nrow=5,ncol=5)
p_value7=matrix(0,nrow=5,ncol=5)

for (i in 1:5){
  for (j in 1:5){
    # cause first and effet later in grangertest()
    if (i !=j){
      g=grangertest(Y[,j],Y[,i],order=1)
      p_value1[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
      g=grangertest(Y[,j],Y[,i],order=5)
      p_value5[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
      g=grangertest(Y[,j],Y[,i],order=10)
      p_value10[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
      g=grangertest(Y[,j],Y[,i],order=15)
      p_value15[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
      g=grangertest(Y[,j],Y[,i],order=7)
      p_value7[i,j] =  sprintf("%.5f",g$`Pr(>F)`[2])
    }
  }
}

p_value1
p_value5
p_value10
p_value15
p_value7

G_1=(p_value5<0.05)
G_5=(p_value5<0.05)
G_10=(p_value10<0.05) 
G_15=(p_value15<0.05) 
G_7=(p_value7<0.1) # In that case works approx OK, but some tweaking is done...

### The most amazing thing is that we're not even looking at conditional GC there, 
### still the method gets most of the strongly interacting species.  

causality_matrix = rbind(c(1,1,1,0,0),c(1,1,1,0,0),c(1,1,1,0,0),c(1,1,1,1,0),c(1,1,1,0,1))
causality_matrix # true causality matrix here. 

### Try to find the right number of lags?
y5sp=data.frame(y)
var5sp<-VAR(y=y5sp, type="none",lag.max=15,ic="FPE")
summary(var5sp)#7 lags

var5sp<-VAR(y=y5sp, type="none",lag.max=15,ic="AIC")
summary(var5sp)#7 lags

#### Redo analysis with 7 lags (see above). 

