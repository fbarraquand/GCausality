########################################################################################################################
########### CP 18/04/2019 - Diagnostic of p-values and thresholds for both GC and CCM, can be used for chaotic and stochastic sim. ###########
########################################################################################################################

rm(list=ls())
graphics.off()

###Define Matthews correlation phi. Not helpful but living that here just in case we need that later on
phi=function(tableau){

        if(("1" %in% rownames(tableau))&("1" %in% colnames(tableau))){
                n11=tableau["1","1"]
        }else{
                n11=0
        }
        if(("1" %in% rownames(tableau))&("0" %in% colnames(tableau))){
                n10=tableau["1","0"]
        }else{
                n10=0
        }
        if(("0" %in% rownames(tableau))&("1" %in% colnames(tableau))){
                n01=tableau["0","1"]
        }else{
                n01=0
        }
        if(("0" %in% rownames(tableau))&("0" %in% colnames(tableau))){
                n00=tableau["0","0"]
        }else{
                n00=0
        }

        n1p=n11+n10
        n0p=n01+n00
        np1=n11+n01
        np0=n00+n10

        n=n11+n10+n01+n00

        den=1.0*n1p*n0p*np0*n1p
        if(den==0.0){den=1.0}

        phi=(n11*n00-n10*n01)/sqrt(den) #the 1.0 to avoid integer overflow

        if(is.na(phi)){
                phi=(n*n11-n1p*np1)/sqrt(n1p*np1*(n-n1p)*(n-np1))
        }

        return(phi)
}

### Define Sokal Michener
sk_index=function(tableau){
        if(("1" %in% rownames(tableau))&("1" %in% colnames(tableau))){
                n11=tableau["1","1"]
        }else{
                n11=0
        }
        if(("0" %in% rownames(tableau))&("0" %in% colnames(tableau))){
                n00=tableau["0","0"]
        }else{
                n00=0
        }
	n=sum(tableau)

	ind=(n11+n00)/n

	return(ind)
}

### Define Yule's Q
yule_index=function(tableau){


        if(("1" %in% rownames(tableau))&("1" %in% colnames(tableau))){
                n11=tableau["1","1"]
        }else{
                n11=0
        }
        if(("1" %in% rownames(tableau))&("0" %in% colnames(tableau))){
                n10=tableau["1","0"]
        }else{
                n10=0
        }
        if(("0" %in% rownames(tableau))&("1" %in% colnames(tableau))){
                n01=tableau["0","1"]
        }else{
                n01=0
        }
        if(("0" %in% rownames(tableau))&("0" %in% colnames(tableau))){
                n00=tableau["0","0"]
        }else{
                n00=0
        }

	if(n11==0|n00==0){
		id=1
	}else{
		id=(n11*n00-n10*n01)/(n11*n00+n10*n01)
	}

	return(id)
}


#tab_inter=read.csv("results/DataCompet_stochModel_inter_withRhoMaxSpec.csv")
#tab_nointer=read.csv("results/DataCompet_stochModel_noInter_withRhoMaxSpec.csv")
tab_inter=read.csv("results/DataCompet_CHAOS_inter_withRhoMaxSpec.csv")
tab_nointer=read.csv("results/DataCompet_CHAOS_noInter_withRhoMaxSpec.csv")

colo=c("red","blue","orange","cyan")

alpha=0.1

###Work on GC
pdf("explo_CHAOS_GC.pdf",width=10,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,2,3,1))

logz=data.frame(tab_inter$log_12_inter,tab_nointer$log_12_noInter,tab_inter$log_21_inter,tab_nointer$log_21_noInter)
boxplot(logz,col=colo,range=0,main="log ratio",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"))
mtext("a)",side=2,las=2,at=max(logz)*1.1)

lagz=data.frame(tab_inter$lag_order_inter_GC,tab_nointer$lag_order_noInter_GC)

effectz=data.frame(tab_inter$effect_12_inter/tab_inter$lag_order_inter_GC,tab_nointer$effect_12_noInter/tab_nointer$lag_order_noInter_GC,tab_inter$effect_21_inter/tab_inter$lag_order_inter_GC,tab_nointer$effect_21_noInter/tab_nointer$lag_order_noInter_GC)
boxplot(effectz,col=colo,range=0,main="mean effect",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"))
mtext("b)",side=2,las=2,at=max(effectz)*1.1)

z=data.frame(tab_inter$Pval_12_inter_GC,tab_nointer$Pval_12_noInter_GC,tab_inter$Pval_21_inter_GC,tab_nointer$Pval_21_noInter_GC)
boxplot(log10(z),col=colo,range=0,main="log10(p-value)",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-5,0))
mtext("c)",side=2,las=2,at=5*0.1)
abline(h=log10(alpha))

#Plot false negatives and false positives as a function of the threshold
seq_test=seq(0.01,0.3,0.01)
perc_fn_12=matrix(NA,length(seq_test),2) #1=log, 2=effect
perc_fp_12=matrix(NA,length(seq_test),2)
perc_fn_21=matrix(NA,length(seq_test),2) #1=log, 2=effect
perc_fp_21=matrix(NA,length(seq_test),2)
for (i in 1:length(seq_test)){
	perc_fn_12[i,1]=sum((tab_inter$Pval_12_inter_GC>alpha)|(tab_inter$log_12_inter<seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fn_21[i,1]=sum((tab_inter$Pval_21_inter_GC>alpha)|(tab_inter$log_21_inter<seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fp_12[i,1]=sum((tab_nointer$Pval_12_noInter_GC<alpha)&(tab_nointer$log_12_noInter>seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fp_21[i,1]=sum((tab_nointer$Pval_21_noInter_GC<alpha)&(tab_nointer$log_21_noInter>seq_test[i]),na.rm=T)/nrow(tab_inter)
	
	perc_fn_12[i,2]=sum((tab_inter$Pval_12_inter_GC>alpha)|(tab_inter$effect_12_inter<seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fn_21[i,2]=sum((tab_inter$Pval_21_inter_GC>alpha)|(tab_inter$effect_21_inter<seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fp_12[i,2]=sum((tab_nointer$Pval_12_noInter_GC<alpha)&(tab_nointer$effect_12_noInter>seq_test[i]),na.rm=T)/nrow(tab_inter)
	perc_fp_21[i,2]=sum((tab_nointer$Pval_21_noInter_GC<alpha)&(tab_nointer$effect_21_noInter>seq_test[i]),na.rm=T)/nrow(tab_inter)
}
par(lwd=1.5)
plot(seq_test,perc_fn_12[,1],col="red",ylim=c(0,0.2),lty=1,t="l",xlab="Threshold",ylab="%",main="Performance=f(val)")
lines(seq_test,perc_fn_21[,1],col="orange")
lines(seq_test,perc_fn_12[,2],col="red",lty=2)
lines(seq_test,perc_fn_21[,2],col="orange",lty=2)

lines(seq_test,perc_fp_12[,1],col="blue")
lines(seq_test,perc_fp_21[,1],col="cyan")
lines(seq_test,perc_fp_12[,2],col="blue",lty=2)
lines(seq_test,perc_fp_21[,2],col="cyan",lty=2)
abline(v=0.04)
mtext("d)",side=2,las=2,at=0.2*1.1)
legend("top",c("1->2 false neg.","2->1 false neg.","1->2 false pos.","2->1 false pos.","log-ratio","mean effect"),col=c("red","orange","blue","cyan","black","black"),lty=c(1,1,1,1,1,2),bty="n")

dev.off()

######FOR table
print("GC BASED")
alpha=0.1
threshold=0.04
print("1 causes 2, with")
print(sum(tab_inter$Pval_12_inter_GC<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_12_inter>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_GC<alpha)&(tab_inter$log_12_inter>threshold))/nrow(tab_inter))
tab_inter$index_1cause2_inter_GC=(tab_inter$Pval_12_inter_GC<alpha)*(tab_inter$log_12_inter>threshold)
print("2 causes 1, with")
print(sum(tab_inter$Pval_21_inter_GC<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_21_inter>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_GC<alpha)&(tab_inter$log_21_inter>threshold))/nrow(tab_inter))
tab_inter$index_2cause1_inter_GC=(tab_inter$Pval_21_inter_GC<alpha)*(tab_inter$log_21_inter>threshold)
print("1 causes 2, without")
print(sum(tab_nointer$Pval_12_noInter_GC<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_12_noInter>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_noInter_GC<alpha)&(tab_nointer$log_12_noInter>threshold),na.rm=T)/nrow(tab_nointer))
tab_nointer$index_1cause2_inter_GC=(tab_nointer$Pval_12_noInter_GC<alpha)*(tab_nointer$log_12_noInter>threshold)
print("2 causes 1, without")
print(sum(tab_nointer$Pval_21_noInter_GC<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_21_noInter>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_noInter_GC<alpha)&(tab_nointer$log_21_noInter>threshold),na.rm=T)/nrow(tab_nointer))
tab_nointer$index_2cause1_inter_GC=(tab_nointer$Pval_21_noInter_GC<alpha)*(tab_nointer$log_21_noInter>threshold)




###Work on CCM

pdf("explo_CHAOS_CCM_pval_tmp.pdf",width=10,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,2,3,1))
z=data.frame(tab_inter$Pval_12_inter_CCM,tab_nointer$Pval_12_noInter_CCM,tab_inter$Pval_21_inter_CCM,tab_nointer$Pval_21_noInter_CCM)
boxplot(log10(z),col=colo,range=0,main="CobeyBaskerville",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-3.5,0))
p0=lapply(z,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=log10(alpha))


z=data.frame(tab_inter$Pval_12_inter_CCM_surr,tab_nointer$Pval_12_noInter_CCM_surr,tab_inter$Pval_21_inter_CCM_surr,tab_nointer$Pval_21_noInter_CCM_surr)
boxplot(log10(z),col=colo,range=0,main="Permutation",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-3.5,0))
p0=lapply(z,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=log10(alpha))

z=data.frame(tab_inter$Pval_12_inter_CCM_surr_twin,tab_nointer$Pval_12_noInter_CCM_surr_twin,tab_inter$Pval_21_inter_CCM_surr_twin,tab_nointer$Pval_21_noInter_CCM_surr_twin)
boxplot(log10(z),col=colo,range=0,main="Twin",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-3.5,0))
p0=lapply(z,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=log10(alpha))

z=data.frame(tab_inter$Pval_12_inter_CCM_surr_ebi,tab_nointer$Pval_12_noInter_CCM_surr_ebi,tab_inter$Pval_21_inter_CCM_surr_ebi,tab_nointer$Pval_21_noInter_CCM_surr_ebi)
boxplot(log10(z),col=colo,range=0,main="Ebisuzaki",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-3.5,0))
abline(h=log10(alpha))
p0=lapply(z,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
dev.off()

pdf("explo_CHAOS_rho_val_CCM.pdf",width=10,height=5)
par(mfrow=c(1,3))
rhoz1=data.frame(tab_inter$RhoLMax_12_inter_v1,tab_nointer$RhoLMax_12_noInter_v1,tab_inter$RhoLMax_21_inter_v1,tab_nointer$RhoLMax_12_noInter_v1)
boxplot(rhoz1,col=colo,range=0,main="Rho max bootstrap",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"))
rhoz2=data.frame(tab_inter$RhoLMax_12_inter_v2,tab_nointer$RhoLMax_12_noInter_v2,tab_inter$RhoLMax_21_inter_v2,tab_nointer$RhoLMax_12_noInter_v2)
boxplot(rhoz2,col=colo,range=0,main="Rho max real data",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"))
plot(0,0,t="n",xlim=c(0,1.0),ylim=c(0,1.0),xlab="bootstrap",ylab="real")
for(i in 1:length(rhoz1)){
	points(rhoz1[[i]],rhoz2[[i]],col=colo[i],pch=16)
}
abline(a=0,b=1)
dev.off()

pdf("threshold_according_to_pval_and_rho_CCM_CHAOS.pdf",width=10,height=10)
par(mfrow=c(2,2),mar=c(4,2,3,1),lwd=2)
endj=c("","_surr","_surr_twin","_surr_ebi")
mainj=c("CobeyBaskerville","Permutation","Twin","Ebisuzaki")
seq_test=seq(0.01,0.3,0.01)
#Plot false negatives and false positives as a function of the threshold
for(j in 1:length(endj)){
perc_fn_12=matrix(NA,length(seq_test),2) #1=rhov1,rhov2
perc_fp_12=matrix(NA,length(seq_test),2)
for (i in 1:length(seq_test)){
        perc_fn_12[i,1]=sum((tab_inter[,paste("Pval_12_inter_CCM",endj[j],sep="")]>alpha)|(tab_inter$RhoLMax_12_inter_v1<seq_test[i]))/nrow(tab_inter)
        perc_fn_21[i,1]=sum((tab_inter[,paste("Pval_21_inter_CCM",endj[j],sep="")]>alpha)|(tab_inter$RhoLMax_21_inter_v1<seq_test[i]))/nrow(tab_inter)
        perc_fp_12[i,1]=sum((tab_nointer[,paste("Pval_12_noInter_CCM",endj[j],sep="")]<alpha)&(tab_nointer$RhoLMax_12_noInter_v1>seq_test[i]))/nrow(tab_inter)
        perc_fp_21[i,1]=sum((tab_nointer[,paste("Pval_21_noInter_CCM",endj[j],sep="")]<alpha)&(tab_nointer$RhoLMax_12_noInter_v1>seq_test[i]))/nrow(tab_inter)

        perc_fn_12[i,2]=sum((tab_inter[,paste("Pval_12_inter_CCM",endj[j],sep="")]>alpha)|(tab_inter$RhoLMax_12_inter_v2<seq_test[i]))/nrow(tab_inter)
        perc_fn_21[i,2]=sum((tab_inter[,paste("Pval_21_inter_CCM",endj[j],sep="")]>alpha)|(tab_inter$RhoLMax_21_inter_v2<seq_test[i]))/nrow(tab_inter)
        perc_fp_12[i,2]=sum((tab_nointer[,paste("Pval_12_noInter_CCM",endj[j],sep="")]<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>seq_test[i]))/nrow(tab_inter)
        perc_fp_21[i,2]=sum((tab_nointer[,paste("Pval_21_noInter_CCM",endj[j],sep="")]<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>seq_test[i]))/nrow(tab_inter)
}
if(j==1){
yl=0.25
}else{
yl=0.15
}
plot(seq_test,perc_fn_12[,1],col="red",ylim=c(0,yl),lty=1,t="l",xlab="Threshold",ylab="%",main=mainj[j])
lines(seq_test,perc_fn_21[,1],col="orange")
lines(seq_test,perc_fn_12[,2],col="red",lty=2)
lines(seq_test,perc_fn_21[,2],col="orange",lty=2)

lines(seq_test,perc_fp_12[,1],col="blue")
lines(seq_test,perc_fp_21[,1],col="cyan")
lines(seq_test,perc_fp_12[,2],col="blue",lty=2)
lines(seq_test,perc_fp_21[,2],col="cyan",lty=2)
abline(v=0.1)
 }
legend("topright",c("1->2 false neg.","2->1 false neg.","1->2 false pos.","2->1 false pos.","rho with replacement","rho without replacement"),col=c("red","orange","blue","cyan","black","black"),lty=c(1,1,1,1,1,2),bty="n")

dev.off()

print("FOR CCM table")
alpha=0.1
print("1 causes 2, with")
print(sum(tab_inter$Pval_12_inter_CCM_surr<alpha)/nrow(tab_inter))
print(sum(tab_inter$RhoLMax_12_inter_v2>0.1)/nrow(tab_inter))
print(sum(tab_inter$RhoLMax_12_inter_v2>0.2)/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_12_inter_v2>0.1))/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_12_inter_v2>0.2))/nrow(tab_inter))
tab_inter$index_1cause2_inter_CCM=(tab_inter$Pval_12_inter_CCM_surr<alpha)*(tab_inter$RhoLMax_12_inter_v2>0.2)
print("2 causes 1, with")
print(sum(tab_inter$Pval_21_inter_CCM_surr<alpha)/nrow(tab_inter))
print(sum(tab_inter$RhoLMax_21_inter_v2>0.1)/nrow(tab_inter))
print(sum(tab_inter$RhoLMax_21_inter_v2>0.2)/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_21_inter_v2>0.1))/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_21_inter_v2>0.2))/nrow(tab_inter))
tab_inter$index_2cause1_inter_CCM=(tab_inter$Pval_21_inter_CCM_surr<alpha)*(tab_inter$RhoLMax_21_inter_v2>0.2)
print("1 causes 2, without")
print(sum(tab_nointer$Pval_12_noInter_CCM_surr<alpha)/nrow(tab_nointer))
print(sum(tab_nointer$RhoLMax_12_noInter_v2>0.1)/nrow(tab_nointer))
print(sum(tab_nointer$RhoLMax_12_noInter_v2>0.2)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>0.1))/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>0.2))/nrow(tab_nointer))
tab_nointer$index_1cause2_inter_CCM=(tab_nointer$Pval_12_noInter_CCM_surr<alpha)*(tab_nointer$RhoLMax_12_noInter_v2>0.2)
print("2 causes 1, without")
print(sum(tab_nointer$Pval_21_noInter_CCM_surr<alpha)/nrow(tab_nointer))
print(sum(tab_nointer$RhoLMax_21_noInter_v2>0.1)/nrow(tab_nointer))
print(sum(tab_nointer$RhoLMax_21_noInter_v2>0.2)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_21_noInter_v2>0.1))/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_21_noInter_v2>0.2))/nrow(tab_nointer))
tab_nointer$index_2cause1_inter_CCM=(tab_nointer$Pval_21_noInter_CCM_surr<alpha)*(tab_nointer$RhoLMax_21_noInter_v2>0.2)



### For phi
print("######################################### PHI ################################")
print("1 causes 2, with")
plou=table(tab_inter$index_1cause2_inter_GC,tab_inter$index_1cause2_inter_CCM)
#print(phi(plou))
print(yule_index(plou))
print(sk_index(plou))
print("2 causes 1, with")
plou=table(tab_inter$index_2cause1_inter_GC,tab_inter$index_2cause1_inter_CCM)
#print(phi(plou))
print(sk_index(plou))
print(yule_index(plou))
print("1 causes 2, without")
plou=table(tab_nointer$index_1cause2_inter_GC,tab_nointer$index_1cause2_inter_CCM)
#print(phi(plou))
print(sk_index(plou))
print(yule_index(plou))
print("2 causes 1, without")
plou=table(tab_nointer$index_2cause1_inter_GC,tab_nointer$index_2cause1_inter_CCM)
#print(phi(plou))
print(sk_index(plou))
print(yule_index(plou))
