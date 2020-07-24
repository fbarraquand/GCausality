########################################################################################################################
########### CP 18/04/2019 - Diagnostic of p-values and thresholds for both GC and CCM, can be used for chaotic and stochastic sim. ###########
########### CP 24/05/2019 - Added Matthews correlation, Sokal-Michener and Yule's indices for similarity
########### CP 22/07/2019 - Added numbering to figures and automatic switch from chaos to stochastic data
########### CP 22/07/2019 - Removed Yule's Q index from the final table, too hard to interprete
########################################################################################################################

rm(list=ls())
graphics.off()

#type="CHAOS"
type="stochModel"

library(xtable)

###Define Matthews correlation phi. Not helpful but leaving that here just in case we need that later on
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


tab_inter=read.csv(paste("results/DataCompet",type,"inter_withRhoMaxSpec.csv",sep="_"))
tab_nointer=read.csv(paste("results/DataCompet",type,"noInter_withRhoMaxSpec.csv",sep="_"))

#colo=c("red","blue","orange","cyan")
colo=c("blue","cyan","red","orange")

alpha=0.1

##################################FIGS2#############################
###Work on GC
pdf(paste("fig/explo",type,"GC.pdf",sep="_"),width=10,height=10)
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
plot(seq_test,perc_fn_12[,1],col="blue",ylim=c(0,0.2),lty=1,t="l",xlab="Threshold",ylab="%",main="Performance=f(val)")
lines(seq_test,perc_fn_21[,1],col="red")
lines(seq_test,perc_fn_12[,2],col="blue",lty=2)
lines(seq_test,perc_fn_21[,2],col="red",lty=2)

lines(seq_test,perc_fp_12[,1],col="cyan")
lines(seq_test,perc_fp_21[,1],col="orange")
lines(seq_test,perc_fp_12[,2],col="cyan",lty=2)
lines(seq_test,perc_fp_21[,2],col="orange",lty=2)
abline(v=0.04)
mtext("d)",side=2,las=2,at=0.2*1.1)
legend("top",c("1->2 false neg.","2->1 false neg.","1->2 false pos.","2->1 false pos.","log-ratio","mean effect"),col=c("blue","red","cyan","orange","black","black"),lty=c(1,1,1,1,1,2),bty="n")

dev.off()
#####################################END FIGS2#######################

######FOR table
#table_to_write=matrix(NA,4,10)
table_to_write=matrix(NA,4,9)
rownames(table_to_write)=c("Inter12","Inter21","NoInter12","NoInter21")
#colnames(table_to_write)=c("GCpval","GCLR","GCboth","CCMpval","CCMrho1","CCMrho2","both1","both2","ISM","Q")
colnames(table_to_write)=c("GCpval","GCLR","GCboth","CCMpval","CCMrho1","CCMrho2","both1","both2","ISM")
alpha=0.1
threshold=0.04

table_to_write[1,1]=sum(tab_inter$Pval_12_inter_GC<alpha)/nrow(tab_inter)
table_to_write[1,2]=sum(tab_inter$log_12_inter>threshold)/nrow(tab_inter)
table_to_write[1,3]=sum((tab_inter$Pval_12_inter_GC<alpha)&(tab_inter$log_12_inter>threshold))/nrow(tab_inter)
tab_inter$index_1cause2_inter_GC=(tab_inter$Pval_12_inter_GC<alpha)*(tab_inter$log_12_inter>threshold)

table_to_write[2,1]=sum(tab_inter$Pval_21_inter_GC<alpha)/nrow(tab_inter)
table_to_write[2,2]=sum(tab_inter$log_21_inter>threshold)/nrow(tab_inter)
table_to_write[2,3]=sum((tab_inter$Pval_21_inter_GC<alpha)&(tab_inter$log_21_inter>threshold))/nrow(tab_inter)
tab_inter$index_2cause1_inter_GC=(tab_inter$Pval_21_inter_GC<alpha)*(tab_inter$log_21_inter>threshold)

table_to_write[3,1]=sum(tab_nointer$Pval_12_noInter_GC<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[3,2]=sum(tab_nointer$log_12_noInter>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[3,3]=sum((tab_nointer$Pval_12_noInter_GC<alpha)&(tab_nointer$log_12_noInter>threshold),na.rm=T)/nrow(tab_nointer)
tab_nointer$index_1cause2_inter_GC=(tab_nointer$Pval_12_noInter_GC<alpha)*(tab_nointer$log_12_noInter>threshold)

table_to_write[4,1]=sum(tab_nointer$Pval_21_noInter_GC<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[4,2]=sum(tab_nointer$log_21_noInter>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[4,3]=sum((tab_nointer$Pval_21_noInter_GC<alpha)&(tab_nointer$log_21_noInter>threshold),na.rm=T)/nrow(tab_nointer)
tab_nointer$index_2cause1_inter_GC=(tab_nointer$Pval_21_noInter_GC<alpha)*(tab_nointer$log_21_noInter>threshold)




###Work on CCM
#########################################FIGS3####################################
pdf(paste("fig/explo",type,"CCM_pval_tmp.pdf",sep="_"),width=10,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,2,3,1))
z=data.frame(tab_inter$Pval_12_inter_CCM,tab_nointer$Pval_12_noInter_CCM,tab_inter$Pval_21_inter_CCM,tab_nointer$Pval_21_noInter_CCM)
boxplot(log10(z),col=colo,range=0,main="CobeyBaskerville",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-2.25,0))
mtext("a)",side=2,las=2,at=0.25)
p0=lapply(z,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-2.2))
abline(h=log10(alpha))


z=data.frame(tab_inter$Pval_12_inter_CCM_surr,tab_nointer$Pval_12_noInter_CCM_surr,tab_inter$Pval_21_inter_CCM_surr,tab_nointer$Pval_21_noInter_CCM_surr)
a=boxplot(log10(z),col=colo,range=0,main="Permutation",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-2.25,0))
mtext("b)",side=2,las=2,at=0.25)
p0=lapply(z,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=log10(alpha))

z=data.frame(tab_inter$Pval_12_inter_CCM_surr_twin,tab_nointer$Pval_12_noInter_CCM_surr_twin,tab_inter$Pval_21_inter_CCM_surr_twin,tab_nointer$Pval_21_noInter_CCM_surr_twin)
a=boxplot(log10(z),col=colo,range=0,main="Twin",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-2.25,0))
mtext("c)",side=2,las=2,at=0.25)
p0=lapply(z,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=log10(alpha))

z=data.frame(tab_inter$Pval_12_inter_CCM_surr_ebi,tab_nointer$Pval_12_noInter_CCM_surr_ebi,tab_inter$Pval_21_inter_CCM_surr_ebi,tab_nointer$Pval_21_noInter_CCM_surr_ebi)
a=boxplot(log10(z),col=colo,range=0,main="Ebisuzaki",names=c("1->2 inter","1->2 no inter","2->1 inter","2->1 no inter"),ylim=c(-2.25,0))
abline(h=log10(alpha))
mtext("d)",side=2,las=2,at=0.25)
p0=lapply(z,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
dev.off()
#########################################END FIGS3####################################


pdf(paste("fig/explo",type,"rho_val_CCM.pdf",sep="_"),width=10,height=5)
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


####################### FIGS4 ###############################
pdf(paste("fig/threshold_according_to_pval_and_rho_CCM_",type,".pdf",sep=""),width=10,height=10)
par(mfrow=c(2,2),mar=c(4,4,3,1),lwd=2)
endj=c("","_surr","_surr_twin","_surr_ebi")
mainj=c("CobeyBaskerville","Permutation","Twin","Ebisuzaki")
seq_test=seq(0.01,0.3,0.01)
margin=c("a)","b)","c)","d)")
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
plot(seq_test,perc_fn_12[,1],col="blue",ylim=c(0,yl),lty=1,t="l",xlab="Threshold",ylab="% errors",main=mainj[j])
mtext(margin[j],side=2,las=2,at=yl*1.1)
lines(seq_test,perc_fn_21[,1],col="red")
lines(seq_test,perc_fn_12[,2],col="blue",lty=2)
lines(seq_test,perc_fn_21[,2],col="red",lty=2)

lines(seq_test,perc_fp_12[,1],col="cyan")
lines(seq_test,perc_fp_21[,1],col="orange")
lines(seq_test,perc_fp_12[,2],col="cyan",lty=2)
lines(seq_test,perc_fp_21[,2],col="orange",lty=2)
abline(v=0.1)
 }
legend("topright",c("1->2 false neg.","2->1 false neg.","1->2 false pos.","2->1 false pos.","rho with replacement","rho without replacement"),col=c("blue","red","cyan","orange","black","black"),lty=c(1,1,1,1,1,2),bty="n")

dev.off()
####################################"END FIGS4############################

alpha=0.1
threshold=0.1

table_to_write[1,4]=sum(tab_inter$Pval_12_inter_CCM_surr<alpha)/nrow(tab_inter)
table_to_write[1,5]=sum(tab_inter$RhoLMax_12_inter_v2>0.1)/nrow(tab_inter)
table_to_write[1,6]=sum(tab_inter$RhoLMax_12_inter_v2>0.2)/nrow(tab_inter)
table_to_write[1,7]=sum((tab_inter$Pval_12_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_12_inter_v2>0.1))/nrow(tab_inter)
table_to_write[1,8]=sum((tab_inter$Pval_12_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_12_inter_v2>0.2))/nrow(tab_inter)
tab_inter$index_1cause2_inter_CCM=(tab_inter$Pval_12_inter_CCM_surr<alpha)*(tab_inter$RhoLMax_12_inter_v2>threshold)

table_to_write[2,4]=sum(tab_inter$Pval_21_inter_CCM_surr<alpha)/nrow(tab_inter)
table_to_write[2,5]=sum(tab_inter$RhoLMax_21_inter_v2>0.1)/nrow(tab_inter)
table_to_write[2,6]=sum(tab_inter$RhoLMax_21_inter_v2>0.2)/nrow(tab_inter)
table_to_write[2,7]=sum((tab_inter$Pval_21_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_21_inter_v2>0.1))/nrow(tab_inter)
table_to_write[2,8]=sum((tab_inter$Pval_21_inter_CCM_surr<alpha)&(tab_inter$RhoLMax_21_inter_v2>0.2))/nrow(tab_inter)
tab_inter$index_2cause1_inter_CCM=(tab_inter$Pval_21_inter_CCM_surr<alpha)*(tab_inter$RhoLMax_21_inter_v2>threshold)

table_to_write[3,4]=sum(tab_nointer$Pval_12_noInter_CCM_surr<alpha)/nrow(tab_nointer)
table_to_write[3,5]=sum(tab_nointer$RhoLMax_12_noInter_v2>0.1)/nrow(tab_nointer)
table_to_write[3,6]=sum(tab_nointer$RhoLMax_12_noInter_v2>0.2)/nrow(tab_nointer)
table_to_write[3,7]=sum((tab_nointer$Pval_12_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>0.1))/nrow(tab_nointer)
table_to_write[3,8]=sum((tab_nointer$Pval_12_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_12_noInter_v2>0.2))/nrow(tab_nointer)
tab_nointer$index_1cause2_inter_CCM=(tab_nointer$Pval_12_noInter_CCM_surr<alpha)*(tab_nointer$RhoLMax_12_noInter_v2>threshold)

table_to_write[4,4]=sum(tab_nointer$Pval_21_noInter_CCM_surr<alpha)/nrow(tab_nointer)
table_to_write[4,5]=sum(tab_nointer$RhoLMax_21_noInter_v2>0.1)/nrow(tab_nointer)
table_to_write[4,6]=sum(tab_nointer$RhoLMax_21_noInter_v2>0.2)/nrow(tab_nointer)
table_to_write[4,7]=sum((tab_nointer$Pval_21_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_21_noInter_v2>0.1))/nrow(tab_nointer)
table_to_write[4,8]=sum((tab_nointer$Pval_21_noInter_CCM_surr<alpha)&(tab_nointer$RhoLMax_21_noInter_v2>0.2))/nrow(tab_nointer)
tab_nointer$index_2cause1_inter_CCM=(tab_nointer$Pval_21_noInter_CCM_surr<alpha)*(tab_nointer$RhoLMax_21_noInter_v2>threshold)


### For phi
plou=table(tab_inter$index_1cause2_inter_GC,tab_inter$index_1cause2_inter_CCM)
table_to_write[1,9]=sk_index(plou)
#table_to_write[1,10]=yule_index(plou)
plou=table(tab_inter$index_2cause1_inter_GC,tab_inter$index_2cause1_inter_CCM)
table_to_write[2,9]=sk_index(plou)
#table_to_write[2,10]=yule_index(plou)
plou=table(tab_nointer$index_1cause2_inter_GC,tab_nointer$index_1cause2_inter_CCM)
table_to_write[3,9]=sk_index(plou)
#table_to_write[3,10]=yule_index(plou)
plou=table(tab_nointer$index_2cause1_inter_GC,tab_nointer$index_2cause1_inter_CCM)
table_to_write[4,9]=sk_index(plou)
#table_to_write[4,10]=yule_index(plou)

table_to_write[,1:8]=100*table_to_write[,1:8]

print.xtable(xtable(table_to_write,digits=c(1,rep(1,8),2)),paste("results/pval_threshold_",type,".tex",sep=""),type="latex")
