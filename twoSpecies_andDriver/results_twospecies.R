########################################################################################################################
########### CP 19/04/2019 - Diagnostic of p-values and thresholds for both GC and CCM, can be used for sim with driver. ###########
########################################################################################################################

rm(list=ls())
graphics.off()

colo=c("red","blue","orange","cyan")

pdf("explo_with_driver_GC.pdf",width=10,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,2,3,1))
#let start by GC
tab_GC=read.csv('DataCompet_driver_inter_factorized_GC_otf.csv')
tab_inter=tab_GC[1:500,]
tab_nointer=tab_GC[501:1000,]
logz=data.frame(tab_inter$log_12_inter_exo,tab_inter$log_12_inter_noexo,tab_nointer$log_12_inter_exo,tab_nointer$log_12_inter_noexo)
boxplot(logz,col=colo,range=0,main="log ratio 1->2",names=c("inter conditional","inter pairwise","no inter conditional","no inter pairwise"),ylim=c(0,0.25))
abline(h=0.04)
#mtext("a)",side=2,las=2,at=max(logz)*1.1)

logz=data.frame(tab_inter$log_21_inter_exo,tab_inter$log_21_inter_noexo,tab_nointer$log_21_inter_exo,tab_nointer$log_21_inter_noexo)
boxplot(logz,col=colo,range=0,main="log ratio 2->1",names=c("inter conditional","inter pairwise","no inter conditional","no inter pairwise"),ylim=c(0,0.25))
abline(h=0.04)

pvalz=data.frame(tab_inter$Pval_12_inter_GC_exo,tab_inter$Pval_12_inter_GC_noexo,tab_nointer$Pval_12_inter_GC_exo,tab_nointer$Pval_12_inter_GC_noexo)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal ratio 1->2",names=c("inter conditional","inter pairwise","no inter conditional","no inter pairwise"))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_21_inter_GC_exo,tab_inter$Pval_21_inter_GC_no_exo,tab_nointer$Pval_21_inter_GC_exo,tab_nointer$Pval_21_inter_GC_no_exo)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal ratio 2->1",names=c("inter conditional","inter pairwise","no inter conditional","no inter pairwise"))
abline(h=-1)
dev.off()

print("GC pairwise")
alpha=0.1
threshold=0.04
print("1 causes 2, with")
print(sum(tab_inter$Pval_12_inter_GC_noexo<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_12_inter_noexo>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_GC_noexo<alpha)&(tab_inter$log_12_inter_noexo>threshold))/nrow(tab_inter))
print("2 causes 1, with")
print(sum(tab_inter$Pval_21_inter_GC_no_exo<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_21_inter_noexo>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_GC_no_exo<alpha)&(tab_inter$log_21_inter_noexo>threshold))/nrow(tab_inter))
print("1 causes 2, without")
print(sum(tab_nointer$Pval_12_inter_GC_noexo<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_12_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_inter_GC_noexo<alpha)&(tab_nointer$log_12_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer))
print("1 causes 2, without")
print(sum(tab_nointer$Pval_21_inter_GC_no_exo<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_21_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_inter_GC_no_exo<alpha)&(tab_nointer$log_21_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer))

print("GC conditional")
alpha=0.1
threshold=0.04
print("1 causes 2, with")
print(sum(tab_inter$Pval_12_inter_GC_exo<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_12_inter_exo>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_GC_exo<alpha)&(tab_inter$log_12_inter_exo>threshold))/nrow(tab_inter))
print("2 causes 1, with")
print(sum(tab_inter$Pval_21_inter_GC_exo<alpha)/nrow(tab_inter))
print(sum(tab_inter$log_21_inter_exo>threshold)/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_GC_exo<alpha)&(tab_inter$log_21_inter_exo>threshold))/nrow(tab_inter))
print("1 causes 2, without")
print(sum(tab_nointer$Pval_12_inter_GC_exo<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_12_inter_exo>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_inter_GC_exo<alpha)&(tab_nointer$log_12_inter_exo>threshold),na.rm=T)/nrow(tab_nointer))
print("1 causes 2, without")
print(sum(tab_nointer$Pval_21_inter_GC_exo<alpha,na.rm=T)/nrow(tab_nointer))
print(sum(tab_nointer$log_21_inter_exo>threshold,na.rm=T)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_inter_GC_exo<alpha)&(tab_nointer$log_21_inter_exo>threshold),na.rm=T)/nrow(tab_nointer))


########For CCM
tab_inter=read.csv("DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp1sp2factorized_CCM_otf.csv")

pdf("explo_with_driver_CCM_sp1sp2.pdf",width=15,height=10)
par(mfrow=c(1,3))
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal seasonal surr",names=c("sp1->sp2 inter","sp1->sp2 nointer","sp2->sp1 inter","sp2->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal permutation",names=c("sp1->sp2 inter","sp1->sp2 nointer","sp2->sp1 inter","sp2->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=colo,range=0,main="Rho",names=c("sp1->sp2 inter","sp1->sp2 nointer","sp2->sp1 inter","sp2->sp1 nointer"))
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)
dev.off()

print("FOR CCM table")
alpha=0.1
print("1 causes 2, with")
print(sum(tab_inter$Pval_12_inter_CCM_surr_season<alpha)/nrow(tab_inter))
print(sum(tab_inter$Rho_12>0.2)/nrow(tab_inter))
print(sum((tab_inter$Pval_12_inter_CCM_surr_season<alpha)&(tab_inter$Rho_12>0.2))/nrow(tab_inter))
print("2 causes 1, with")
print(sum(tab_inter$Pval_21_inter_CCM_surr_season<alpha)/nrow(tab_inter))
print(sum(tab_inter$Rho_21>0.2)/nrow(tab_inter))
print(sum((tab_inter$Pval_21_inter_CCM_surr_season<alpha)&(tab_inter$Rho_21>0.2))/nrow(tab_inter))
print("1 causes 2, without")
print(sum(tab_nointer$Pval_12_inter_CCM_surr_season<alpha)/nrow(tab_nointer))
print(sum(tab_nointer$Rho_12>0.2)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_12_inter_CCM_surr_season<alpha)&(tab_nointer$Rho_12>0.2))/nrow(tab_nointer))
print("2 causes 1, without")
print(sum(tab_nointer$Pval_21_inter_CCM_surr_season<alpha)/nrow(tab_nointer))
print(sum(tab_nointer$Rho_21>0.2)/nrow(tab_nointer))
print(sum((tab_nointer$Pval_21_inter_CCM_surr_season<alpha)&(tab_nointer$Rho_21>0.2))/nrow(tab_nointer))

########For CCM

pdf("explo_with_driver_CCM_sp1sp2temp_dummy.pdf",width=15,height=20)
tab_inter=read.csv("DataCompet_driver_intersp1tempfactorized_CCM_otf.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp1tempfactorized_CCM_otf.csv")
par(mfrow=c(2,3))
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal seasonal surr",names=c("sp1->temp inter","sp1->temp nointer","temp->sp1 inter","temp->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal permutation",names=c("sp1->temp inter","sp1->temp nointer","temp->sp1 inter","temp->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=colo,range=0,main="Rho",names=c("sp1->temp inter","sp1->temp nointer","temp->sp1 inter","temp->sp1 nointer"))
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)

tab_inter=read.csv("DataCompet_driver_intersp2tempfactorized_CCM_otf.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp2tempfactorized_CCM_otf.csv")
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal seasonal surr",names=c("sp2->temp inter","sp2->temp nointer","temp->sp2 inter","temp->sp2 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal permutation",names=c("sp2->temp inter","sp2->temp nointer","temp->sp2 inter","temp->sp2 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=colo,range=0,main="Rho",names=c("sp2->temp inter","sp2->temp nointer","temp->sp2 inter","temp->sp2 nointer"))
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)
dev.off()
