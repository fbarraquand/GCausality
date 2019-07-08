########################################################################################################################
########### CP 19/04/2019 - Diagnostic of p-values and thresholds for both GC and CCM, can be used for sim with driver. ###########
########### CP 24/05/2019 - Added Sokal-Michener and Yule's indices for similarity
########### CP 08/07/2019 - Write results in a real table
########################################################################################################################

rm(list=ls())
graphics.off()

library(xtable)

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


colo=c("red","blue","orange","cyan")

#pdf("explo_with_driver_GC_tmax800.pdf",width=10,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,2,3,1))
#let start by GC
tab_GC=read.csv('DataCompet_driver_inter_factorized_GC_otf.csv')
#tab_GC=read.csv('DataCompet_driver_inter_factorized_GC_otf_tmax800.csv')
#tab_GC=read.csv('DataCompet_driver_inter_factorized_GC_otf_tmax800_withtherightmeany1.csv')
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
#dev.off()

table_to_write=matrix(NA,4,11)
rownames(table_to_write)=c("Inter12","Inter21","NoInter12","NoInter21")
colnames(table_to_write)=c("GCpvalpair","GCLRpair","GCbothpair","GCpvalcond","GCLRcond","GCbothcond","CCMpval","CCMrho2","both2","ISM","Q")

alpha=0.1
threshold=0.04

table_to_write[1,1]=sum(tab_inter$Pval_12_inter_GC_noexo<alpha)/nrow(tab_inter)
table_to_write[1,2]=sum(tab_inter$log_12_inter_noexo>threshold)/nrow(tab_inter)
table_to_write[1,3]=sum((tab_inter$Pval_12_inter_GC_noexo<alpha)&(tab_inter$log_12_inter_noexo>threshold))/nrow(tab_inter)

table_to_write[2,1]=sum(tab_inter$Pval_21_inter_GC_no_exo<alpha)/nrow(tab_inter)
table_to_write[2,2]=sum(tab_inter$log_21_inter_noexo>threshold)/nrow(tab_inter)
table_to_write[2,3]=sum((tab_inter$Pval_21_inter_GC_no_exo<alpha)&(tab_inter$log_21_inter_noexo>threshold))/nrow(tab_inter)

table_to_write[3,1]=sum(tab_nointer$Pval_12_inter_GC_noexo<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[3,2]=sum(tab_nointer$log_12_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[3,3]=sum((tab_nointer$Pval_12_inter_GC_noexo<alpha)&(tab_nointer$log_12_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer)

table_to_write[4,1]=sum(tab_nointer$Pval_21_inter_GC_no_exo<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[4,2]=sum(tab_nointer$log_21_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[4,3]=sum((tab_nointer$Pval_21_inter_GC_no_exo<alpha)&(tab_nointer$log_21_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer)

table_to_write[1,4]=sum(tab_inter$Pval_12_inter_GC_exo<alpha)/nrow(tab_inter)
table_to_write[1,5]=sum(tab_inter$log_12_inter_exo>threshold)/nrow(tab_inter)
table_to_write[1,6]=sum((tab_inter$Pval_12_inter_GC_exo<alpha)&(tab_inter$log_12_inter_exo>threshold))/nrow(tab_inter)
index_1cause2_inter_GC=(tab_inter$Pval_12_inter_GC_exo<alpha)*(tab_inter$log_12_inter_exo>threshold)

table_to_write[2,4]=sum(tab_inter$Pval_21_inter_GC_exo<alpha)/nrow(tab_inter)
table_to_write[2,5]=sum(tab_inter$log_21_inter_exo>threshold)/nrow(tab_inter)
table_to_write[2,6]=sum((tab_inter$Pval_21_inter_GC_exo<alpha)&(tab_inter$log_21_inter_exo>threshold))/nrow(tab_inter)
index_2cause1_inter_GC=(tab_inter$Pval_21_inter_GC_exo<alpha)*(tab_inter$log_21_inter_exo>threshold)

table_to_write[3,4]=sum(tab_nointer$Pval_12_inter_GC_exo<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[3,5]=sum(tab_nointer$log_12_inter_exo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[3,6]=sum((tab_nointer$Pval_12_inter_GC_exo<alpha)&(tab_nointer$log_12_inter_exo>threshold),na.rm=T)/nrow(tab_nointer)
index_1cause2_nointer_GC=(tab_nointer$Pval_12_inter_GC_exo<alpha)*(tab_nointer$log_12_inter_exo>threshold)

table_to_write[4,4]=sum(tab_nointer$Pval_21_inter_GC_exo<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[4,5]=sum(tab_nointer$log_21_inter_exo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[4,6]=sum((tab_nointer$Pval_21_inter_GC_exo<alpha)&(tab_nointer$log_21_inter_exo>threshold),na.rm=T)/nrow(tab_nointer)
index_2cause1_nointer_GC=(tab_nointer$Pval_21_inter_GC_exo<alpha)*(tab_nointer$log_21_inter_exo>threshold)

########For CCM
tab_inter=read.csv("DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp1sp2factorized_CCM_otf.csv")
#tab_inter=read.csv("DataCompet_driver_intersp1sp2factorized_CCM_otf_tmax800.csv")
#tab_nointer=read.csv("DataCompet_driver_noIntersp1sp2factorized_CCM_otf_tmax800.csv")
#tab_inter=read.csv("DataCompet_driver_intersp1sp2factorized_CCM_otf_tmax800_withtherightmeany1.csv")
#tab_nointer=read.csv("DataCompet_driver_noIntersp1sp2factorized_CCM_otf_tmax800_withtherightmeany1.csv")

pdf("explo_with_driver_CCM_sp1sp2_tmax800.pdf",width=15,height=10)
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

alpha=0.1
threshold=0.2

table_to_write[1,7]=sum(tab_inter$Pval_12_inter_CCM_surr_season<alpha)/nrow(tab_inter)
table_to_write[1,8]=sum(tab_inter$Rho_12>0.2)/nrow(tab_inter)
table_to_write[1,9]=sum((tab_inter$Pval_12_inter_CCM_surr_season<alpha)&(tab_inter$Rho_12>0.2))/nrow(tab_inter)
index_1cause2_inter_CCM=(tab_inter$Pval_12_inter_CCM_surr_season<alpha)*(tab_inter$Rho_12>threshold)

table_to_write[2,7]=sum(tab_inter$Pval_21_inter_CCM_surr_season<alpha)/nrow(tab_inter)
table_to_write[2,8]=sum(tab_inter$Rho_21>0.2)/nrow(tab_inter)
table_to_write[2,9]=sum((tab_inter$Pval_21_inter_CCM_surr_season<alpha)&(tab_inter$Rho_21>0.2))/nrow(tab_inter)
index_2cause1_inter_CCM=(tab_inter$Pval_21_inter_CCM_surr_season<alpha)*(tab_inter$Rho_21>threshold)

table_to_write[3,7]=sum(tab_nointer$Pval_12_inter_CCM_surr_season<alpha)/nrow(tab_nointer)
table_to_write[3,8]=sum(tab_nointer$Rho_12>0.2)/nrow(tab_nointer)
table_to_write[3,9]=sum((tab_nointer$Pval_12_inter_CCM_surr_season<alpha)&(tab_nointer$Rho_12>0.2))/nrow(tab_nointer)
index_1cause2_nointer_CCM=(tab_nointer$Pval_12_inter_CCM_surr_season<alpha)*(tab_nointer$Rho_12>threshold)

table_to_write[4,7]=sum(tab_nointer$Pval_21_inter_CCM_surr_season<alpha)/nrow(tab_nointer)
table_to_write[4,8]=sum(tab_nointer$Rho_21>0.2)/nrow(tab_nointer)
table_to_write[4,9]=sum((tab_nointer$Pval_21_inter_CCM_surr_season<alpha)&(tab_nointer$Rho_21>0.2))/nrow(tab_nointer)
index_2cause1_nointer_CCM=(tab_nointer$Pval_21_inter_CCM_surr_season<alpha)*(tab_nointer$Rho_21>threshold)

########For CCM

pdf("explo_with_driver_CCM_sp1sp2temp_dummy_tmax800.pdf",width=16,height=10)
tab_inter=read.csv("DataCompet_driver_intersp1tempfactorized_CCM_otf_tmax800.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp1tempfactorized_CCM_otf_tmax800.csv")
par(mfrow=c(2,3),cex=1.25,mar=c(3,3,1.5,0.5))
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal seasonal surr",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85)
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal permutation",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85)
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=colo,range=0,main="Rho",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85)
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)

tab_inter=read.csv("DataCompet_driver_intersp2tempfactorized_CCM_otf_tmax800.csv")
tab_nointer=read.csv("DataCompet_driver_noIntersp2tempfactorized_CCM_otf_tmax800.csv")
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal seasonal surr",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85)
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal permutation",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85)
p0=lapply(pvalz,function(x) sum(x==0))
text(p0,x=1:4,y=rep(-3))
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=colo,range=0,main="Rho",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85)
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)
dev.off()

### For phi
plou=table(index_1cause2_inter_GC,index_1cause2_inter_CCM)
table_to_write[1,10]=sk_index(plou)
table_to_write[1,11]=yule_index(plou)

plou=table(index_2cause1_inter_GC,index_2cause1_inter_CCM)
table_to_write[2,10]=sk_index(plou)
table_to_write[2,11]=yule_index(plou)

plou=table(index_1cause2_nointer_GC,index_1cause2_nointer_CCM)
table_to_write[3,10]=sk_index(plou)
table_to_write[3,11]=yule_index(plou)

plou=table(index_2cause1_nointer_GC,index_2cause1_nointer_CCM)
table_to_write[4,10]=sk_index(plou)
table_to_write[4,11]=yule_index(plou)

table_to_write[,1:9]=100*table_to_write[,1:9]

print.xtable(xtable(table_to_write,digits=c(1,rep(1,9),2,2)),"pval_threshold_2spdriver.tex" ,type="latex")
