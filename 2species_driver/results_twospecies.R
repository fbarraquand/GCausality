########################################################################################################################
########### CP 19/04/2019 - Diagnostic of p-values and thresholds for both GC and CCM, can be used for sim with driver. ###########
########### CP 24/05/2019 - Added Sokal-Michener and Yule's indices for similarity
########### CP 08/07/2019 - Write results in a real table
########### CP 22/07/2019 - Remove the Q index
########### CP 18/06/2020 - Adapted the code to the new files (with Wald and F-test for Granger-Causality) and added the comparison between Wald and Ftest
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


#tab_GC=read.csv('results/DataCompet_driver_inter_factorized_GC_otf.csv')
#tab_GC=read.csv('results/DataCompet_driver_intersp1sp2factorized_GC_otf_with_F_Wald_test.csv')
tab_GC=read.csv('results/DataCompet_driver_intersp1sp2factorized_GC_otf_with_F_Wald_test.csv')
#tab_GC=read.csv('results/explo/DataCompet_driver_inter_factorized_GC_otf_tmax800.csv')
tab_inter=tab_GC[1:500,]
tab_nointer=tab_GC[501:1000,]

colo=c("red","blue","orange","cyan")

pdf("fig/explo_with_driver_GC.pdf",width=12,height=10)
par(mfrow=c(2,2),cex=1.,mar=c(4,4.2,2,1))
#let start by GC
logz=data.frame(tab_inter$log_12_inter_exo,tab_inter$log_12_inter_noexo,tab_nointer$log_12_inter_exo,tab_nointer$log_12_inter_noexo)
boxplot(logz,col=c("white","blue","white","cyan"),border=c("blue","black","cyan","black"),range=0,ylab=expression("G"["1->2"]),main="",names=c("inter cond.","inter pair.","no inter cond.","no inter pair."),ylim=c(0,0.25),cex.lab=1.2)
abline(h=0.04)
mtext("a)",side=2,las=2,at=0.25*1.1)

logz=data.frame(tab_inter$log_21_inter_exo,tab_inter$log_21_inter_noexo,tab_nointer$log_21_inter_exo,tab_nointer$log_21_inter_noexo)
boxplot(logz,col=c("white","red","white","orange"),border=c("red","black","orange","black"),range=0,ylab=expression("G"["2->1"]),main="",names=c("inter cond.","inter pair.","no inter cond.","no inter pair."),ylim=c(0,0.25),cex.lab=1.2)
abline(h=0.04)
mtext("b)",side=2,las=2,at=0.25*1.1)

pvalz=data.frame(tab_inter$Pval_12_inter_GC_exo,tab_inter$Pval_12_inter_GC_noexo_Ftest,tab_inter$Pval_12_inter_GC_noexo_Wald,tab_nointer$Pval_12_inter_GC_exo,tab_nointer$Pval_12_inter_GC_noexo_Ftest,tab_nointer$Pval_12_inter_GC_noexo_Wald)
#boxplot(log10(pvalz),ylim=c(-5,0),col=colo,range=0,main="PVal ratio 1->2",names=c("inter conditional","inter pairwise F","interpairwise W","no inter conditional","no inter pairwise F","no inter pairwise W"))
p0=lapply(pvalz,function(x) sum(x==0))
a=as.matrix(log10(pvalz))
a[is.infinite(a)]=NA
#boxplot(a,ylim=c(-5,0),col=colo,range=0,main="PVal ratio 1->2",names=c("inter conditional","inter pairwise F","interpairwise W","no inter conditional","no inter pairwise F","no inter pairwise W"),width=c(0.75,0.35,0.35,0.75,0.35,0.35),at=c(1,1.8,2.2,3,3.8,4.2))
boxplot(a,ylim=c(-5,0),col=c("white","darkblue","deepskyblue3","white","cyan3","cyan4"),border=c("blue","black","black","cyan","black","black"),range=0,main="",ylab=expression(paste("log"[10],"(p-value"["1->2"],")")),names=c("inter","inter","inter","no inter","no inter","no inter"),width=c(0.75,0.35,0.35,0.75,0.35,0.35),at=c(1,1.8,2.2,3,3.8,4.2),cex.lab=1.2)
mtext(text=c("cond.","pairF","pairW","cond.","pairF","pairW"),side=1,line=2,at=c(1,1.75,2.25,3,3.75,4.25))
abline(h=-1)
mtext("c)",side=2,las=2,at=0.5)

if(Reduce("+",p0)>0){
text(p0,x=c(1,1.8,2.2,3,3.8,4.2)+0.05,y=rep(-5))
}


#pvalz=data.frame(tab_inter$Pval_21_inter_GC_exo,tab_inter$Pval_21_inter_GC_no_exo,tab_nointer$Pval_21_inter_GC_exo,tab_nointer$Pval_21_inter_GC_no_exo)
pvalz=data.frame(tab_inter$Pval_21_inter_GC_exo,tab_inter$Pval_21_inter_GC_no_exo_Ftest,tab_inter$Pval_21_inter_GC_no_exo_Wald,tab_nointer$Pval_21_inter_GC_exo,tab_nointer$Pval_21_inter_GC_no_exo_Ftest,tab_nointer$Pval_21_inter_GC_no_exo_Wald)
p0=lapply(pvalz,function(x) sum(x==0))
a=as.matrix(log10(pvalz))
a[is.infinite(a)]=NA
boxplot(a,ylim=c(-5,0),col=c("white","darkred","firebrick3","white","darkorange","darkorange3"),border=c("red","black","black","orange","black","black"),range=0,main="",ylab=expression(paste("log"[10],"(p-value"["2->1"],")")),names=c("inter","inter","inter","no inter","no inter","no inter"),width=c(0.75,0.35,0.35,0.75,0.35,0.35),at=c(1,1.8,2.2,3,3.8,4.2),cex.lab=1.2)
mtext(text=c("cond.","pairF","pairW","cond.","pairF","pairW"),side=1,line=2,at=c(1,1.75,2.25,3,3.75,4.25))
mtext("d)",side=2,las=2,at=0.5)
abline(h=-1)
if(Reduce("+",p0)>0){
text(p0,x=c(1,1.8,2.2,3,3.8,4.2)+0.05,y=rep(-5))
}
dev.off()

pdf("fig/explo_with_driver_GC_temperature_as_exogen_onepanel.pdf",width=8,height=6)
par(mfrow=c(1,1),cex=1.,mar=c(4,4,1,1))
effect_1_signif_inter=tab_inter$effect_exo1[tab_inter$Pval_exo1<=0.1]
effect_1_notsignif_inter=tab_inter$effect_exo1[tab_inter$Pval_exo1>0.1]
effect_2_signif_inter=tab_inter$effect_exo2[tab_inter$Pval_exo2<=0.1]
effect_2_notsignif_inter=tab_inter$effect_exo2[tab_inter$Pval_exo2>0.1]

effect_1_signif_nointer=tab_nointer$effect_exo1[tab_nointer$Pval_exo1<=0.1]
effect_1_notsignif_nointer=tab_nointer$effect_exo1[tab_nointer$Pval_exo1>0.1]
effect_2_signif_nointer=tab_nointer$effect_exo2[tab_nointer$Pval_exo2<=0.1]
effect_2_notsignif_nointer=tab_nointer$effect_exo2[tab_nointer$Pval_exo2>0.1]

effect=list(effect_1_signif_inter,effect_1_notsignif_inter,effect_2_signif_inter,effect_2_notsignif_inter,effect_1_signif_nointer,effect_1_notsignif_nointer,effect_2_signif_nointer,effect_2_notsignif_nointer)
boxplot(effect,range=0,ylim=range(c(unlist(effect)),0.5),border=c("black","blue","black","red","black","cyan","black","orange"),col=c("blue","white","red","white","cyan","white","orange","white"),names=c("T->1 inter","T->1 inter","T->2 inter","T->2 inter","T->1 no inter","T->1 no inter","T->2 no inter","T->2 no inter"),ylab="Temperature regression coefficient",xlab="")
mtext(rep(c("p<0.1","p>0.1"),4),side=1,line=2,at=1:8)
abline(h=0.5,lwd=2,col="black")
abline(h=0,lty=2)
dev.off()


pdf("fig/explo_with_driver_GC_temperature_as_exogen.pdf",width=10,height=6)
par(mfrow=c(1,2),cex=1.,mar=c(4,2,3,1))
effect=data.frame(tab_inter$effect_exo1,tab_inter$effect_exo2,tab_nointer$effect_exo1,tab_nointer$effect_exo2)
boxplot(effect,col=c(colT1I,colT2I,colT1NI,colT2NI),range=0,main="effect of temperature",names=c("temp->1 inter.","temp->2 inter","temp->1 no inter.","temp->2 no inter."))
#mtext("a)",side=2,las=2,at=max(logz)*1.1)

pvalz=data.frame(tab_inter$Pval_exo1,tab_inter$Pval_exo2,tab_nointer$Pval_exo1,tab_nointer$Pval_exo2)
p0=lapply(pvalz,function(x) sum(x==0))
a=as.matrix(log10(pvalz))
a[is.infinite(a)]=NA
boxplot(a,col=c(colT1I,colT2I,colT1NI,colT2NI),range=0,main="Pvalue of temperature",names=c("temp->1 inter.","temp->2 inter","temp->1 no inter.","temp->2 no inter."),ylim=c(-5,0.))
abline(h=-1)
dev.off()

###Effect of temperature with GC tests
#sp1temp
#tab_GC_sp1temp=read.csv('results/DataCompet_driver_intersp1tempfactorized_GC_otf_with_F_Wald_test.csv')
tab_GC_sp1temp=read.csv('results/DataCompet_driver_intersp1tempfactorized_GC_otf_with_F_Wald_test.csv')
sp1tempF_inter=tab_GC_sp1temp$Pval_21_inter_GC_no_exo_Ftest[1:500]
sp1tempW_inter=tab_GC_sp1temp$Pval_21_inter_GC_no_exo_Wald[1:500]
sp1tempF_nointer=tab_GC_sp1temp$Pval_21_inter_GC_no_exo_Ftest[501:1000]
sp1tempW_nointer=tab_GC_sp1temp$Pval_21_inter_GC_no_exo_Wald[501:1000]

#sp2temp
#tab_GC_sp2temp=read.csv('results/DataCompet_driver_intersp2tempfactorized_GC_otf_with_F_Wald_test.csv')
tab_GC_sp2temp=read.csv('results/DataCompet_driver_intersp2tempfactorized_GC_otf_with_F_Wald_test.csv')
sp2tempF_inter=tab_GC_sp2temp$Pval_21_inter_GC_no_exo_Ftest[1:500]
sp2tempW_inter=tab_GC_sp2temp$Pval_21_inter_GC_no_exo_Wald[1:500]
sp2tempF_nointer=tab_GC_sp2temp$Pval_21_inter_GC_no_exo_Ftest[501:1000]
sp2tempW_nointer=tab_GC_sp2temp$Pval_21_inter_GC_no_exo_Wald[501:1000]

pdf("fig/explo_GC_temperature_on_sp1_sp2_test.pdf",width=10,height=10)
pvalz=data.frame(sp1tempF_inter,sp1tempW_inter,sp1tempF_nointer,sp1tempW_nointer,sp2tempF_inter,sp2tempW_inter,sp2tempF_nointer,sp2tempW_nointer)
boxplot(log10(pvalz),ylim=c(-10,0),col=c("darkblue","deepskyblue3","cyan3","cyan4"),range=0,main="PVal temp->sp",names=c("sp1 inter","sp1 inter","sp1 nointer","sp1 nointer","sp2 inter","sp2 inter","sp2 nointer","sp2 nointer"))
mtext(text=c("Ftest","Waldtest","Ftest","Waldtest","Ftest","Waldtest","Ftest","Waldtest"),side=1,line=2,at=1:8)
abline(h=-1)
abline(v=4.5)
dev.off()


table_to_write=matrix(NA,4,10)
rownames(table_to_write)=c("Inter12","Inter21","NoInter12","NoInter21")
colnames(table_to_write)=c("GCpvalpair","GCLRpair","GCbothpair","GCpvalcond","GCLRcond","GCbothcond","CCMpval","CCMrho2","both2","ISM")

alpha=0.1
threshold=0.04

table_to_write[1,1]=sum(tab_inter$Pval_12_inter_GC_noexo_Ftest<alpha)/nrow(tab_inter)
table_to_write[1,2]=sum(tab_inter$log_12_inter_noexo>threshold)/nrow(tab_inter)
table_to_write[1,3]=sum((tab_inter$Pval_12_inter_GC_noexo_Ftest<alpha)&(tab_inter$log_12_inter_noexo>threshold))/nrow(tab_inter)

table_to_write[2,1]=sum(tab_inter$Pval_21_inter_GC_no_exo_Ftest<alpha)/nrow(tab_inter)
table_to_write[2,2]=sum(tab_inter$log_21_inter_noexo>threshold)/nrow(tab_inter)
table_to_write[2,3]=sum((tab_inter$Pval_21_inter_GC_no_exo_Ftest<alpha)&(tab_inter$log_21_inter_noexo>threshold))/nrow(tab_inter)

table_to_write[3,1]=sum(tab_nointer$Pval_12_inter_GC_noexo_Ftest<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[3,2]=sum(tab_nointer$log_12_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[3,3]=sum((tab_nointer$Pval_12_inter_GC_noexo_Ftest<alpha)&(tab_nointer$log_12_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer)

table_to_write[4,1]=sum(tab_nointer$Pval_21_inter_GC_no_exo_Ftest<alpha,na.rm=T)/nrow(tab_nointer)
table_to_write[4,2]=sum(tab_nointer$log_21_inter_noexo>threshold,na.rm=T)/nrow(tab_nointer)
table_to_write[4,3]=sum((tab_nointer$Pval_21_inter_GC_no_exo_Ftest<alpha)&(tab_nointer$log_21_inter_noexo>threshold),na.rm=T)/nrow(tab_nointer)

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
#tab_inter=read.csv("results/DataCompet_driver_intersp1sp2factorized_CCM_otf.csv")
#tab_nointer=read.csv("results/DataCompet_driver_noIntersp1sp2factorized_CCM_otf.csv")
#tab_inter=read.csv("results/DataCompet_driver_intersp1sp2factorized_CCM_otf_test.csv")
#tab_nointer=read.csv("results/DataCompet_driver_noIntersp1sp2factorized_CCM_otf_test.csv")
#tab_inter=read.csv("results/explo/DataCompet_driver_intersp1sp2factorized_CCM_otf_tmax800.csv")
#tab_nointer=read.csv("results/explo/DataCompet_driver_noIntersp1sp2factorized_CCM_otf_tmax800.csv")

tab_inter=read.csv("results/DataCompet_driver_intersp1sp2factorized_CCM_otf_test.csv")
tab_nointer=read.csv("results/DataCompet_driver_noIntersp1sp2factorized_CCM_otf_test.csv")


pdf("fig/explo_with_driver_CCM_sp1sp2.pdf",width=15,height=10)
par(mfrow=c(1,3))
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=colo,range=0,main="PVal seasonal surr",names=c("sp1->sp2 inter","sp1->sp2 nointer","sp2->sp1 inter","sp2->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=colo,range=0,main="PVal permutation",names=c("sp1->sp2 inter","sp1->sp2 nointer","sp2->sp1 inter","sp2->sp1 nointer"))
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
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

pdf("fig/explo_with_driver_CCM_sp1sp2temp_dummy.pdf",width=16,height=10)


#####COLORS FROM FIG S9
col1TI="blue"
col1TNI=rgb(155/256,79/256,150/256,1)
col2TI="blue"
col2TNI=rgb(155/256,79/256,150/256,1)

colT1I="red"
colT1NI=rgb(255/256,165/256,0,1)
colT2I="red"
colT2NI=rgb(255/256,165/256,0,1)

#tab_inter=read.csv("results/DataCompet_driver_intersp1tempfactorized_CCM_otf_test.csv")
#tab_nointer=read.csv("results/DataCompet_driver_noIntersp1tempfactorized_CCM_otf_test.csv")

tab_inter=read.csv("results/DataCompet_driver_intersp1tempfactorized_CCM_otf_test.csv")
tab_nointer=read.csv("results/DataCompet_driver_noIntersp1tempfactorized_CCM_otf_test.csv")

par(mfrow=c(2,3),cex=1.25,mar=c(3,4,1.5,0.5))
pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=c(col1TI,col1TNI,colT1I,colT1NI),range=0,main="Seasonal surr",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85,ylab=expression(paste("log"[10],"(p-value)")))
mtext("a)",side=2,las=2,at=0.2,line=1.5)
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=c(col1TI,col1TNI,colT1I,colT1NI),range=0,main="Permutation",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85)
mtext("b)",side=2,las=2,at=0.2,line=1.5)
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=c(col1TI,col1TNI,colT1I,colT1NI),range=0,main="Rho",names=c("s1->T +I","s1->T -I","T->s1 +I","T->s1 -I"),cex.axis=0.85)
mtext("c)",side=2,las=2,at=max(rhoz)*1.1,line=1.5)
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)

#tab_inter=read.csv("results/DataCompet_driver_intersp2tempfactorized_CCM_otf_test.csv")
#tab_nointer=read.csv("results/DataCompet_driver_noIntersp2tempfactorized_CCM_otf_test.csv")

tab_inter=read.csv("results/DataCompet_driver_intersp2tempfactorized_CCM_otf_test.csv")
tab_nointer=read.csv("results/DataCompet_driver_noIntersp2tempfactorized_CCM_otf_test.csv")

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_season,tab_nointer$Pval_12_inter_CCM_surr_season,tab_inter$Pval_21_inter_CCM_surr_season,tab_nointer$Pval_21_inter_CCM_surr_season)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=c(col2TI,col2TNI,colT2I,colT2NI),range=0,main="",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85,ylab=expression(paste("log"[10],"(p-value)")))
mtext("d)",side=2,las=2,at=0.2,line=1.5)
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=-1)

pvalz=data.frame(tab_inter$Pval_12_inter_CCM_surr_sample,tab_nointer$Pval_12_inter_CCM_surr_sample,tab_inter$Pval_21_inter_CCM_surr_sample,tab_nointer$Pval_21_inter_CCM_surr_sample)
boxplot(log10(pvalz),ylim=c(-2.25,0),col=c(col2TI,col2TNI,colT2I,colT2NI),range=0,main="",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85)
mtext("e)",side=2,las=2,at=0.2,line=1.5)
p0=lapply(pvalz,function(x) sum(x==0))
if(Reduce("+",p0)>0){
text(p0,x=1:4,y=rep(-3))
}
abline(h=-1)

rhoz=data.frame(tab_inter$Rho_12,tab_nointer$Rho_12,tab_inter$Rho_21,tab_nointer$Rho_21)
boxplot(rhoz,col=c(col2TI,col2TNI,colT2I,colT2NI),range=0,main="",names=c("s2->T +I","s2->T -I","T->s2 +I","T->s2 -I"),cex.axis=0.85)
mtext("f)",side=2,las=2,at=max(rhoz)*1.1,line=1.5)
abline(h=0.1,lty=3)
abline(h=0.2,lty=3)
dev.off()

### For phi
plou=table(index_1cause2_inter_GC,index_1cause2_inter_CCM)
table_to_write[1,10]=sk_index(plou)
#table_to_write[1,11]=yule_index(plou)

plou=table(index_2cause1_inter_GC,index_2cause1_inter_CCM)
table_to_write[2,10]=sk_index(plou)
#table_to_write[2,11]=yule_index(plou)

plou=table(index_1cause2_nointer_GC,index_1cause2_nointer_CCM)
table_to_write[3,10]=sk_index(plou)
#table_to_write[3,11]=yule_index(plou)

plou=table(index_2cause1_nointer_GC,index_2cause1_nointer_CCM)
table_to_write[4,10]=sk_index(plou)
#table_to_write[4,11]=yule_index(plou)

table_to_write[,1:9]=100*table_to_write[,1:9]

print.xtable(xtable(table_to_write,digits=c(1,rep(1,9),2)),"results/pval_threshold_2spdriver_test.tex" ,type="latex")
