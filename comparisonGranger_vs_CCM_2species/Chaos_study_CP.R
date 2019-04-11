rm(list=ls())
graphics.off()


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

	phi=(n11*n00-n10*n01)/sqrt(1.0*n1p*n0p*np0*n1p) #the 1.0 to avoid integer overflow

	if(is.na(phi)){
		phi=(n*n11-n1p*np1)/sqrt(n1p*np1*(n-n1p)*(n-np1))
	}

	return(phi)


}

#table_inter=read.csv("results/DataCompet_chaos_withinter.csv")

#table_inter=read.csv("results/DataCompet_stochModel_inter.csv")

table_inter=read.csv("../twoSpecies_andDriver/DataCompet_driver_intersp1sp2.csv")


#1 cause 
plou=table(table_inter$index_1cause2_inter_GC,table_inter$index_1cause2_inter_CCM)
phi(plou)
#1 always causes 2 in this case

#2 cause
plou=table_inter$index_2cause1_inter_GC+table_inter$index_2cause1_inter_CCM
print(table(plou))
GCok=sum(table_inter$index_2cause1_inter_GC==1&table_inter$index_2cause1_inter_CCM==0)
CCMok=sum(table_inter$index_2cause1_inter_GC==0&table_inter$index_2cause1_inter_CCM==1)





table_no_inter=read.csv("results/DataCompet_chaos_withoutinter.csv")
#table_no_inter=read.csv("results/DataCompet_stochModel_noInter.csv")

#1 cause 
plou=table_no_inter$index_1cause2_inter_GC+table_no_inter$index_1cause2_inter_CCM
CCMok=sum(table_no_inter$index_1cause2_inter_GC==1&table_no_inter$index_1cause2_inter_CCM==0)
GCok=sum(table_no_inter$index_1cause2_inter_GC==0&table_no_inter$index_1cause2_inter_CCM==1)

#2 cause
plou=table_no_inter$index_2cause1_inter_GC+table_no_inter$index_2cause1_inter_CCM
print(table(plou))
CCMok=sum(table_no_inter$index_2cause1_inter_GC==1&table_no_inter$index_2cause1_inter_CCM==0)
GCok=sum(table_no_inter$index_2cause1_inter_GC==0&table_no_inter$index_2cause1_inter_CCM==1)


#Compute binary correlation
a=table(table_inter$index_2cause1_inter_GC,table_inter$index_2cause1_inter_CCM)
plou=chisq.test(a)
phi=sqrt(plou$statistic/500)
print(phi)


a=table(table_no_inter$index_1cause2_inter_GC,table_no_inter$index_1cause2_inter_CCM)
plou=chisq.test(a)
phi=sqrt(plou$statistic/500)
print(phi)

a=table(table_no_inter$index_2cause1_inter_GC,table_no_inter$index_2cause1_inter_CCM)
plou=chisq.test(a,,simulate.p.value=T)
phi=sqrt(plou$statistic/500)
print(phi)
