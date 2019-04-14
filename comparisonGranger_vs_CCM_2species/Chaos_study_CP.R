### CP April 2019
### We compare the results from Granger-causality and CCM on the same simulations (chaotic, stochastic, with a driver...) and try to compute a binary correlation phi.
### Note on phi : we write a function, but we can also use the library(psych) and the corresponding phi, or the sqrt(chisq.test$statistics/nsample)

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


print("1 causes 2")
plou=table(table_inter$index_1cause2_inter_GC,table_inter$index_1cause2_inter_CCM)
print("GC")
print(sum(table_inter$index_1cause2_inter_GC==1)/500)
print("CCM")
print(sum(table_inter$index_1cause2_inter_CCM==1)/500)
phi(plou)

print("2 causes 1")
plou=table(table_inter$index_2cause1_inter_GC,table_inter$index_2cause1_inter_CCM)
print("GC")
print(sum(table_inter$index_2cause1_inter_GC==1)/500)
print("CCM")
print(sum(table_inter$index_2cause1_inter_CCM==1)/500)
phi(plou)

