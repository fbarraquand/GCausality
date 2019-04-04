rm(list=ls())
graphics.off()

table_inter=read.csv("results/DataCompet_chaos_withinter.csv")

#1 cause 
plou=table_inter$index_1cause2_inter_GC+table_inter$index_1cause2_inter_CCM
#1 always causes 2 in this case

#2 cause
plou=table_inter$index_2cause1_inter_GC+table_inter$index_2cause1_inter_CCM
print(table(plou))
GCok=sum(table_inter$index_2cause1_inter_GC==1&table_inter$index_2cause1_inter_CCM==0)
CCMok=sum(table_inter$index_2cause1_inter_GC==0&table_inter$index_2cause1_inter_CCM==1)





table_no_inter=read.csv("results/DataCompet_chaos_withoutinter.csv")

#1 cause 
plou=table_no_inter$index_1cause2_inter_GC+table_no_inter$index_1cause2_inter_CCM
CCMok=sum(table_no_inter$index_1cause2_inter_GC==1&table_no_inter$index_1cause2_inter_CCM==0)
GCok=sum(table_no_inter$index_1cause2_inter_GC==0&table_no_inter$index_1cause2_inter_CCM==1)

#2 cause
plou=table_no_inter$index_2cause1_inter_GC+table_no_inter$index_2cause1_inter_CCM
print(table(plou))
CCMok=sum(table_no_inter$index_2cause1_inter_GC==1&table_no_inter$index_2cause1_inter_CCM==0)
GCok=sum(table_no_inter$index_2cause1_inter_GC==0&table_no_inter$index_2cause1_inter_CCM==1)



