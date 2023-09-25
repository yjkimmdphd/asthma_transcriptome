###################################################################################
# analyzing overlap of DEG in different analysis in nasal_bronchial_BAL_rnaseq_v2.R
###################################################################################
library(dplyr)
deg.tab<-file.path(getwd(),"output/DEG_table_2023-09-20.csv")%>%read.csv

# vendiagram of DEG of serum Eos abs count and % as predictors, if including 0% = Eos%
deg.tab.serEos<-filter(deg.tab,results=="res6"|results=="res7"|results=="res16"|results=="res17")
deg.tab.serEos$genes%>%table
genelist<-list(
  ser_AEC1=filter(deg.tab,results=="res6")%>%select(genes),
  ser_Eos_p1=filter(deg.tab,results=="res7")%>%select(genes),
  ser_AEC2=filter(deg.tab,results=="res16")%>%select(genes),
  ser_Eos_p2=filter(deg.tab,results=="res17")%>%select(genes)
)

library(gplots)
v.table1<-venn(genelist[1:4])
v.table2<-venn(genelist[1:2])
v.table3<-venn(genelist[3:4])
print(v.table1)
print(v.table2)
print(v.table3)

# vendiagram of DEG of serum Eos abs count and % as predictors, if excluding 0% = Eos%
deg.tab.serEos<-filter(deg.tab,results=="res6"|results=="res7"|results=="res16"|results=="res17"|results=="res25"|results=="res26")
deg.tab.serEos$genes%>%table
genelist<-list(
  ser_AEC1=filter(deg.tab,results=="res6")%>%select(genes),
  ser_Eos_p1=filter(deg.tab,results=="res7")%>%select(genes),
  ser_AEC2=filter(deg.tab,results=="res16")%>%select(genes),
  ser_Eos_p2=filter(deg.tab,results=="res17")%>%select(genes),
  ser_AEC.pos=filter(deg.tab,results=="res25")%>%select(genes),
  ser_Eos_p.pos=filter(deg.tab,results=="res26")%>%select(genes)
  
)
par(mfrow=c(2,2))
library(gplots)
v.table1<-venn(genelist[1:4])
v.table2<-venn(genelist[1:2])
v.table3<-venn(genelist[3:4])
v.table4<-venn(genelist[5:6])
v.table5<-venn(genelist[c(1,3,5)])
v.table6<-venn(genelist[c(2,4,6)])
v.table7<-venn(genelist[c(3:6)])
v.table8<-venn(genelist[c(3,5)])
v.table9<-venn(genelist[c(4,6)])
print(v.table1)
print(v.table2)
print(v.table3)
print(v.table4)
print(v.table5)
print(v.table6)
print(v.table7)
print(v.table8)
print(v.table9)
# v.table 4:9 are interesting

genelist.tab8.9<-list(ser_AEC=attr(v.table8,"intersections")$`ser_AEC2:ser_AEC.pos`,
                      ser_Eos_p=attr(v.table9,"intersections")$`ser_Eos_p2:ser_Eos_p.pos`)
unlist(genelist.tab8.9)%>%as.vector()%>%write.csv(file.path(getwd(),"output/serEos_AEC_perc_gene_list.csv"))
