###################################################################################
# analyzing overlap of DEG in different analysis in nasal_bronchial_BAL_rnaseq_v2.R
###################################################################################
library(dplyr)
sig.deg<-file.path(getwd(),"output/sig_deg_2023-09-28.csv")%>%read.csv
res.table<-file.path(getwd(),"output/res.table_2023-09-28.csv")%>%read.csv%>%na.omit
sig.deg<-left_join(sig.deg,res.table)%>%group_by(fluid_cell,results,units,count_data)
pred<-sig.deg$fluid_cell%>%unique
# vendiagram of DEG of serum Eos abs count and % as predictors, if including 0% = Eos%

library(gplots)

summarize(sig.deg,)

print(list(paste("number of predictors =",length(pred)),
           paste("predictors:", paste(pred,collapse=", ")))
      )

# significant DEG subset by predictor (BAL, serum, Eos, Neut)
deg.bn<-sig.deg%>%filter(fluid_cell==pred[1])
deg.bw<-sig.deg%>%filter(fluid_cell==pred[2])
deg.se<-sig.deg%>%filter(fluid_cell==pred[3])
deg.sn<-sig.deg%>%filter(fluid_cell==pred[4])
deg.be<-sig.deg%>%filter(fluid_cell==pred[5])

# make list of data for venndiagram of DEG 
bn.gl<-ungroup(deg.bn)%>%as.data.frame()%>%split(deg.bn$results)%>%lapply(function(d)select(d,genes)) #for BAL Neut
bw.gl<-ungroup(deg.bw)%>%as.data.frame()%>%split(deg.bw$results)%>%lapply(function(d)select(d,genes)) #for BAL WBC
se.gl<-ungroup(deg.se)%>%as.data.frame()%>%split(deg.se$results)%>%lapply(function(d)select(d,genes)) #for ser Eos
sn.gl<-ungroup(deg.sn)%>%as.data.frame()%>%split(deg.sn$results)%>%lapply(function(d)select(d,genes)) #for ser Neut
be.gl<-ungroup(deg.be)%>%as.data.frame()%>%split(deg.be$results)%>%lapply(function(d)select(d,genes)) #for BAL Eos

# make venndiagram
bn.v<-bn.gl[2:5]%>%venn
bw.v<-bw.gl%>%venn
se.v<-se.gl[c(1,2,5,6)]%>%venn
sn.v<-sn.gl[1:4]%>%venn
be.v<-be.gl%>%venn
attr(bn.v,"intersections")
deg.bn%>%summarise()%>%print #> count data used were when BAL Neut>0

# venndiagram of DEG for BAL Neut
bn.gl<-ungroup(deg.bn)%>%as.data.frame()%>%split(deg.bn$results)%>%lapply(function(d)select(d,genes))
bn.v<-bn.gl[2:5]%>%venn
attr(bn.v,"intersections")
attr(bw.v,"intersections")
attr(se.v,"intersections")
deg.bn%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.bw%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.se%>%summarise()%>%print #> count data used were when BAL Neut>0
deg.sn%>%summarise()%>%print #> count data used were when BAL Neut>0
