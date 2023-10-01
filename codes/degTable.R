########################
## exploring DEG results
########################
library(tidyverse)

# if individual DEG tables of each DEG results are available,
# read them and save them as individual data.frame 
# then summarize the results in res.table
fp<-file.path(getwd(),"output/DEG_2023-09-30")
if(file.exists(fp)){
  flist<-list.files(fp)[1:39]
  fl<-sapply(flist,function(d)str_split(d,"DEG_res"))%>%unlist%>%matrix(ncol=2, byrow=TRUE)
  fl<-fl[,2]%>%str_split("_")%>%unlist%>%matrix(ncol=4,byrow=TRUE)
  fl<-fl[,1]%>%as.numeric%>%na.omit%>%as.numeric
  for(i in fl){
    assign(paste0("res",i),read.csv(file.path(fp,flist[i]))%>%mutate(results=paste0("res",i))%>%rename(X="genes"))
    
    
  }
  # summarize number of significant genes in the 'res.table'
  cont.var<-grep("zero",res.table$fluid_cell, invert=TRUE)
  zero.var<-grep("zero",res.table$fluid_cell)
  
  for(i in c(cont.var,zero.var)){
    res.table[i,"sig.genes"]<-get(paste0("res",i))%>%filter(padj<0.05)%>%nrow
    res.table[i,"results"]<-paste0("res",i)
  }
  print(res.table)
  
  # write sig.table csv
  write.csv(res.table,file.path(getwd(),"output",paste0("res.table_",Sys.Date(),".csv")),row.names = FALSE)
  ##############################
  # make list of significant DEG 
  ##############################
  res.list=c(paste0("res",c(cont.var)))
  dt<-rbindlist((lapply(res.list,get)))
  sig.deg<-dt%>%filter(padj<0.05)%>%group_by(results)
  write.csv(
    sig.deg,file.path(getwd(),"output",paste0("sig_deg_",Sys.Date(),".csv")),row.names=FALSE
  )
} else{
  sig.deg<-file.path(getwd(),"output/sig_deg_2023-09-28.csv")%>%read.csv
  res.table<-file.path(getwd(),"output/res.table_2023-09-28.csv")%>%read.csv%>%na.omit
}

table(sig.deg$genes)%>%as.data.frame()%>%arrange(desc(Freq))%>%print
