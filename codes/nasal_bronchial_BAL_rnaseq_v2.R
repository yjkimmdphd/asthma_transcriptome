##
# RNA seq analysis of nasal/bronchial data with edgeR and limma
# YJK local machine version
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)
library(gplots)

# load phenotype data by sourcing the following code 
source("./codes/phenotype_cleanup_nasal_bronchial_BAL_rnaseq.R")

##############
#MDS by groups
##############

# color coding based on each factor levels of sex, age, BAL Eos%, and batch sourced from the following code:
source("./codes/mds_color_code_nasal_bronchial_BAL_rnaseq.R")

# plot the MDS
pdf(file= file.path(getwd(),"output/mds_eos_batch_age_race.pdf" ))
par(mfrow=c(2,2))
mds.eos<-plotMDS(lcpm.x3,  col= col$eos, labels=p.counts$BAL_Eos_perc, main="Eos %", plot=FALSE)
mds.batch<-plotMDS(lcpm.x3,  col= col$batch, labels=p.counts$Batch, main="batch",plot=FALSE)
mds.age<-plotMDS(lcpm.x3,  col= col$age, labels=p.counts$Age, main = "age", plot=FALSE)
mds.race<-plotMDS(lcpm.x3,  col= col$race, labels=p.counts$Race_corrected, main = "race",plot=FALSE)
dev.off()
# MDS shows significant batch effect

# checking if covariates have batch effect:

cov.data<-cbind(Batch=p.counts$Batch,p.counts[,c(3:23)])
cov.data<-cov.data[,c(1,grep("log",cov.data%>%colnames))]
vari<-colnames(cov.data[-1])

models <- lapply(vari,function(x){y<-as.formula(paste0(x,"~Batch"));return(y)})
aov.results<-lapply(models,function(x){aov(x,data=cov.data)%>%summary}) #> none of the covariates have a significant batch effect
capture.output(aov.results,file=file.path(getwd(),"output/cov_aov_assoc.txt")) #> write aov result

# make a heat map of each predictor variable, and see if they cluster based on batch
pdf(file= file.path(getwd(),"output/cov_batch_heatmap.pdf" ))
par(mfrow=c(1,1))
batch<-p.counts$Batch%>%as.character
hm.matrix<-as.matrix(cov.data[,-1], dimnames=list(batch,vari))
rownames(hm.matrix)<-batch
hm.matrix%>%replace(is.na(.),-1)%>%heatmap.2(srtCol = 25,adjCol = c(0.95,-0.5)) #> shows heatmap
dev.off()

##########################################
# make input dataframe for DEG with DESeq2
##########################################
library(DESeq2)

#### find which are going to be predictor variables for gene expressions, i.e., names of the log transformed cell counts
# save as 'p.counts.var'
p.counts.var<-
  colnames(p.counts)[grep("log",colnames(p.counts))]

##### make dataframe of input values that will be used for DEG with DESeq2
# save as 'df.deseq2input'
# x is original count table
# x2 is x filtered for low gene counts
df.deseq2input<-data.frame(count.data="x",
                           col.data="p.counts", 
                           design = paste0("~ ",p.counts.var," + Batch"), 
                           resoutput = p.counts.var)
df.deseq2input[c(3,4),1:2]<-data.frame(rep("x.BalNeut",2),rep("p.count.BalNeut",2))
df.deseq2input[6:10,1:2]<-data.frame(rep("x.SerCt",5),rep("p.count.SerCt",5))
df.input2<-df.deseq2input%>%mutate(count.data=sub("x","x2",count.data))
df.deseq2input<-rbind(df.deseq2input,df.input2)
print(df.deseq2input)


# update 'df.deseq2input'
input.pos<-data.frame(count.data=rep(c("x.BalEos.pos", "x.BalNeut.pos", "x.serEos.pos", "x.serNeut.pos",
                            "x2.BalEos.pos","x2.BalNeut.pos","x2.serEos.pos","x2.serNeut.pos"),
                          each=2))
input.pos$col.data<-rep(rep(c("p.BalEos.pos","p.BalNeut.pos","p.serEos.pos","p.serNeut.pos"),each=2),2)
input.pos$design<-rep((df.deseq2input[c(1:4,6:9),"design"]),2)
input.pos$resoutput<-rep((df.deseq2input[c(1:4,6:9),"resoutput"]),2)
df.deseq2input<-rbind(df.deseq2input,input.pos)
print(df.deseq2input)

input<-df.deseq2input
z.counts<-c("zero.BALEosnonZero","zero.BALNeutnonZero","zero.serEosnonZero")
z.input<-data.frame(count.data=input[c(11,13,16),"count.data"],
                    col.data=input[c(11,13,16),"col.data"],
                    design = paste0("~ ",z.counts," + Batch"),
                    resoutput = z.counts)
df.deseq2input<-rbind(input,z.input)
print(df.deseq2input)


#### summary of analysis is saved in 'res.table'
fluid.cell.sample<-sapply(df.deseq2input$resoutput,function(d){substring(d,1,7)})
perc<-grep("perc",df.deseq2input$resoutput)
fcp<-data.frame(fluid.cell=fluid.cell.sample)
fcp[perc,"unit"]<-"%"
fcp[-perc,"unit"]<-"abs count"
res.table<-data.frame(software="DESEq2",
                      fluid_cell=fcp[,1],
                      units=fcp[,2],
                      model=df.deseq2input$design,
                      sig.genes=0,
                      count_data=df.deseq2input$count.data)
print(res.table)

#######################################
# run the DEG for continuous predictors
#######################################

#### define function deseq2DEG 
deseq2DEG<-function(countdata,coldata,design,resultname){
  dds<-DESeqDataSetFromMatrix(countData = get(countdata),colData=get(coldata), design=as.formula(design))
  dds<-DESeq(dds)
  res<-results(dds, name=resultname)
  res <- res[order(res$padj),]
  # signficant results with padj<0.05
  res.sig<-res[which(res$padj<0.05),]
  return(list(res,res.sig))
}

zero.var<-grep("zero",res.table$fluid_cell)
cont.var<-grep("zero",res.table$fluid_cell, invert=TRUE)

for(i in cont.var){
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}

# summarize number of significant genes in the 'res.table'
for(i in cont.var){
  res.table[i,"sig.genes"]<-get(paste0("res",i))[[2]]%>%nrow
  res.table[i,"results"]<-paste0("res",i)
}
print(res.table)


# write results

for(i in cont.var){
  a<-get(paste0("res",i))
  write.csv(a[[1]],row.names=TRUE,file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}



print(res.table)

write.csv(res.table,file.path(getwd(),"output",paste0("res.table_",Sys.Date(),".csv")))
##############################
# make list of significant DEG 
##############################

res.list=c(paste0("res",cont.var))
gl<-numeric()
for(i in cont.var){
  gl[i]<-get(res.list[i])[[2]]%>%nrow
}
gl[which(gl==0)]<-1

deg.tab<-data.frame(results=rep(res.list,times=gl))
gn<-sapply(res.list,function(d){
  a<-get(d);b<-a[[2]]%>%rownames();
  b<-unlist(b)%>%as.vector();
  return(b)})
gFC<-sapply(res.list,function(d){
  a<-get(d);
  b<-a[[2]]$log2FoldChange;
  b<-unlist(b)%>%as.vector();
  return(b)})
gpadj<-sapply(res.list,function(d){
  a<-get(d);
  b<-a[[2]]$padj;
  b<-unlist(b)%>%as.vector();
  return(b)})
gnl<-which(lapply(gn,length)==0)
gn[gnl]<-"none"
gFC[gnl]<-NA
gpadj[gnl]<-NA
deg.tab$genes<-unlist(gn)%>%as.vector
deg.tab$log2FC<-unlist(gFC)%>%as.vector
deg.tab$padj<-unlist(gpadj)%>%as.vector

table(deg.tab$genes)%>%as.data.frame()%>%arrange(desc(Freq))

# save DEG table
write.csv(deg.tab,file.path(getwd(),"output",paste0("DEG_table_",Sys.Date(),".csv")))


#Make a basic volcano plot
 
#with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
#with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

########################################################
# DEG comparing zero vs nonZero Eos/Neut in BAL or serum
########################################################

#### define function deseq2DE.disc for ananlysis using model matrix containing discrete variable

deseq2DEG.disc<-function(countdata,coldata,des,resultname){
  dds<-DESeqDataSetFromMatrix(countData = get(countdata),colData=get(coldata), design=des)
  dds<-DESeq(dds)
  res<-results(dds, name=resultname)
  res <- res[order(res$padj),]
  # signficant results with padj<0.05
  res.sig<-res[which(res$padj<0.05),]
  return(list(res,res.sig))
}
####################################
# run the DEG on discrete predictors 
####################################
# DEG

# make model
ml<-list(zeroBalEos=model.matrix(~ zero.BALEos + Batch, p.counts),
         zero.BALNeut=model.matrix(~ zero.BALNeut + Batch, p.counts),
         zero.serEos=model.matrix(~ zero.serEos + Batch, p.counts),
         zero.serNeut=model.matrix(~ zero.serNeut + Batch, p.counts))

for(i in zero.var){
  assign(paste0("res",i),deseq2DEG.disc(df.deseq2input[i,1],df.deseq2input[i,2],ml[[i-length(cont.var)]],df.deseq2input[i,4]))
}

# update the 'res.table' which summarize number of significant genes
fluid.cell.sample<-sapply(df.deseq2input$resoutput,function(d){substring(d,1,7)})
perc<-grep("perc",df.deseq2input$resoutput)
fcp<-data.frame(fluid.cell=fluid.cell.sample)
fcp[perc,"unit"]<-"%"
fcp[-perc,"unit"]<-"abs count"
res.table<-data.frame(software="DESEq2",
                      fluid_cell=fcp[,1],
                      units=fcp[,2],
                      model=df.deseq2input$design,
                      sig.genes=0,
                      count_data=df.deseq2input$count.data)
res.table$fluid_cell[zero.var]<-c("BAL_Eos","BAL_neut","serum_E")
res.table$units[zero.var]<-c("0 vs >0")

print(res.table)

for(i in zero.var){
  res.table[i,"sig.genes"]<-get(paste0("res",i))[[2]]%>%nrow
  res.table[i,"results"]<-paste0("res",i)
}
print(res.table)


######
## write DEG table

res.list=c(paste0("res",1:nrow(res.table)))
gl<-numeric()
for(i in 1:nrow(res.table)){
  gl[i]<-get(res.list[i])[[2]]%>%nrow
}
gl[which(gl==0)]<-1

deg.tab<-data.frame(results=rep(res.list,times=gl))
gn<-sapply(res.list,function(d){
  a<-get(d);b<-a[[2]]%>%rownames();
  b<-unlist(b)%>%as.vector();
  return(b)})
gFC<-sapply(res.list,function(d){
  a<-get(d);
  b<-a[[2]]$log2FoldChange;
  b<-unlist(b)%>%as.vector();
  return(b)})
gpadj<-sapply(res.list,function(d){
  a<-get(d);
  b<-a[[2]]$padj;
  b<-unlist(b)%>%as.vector();
  return(b)})
gnl<-which(lapply(gn,length)==0)
gn[gnl]<-"none"
gFC[gnl]<-NA
gpadj[gnl]<-NA
deg.tab$genes<-unlist(gn)%>%as.vector
deg.tab$log2FC<-unlist(gFC)%>%as.vector
deg.tab$padj<-unlist(gpadj)%>%as.vector

table(deg.tab$genes)%>%as.data.frame()%>%arrange(desc(Freq))

# save results as CSV files
write.csv(deg.tab,file.path(getwd(),"output",paste0("DEG_table_",Sys.Date(),".csv")))
write.csv(res.table,file.path(getwd(),"output",paste0("res.table_",Sys.Date(),".csv")))
for(i in zero.var){
  a<-get(paste0("res",i))
  write.csv(a[[1]],file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}

