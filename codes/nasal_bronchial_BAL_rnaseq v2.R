##
# RNA seq analysis of nasal/bronchial data with edgeR and limma
# YJK local machine version
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)
library(gplots)
############
# count data and phenotype
############

# reading count data for batch1,2,3,4 and saving as 'counts'
file<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.exists(file)
counts<-read.delim(file[1])
counts.ID<-colnames(counts)

# reading count data for batch 5 as 'counts.b5'
file.b5<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")
file.exists(file.b5)
counts.b5<-read.delim(file.b5[1])
counds.b5ID<-colnames(counts.b5)


# load phenotype data by sourcing the following code 
source("./codes/phenotype_cleanup_nasal_bronchial_BAL_rnaseq.R")

# select just the nasal RNAseq counts
b<-counts.ID%in%nasal.ID

x<-counts[,b]
genes<-counts$SampleID
rownames(x)<-genes
x.BalNeut<-x[,c( p.count.BalNeut$SampleID)] # count table for DEG using BAL Neut information as predictor
x.SerCt<-x[,c(p.count.SerCt$SampleID)]# count table for DEG using serum cell counts information as predictor

# convert raw counts to CPM
# average coverage is ~51 million
cpm0 <-cpm(x)
lcpm0<-cpm(x,log=TRUE)
L<-mean(colSums(x))*1e-6 # mean library size
M<-median(colSums(x))*1e-6 # median library size


###########################################
# Normalising gene expression distributions
###########################################

#x3 is the TMM normalized count
# normalize counts with TMM
norm.factor<-calcNormFactors(x, method = "TMM")

x3<-x
sample.size<-length(colnames(x3))
for(i in 1:sample.size){
  x3[,i]<-x3[,i]/norm.factor[i]
}

# calculate lcpm based on TMM normalized counts 
lcpm.x3<-cpm(x3,log=TRUE)


##############
#MDS by groups
##############

# color coding based on each factor levels of sex, age, BAL Eos%, and batch sourced from the following code:
source("./codes/mds_color_code_nasal_bronchial_BAL_rnaseq.R")

# plot the MDS
par(mfrow=c(2,2))
mds.eos<-plotMDS(lcpm.x3,  col= col$eos, labels=p.counts$BAL_Eos_perc, main="Eos %")
mds.batch<-plotMDS(lcpm.x3,  col= col$batch, labels=p.counts$Batch, main="batch")
mds.age<-plotMDS(lcpm.x3,  col= col$age, labels=p.counts$Age, main = "age")
mds.race<-plotMDS(lcpm.x3,  col= col$race, labels=p.counts$Race_corrected, main = "race")

# MDS shows significant batch effect

# checking if covariates have batch effect:

cov.data<-p.counts[,c(33,3:22)]
cov.data<-cov.data[,c(1,grep("log",cov.data%>%colnames))]
vari<-colnames(cov.data[-1])

models <- lapply(vari,function(x){y<-as.formula(paste0(x,"~Batch"));return(y)})
aov.results<-lapply(models,function(x){aov(x,data=cov.data)%>%summary}) #> none of the covariates have a significant batch effect

# make a heat map
batch<-p.counts$Batch%>%as.character
hm.matrix<-as.matrix(cov.data[,-1], dimnames=list(batch,vari))
rownames(hm.matrix)<-batch
hm.matrix%>%replace(is.na(.),-1)%>%heatmap.2(srtCol = 25,adjCol = c(0.95,-0.5)) #> shows heatmap
##########################################
# make input dataframe for DEG with DESeq2
##########################################
library(DESeq2)

#### Removing genes that are lowly expressed

# checking how many genes have 0 count across all samples
table(rowSums(x==0)==45)
# about 4% of the genes have 0 counts across all samples 

# setting a lcpm cutoff for filtering genes with very low counts
lcpm.cutoff <- log2(10/M + 2/L)
dropCutoff<-function(cutoff){
  which(apply(lcpm.x3, 1, max) < cutoff)
}
drop <-dropCutoff(0) 
drop2<-dropCutoff(lcpm.cutoff)
dim(x3[-drop,])
dim(x3[-drop2,])

#define x2 as a new count table after filtering low expressed genes from x
x2<-x[-drop,]

x2.BalNeut<-x2[,c( p.count.BalNeut$SampleID)] # count table for DEG using BAL Neut information as predictor
x2.SerCt<-x2[,c(p.count.SerCt$SampleID)]# count table for DEG using serum cell counts information as predictor

#### find which are going to be predictor variables for gene expressions, i.e., names of the log transformed cell counts
# save as 'p.counts.var'

p.counts.var<-
  colnames(p.counts)[grep("log",colnames(p.counts))]


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
# run the DEG 
###################################################
# DEG
# df.deseq2input sourced from phenotype_cleanup_nasal_bronchial_BAL_rnaseq.R

for(i in 1:nrow(df.deseq2input)){
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}

# summarize number of significant genes in the 'res.table'
for(i in 1:nrow(res.table)){
  res.table[i,"sig.genes"]<-get(paste0("res",i))[[2]]%>%nrow
  res.table[i,"results"]<-paste0("res",i)
}
print(res.table)


# write results

for(i in 1:28){
  a<-get(paste0("res",i))
  write.csv(a[[1]],file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}

print(res.table)

write.csv(res.table,file.path(getwd(),"output",paste0("res.table_",Sys.Date(),".csv")))

# make list of significant DEG 
##############################

res.list=c(paste0("res",1:28))
gl<-numeric()
for(i in 1:28){
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

with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

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

# run the DEG 
###################################################
# DEG

ml<-list(zeroBalEos=model.matrix(~ zero.BALEos + Batch, p.counts),
         zero.BALNeut=model.matrix(~ zero.BALNeut + Batch, p.counts),
         zero.serEos=model.matrix(~ zero.serEos + Batch, p.counts),
         zero.serNeut=model.matrix(~ zero.serNeut + Batch, p.counts))
deseq2DEG.disc(df.deseq2input[29,1], df.deseq2input[29,2],ml[[1]],df.deseq2input[29,4])

for(i in 29:nrow(df.deseq2input)){
  assign(paste0("res",i),deseq2DEG.disc(df.deseq2input[i,1],df.deseq2input[i,2],ml[[i-28]],df.deseq2input[i,4]))
}
deseq2DEG.disc(df.deseq2input[32,1], df.deseq2input[32,2],ml[[4]],df.deseq2input[32,4]) #> not a full rank, because all samples have at least ser.Neut>0

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
res.table<-res.table[-32,]
res.table$fluid_cell[29:31]<-c("BAL_Eos","BAL_neut","serum_E")
res.table$units[29:31]<-c("0 vs >0")

print(res.table)

for(i in 1:nrow(res.table)){
  res.table[i,"sig.genes"]<-get(paste0("res",i))[[2]]%>%nrow
  res.table[i,"results"]<-paste0("res",i)
}
print(res.table)


res.list=c(paste0("res",1:31))
gl<-numeric()
for(i in 1:31){
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
for(i in 29:31){
  a<-get(paste0("res",i))
  write.csv(a[[1]],file.path(getwd(),"output",paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")))
  res.table[i,"output"]<-paste0("DEG_","res",i,"_",res.table[i,"fluid_cell"],"_",Sys.Date(),".csv")
}
