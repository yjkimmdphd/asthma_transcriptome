##
# RNA seq analysis of nasal/bronchial data with edgeR and limma
# YJK local machine version
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)

# reading count data for batch1,2,3,4 and saving as 'counts'
file<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/original_data/MS_asthma/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.exists(file)
counts<-read.delim(file[1])
counts.ID<-colnames(counts)

# reading count data for batch 5 as 'counts.b5'
file.b5<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/original_data/MS_asthma/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")
file.exists(file.b5)
counts.b5<-read.delim(file.b5[1])
counds.b5ID<-colnames(counts.b5)

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/revised_data/asthma-phenotype-filtered-revised.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

setwd("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich")
filename3<-file.path(getwd(),"data/original_data/MS_asthma/MS_asthma_phenotype.batch1234.txt")
file.exists(filename3)
batch.info<-read.delim(filename3)

# calculate absolute Eos from total BAL WBC and BAL Eos %
phenotype<-phenotype%>%mutate(bal_Eos_ct=BAL_Eos_perc*0.01*BAL_WBC)
phenotype$bal_Eos_ct_log<-log(phenotype$bal_Eos_ct+0.001)

# log transform the BAL eos percentage
pseudo.eos<-(phenotype$BAL_Eos_perc+0.001)
phenotype$BAL_Eos_perc_log<-log(pseudo.eos)

# log transform bal wbc count
BAL_WBC_log<-log10(phenotype$BAL_WBC+0.00001)
phenotype$BAL_WBC_log<-BAL_WBC_log

#find which phenotype data columns have eosinophil information
eos.col<-grep("Eos",colnames(phenotype))
head(phenotype[,eos.col])
eos.colnames<-colnames(phenotype[,eos.col])

# subset  'phenotype' with BAL and blood cell count information to 'phenotype.eos'
phenotype.eos<-phenotype[,c("subject_assgn","BAL_WBC", "BAL_WBC_log",eos.colnames)]
phenotype.eos<-phenotype.eos[,c("subject_assgn",colnames(phenotype.eos)[-1]%>%sort)]

# find subject assignment ID with nasal and bronchial cell RNAseq data  
nasal.ID<-phenotype[,"subject_assgn"] # nasal study ID that have been collected
bronch.ID<-sub("N","B",nasal.ID) # bronchial subject ID 
nb.ID<-c(nasal.ID,bronch.ID) # nasal and bronchial subject ID
a<-counts.ID%in%nb.ID  
nb.exist.ID<-counts.ID[a] # subject ID nasal and bronchial samples that are in the count matrix

# subset phenotype data to p.eos
sample.ID<-counts.ID[counts.ID%in%nasal.ID]
p.eos<-phenotype.eos[nasal.ID%in%sample.ID,] # p.eos is the Eos% in BAL 
p.eos$SampleID<-p.eos$subject_assgn


# subset batch info filtered based on nb.exist.ID.
batch.info.ID<-which(batch.info$SampleID%in%nb.exist.ID)
batch.info.BAL<-batch.info[batch.info.ID,]
(batch.info.BAL$Type=="Nasal")%>%table

p.eos<-left_join(p.eos,batch.info.BAL, by="SampleID")
p.eos$Batch<-factor(p.eos$Batch, levels=unique(batch.info.BAL$Batch))
p.eos$Race_corrected<-factor(p.eos$Race_corrected, levels=p.eos$Race_corrected%>%unique)


# select just the nasal RNAseq counts
b<-counts.ID%in%nasal.ID
x<-counts[,b]
genes<-counts$SampleID
rownames(x)<-genes

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

# unsupervised clustering of samples based on arbitrary BAL eos%



e1<-p.eos$BAL_Eos_perc>=3
e2<-p.eos$BAL_Eos_perc<1
e3<-p.eos$BAL_Eos_perc>10

# color code covariates 
col<-data.frame(matrix(nrow=45,ncol=6))
colnames(col)<-c("eos","batch","age","sex","race","AR")
for(i in 1:6) {col[,i]<-rep("black",45)}

col[which(e1==TRUE),"eos"]<-"red" # BAL Eos >= 3% are red
col[which(e2==TRUE),"eos"]<-"green" # BAL Eos = 0% are green
col[which(e3==TRUE),"eos"]<-"blue"# BAL Eos >10% are blue

col[which(p.eos$Batch=="batch1"),"batch"]<-"red" # batch1 are red
col[which(p.eos$Batch=="batch2"),"batch"]<-"green" # batch2 are green
col[which(p.eos$Batch=="batch3"),"batch"]<-"blue"# batch3 are blue. Batch4 are black

col[which(p.eos$Age<=7),"age"]<-"red" # age less than or equal to 7 are red
col[which(p.eos$Age>7 & p.eos$Age<=13),"age"]<-"green" # age 8yo to 13yo are green
col[which(p.eos$Age>13 & p.eos$Age<=19),"age"]<-"blue"# age 14yo to 19yo are blue. age >19 are black

col[which(p.eos$Gender=="Male"),"sex"]<-"blue"
col[which(p.eos$Gender=="Female"),"sex"]<-"red"

Race_corrected<-p.eos$Race_corrected%>%unique
colors<-c("red","blue","green","black","orange","purple")

race.color<-data.frame(Race_corrected,colors)
race.color<-full_join(p.eos,race.color, by="Race_corrected")
col$race<-race.color$colors

par(mfrow=c(2,2))
mds.eos<-plotMDS(lcpm.x3,  col= col$eos, labels=p.eos$BAL_Eos_perc, main="Eos %")
mds.batch<-plotMDS(lcpm.x3,  col= col$batch, labels=p.eos$Batch, main="batch")
mds.age<-plotMDS(lcpm.x3,  col= col$age, labels=p.eos$Age, main = "age")
mds.race<-plotMDS(lcpm.x3,  col= col$race, labels=p.eos$Race_corrected, main = "race")

#########################################
# Removing genes that are lowly expressed
#########################################

# checking how many genes have 0 count across all samples
table(rowSums(x==0)==45)
# about 4% of the genes have 0 counts across all samples 

# setting a lcpm cutoff for filtering genes with very low counts
lcpm.cutoff <- log2(10/M + 2/L)
dropCutoff<-function(cutoff){
  which(apply(lcpm.x3, 1, max) < cutoff)
}
drop <-dropCutoff(lcpm.cutoff) 
dim(x3[-drop,])

########################
#DEG with limma and voom
########################

##### model with BAL Eos percentage 
bal.eos.log<-p.eos$BAL_Eos_perc_log
batch<-factor(p.eos$Batch, levels=unique(batch.info.BAL$Batch))
mm2<-model.matrix(~0+bal.eos.log + batch)
mm2.1<-model.matrix(~bal.eos.log + batch)
# check voom to see if different filtering needed
y<-voom(x3[-drop,],mm2,plot=T)

# may be needs more stringent filtering after define new cutoff point 
drop2<-dropCutoff(0)
y2<-voom(x3[-drop2,], mm2, plot=T)

fit<-lmFit(y2,mm2)
tmp<-contrasts.fit(fit, coef=1)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_balEosPerc_",Sys.Date(),".csv")))

##### model with absolute BAL Eos count 
bal.eos.ct<-p.eos$bal_Eos_ct_log
batch<-factor(p.eos$Batch, levels=unique(batch.info.BAL$Batch))
mm3<-model.matrix(~bal.eos.ct + batch)

# filter using cutoff 0 
drop2<-dropCutoff(0)
y3<-voom(x3[-drop2,], mm3, plot=T)

fit<-lmFit(y3,mm3)
tmp<-contrasts.fit(fit, coef=2)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_balEos_count_",Sys.Date(),".csv")))


##### model with absolute BAL WBC count 
bal.wbc<-p.eos$BAL_WBC_log
batch<-factor(p.eos$Batch, levels=unique(batch.info.BAL$Batch))
mm4<-model.matrix(~bal.wbc + batch)

# filter using cutoff 0 
drop2<-dropCutoff(0)
y4<-voom(x3[-drop2,], mm4, plot=T)

fit<-lmFit(y4,mm4)
tmp<-contrasts.fit(fit, coef=2)
tmp<-eBayes(tmp)
top.table<-topTable(tmp,coef=1,sort.by="p", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
# write.csv(top.table,file.path(getwd(),"output/",paste0("toptable_bal_WBC_",Sys.Date(),".csv")))

########################
#DEG with DEseq2
########################
library(DESeq2)

#### analyze based on log transformed BAL Eos %
dds<-DESeqDataSetFromMatrix(countData = x,colData=p.eos, design=~BAL_Eos_perc_log+Batch)
dds<-DESeq(dds)
res<-results(dds, name=c("BAL_Eos_perc_log"))
res <- res[order(res$padj),]
head(res)

# signficant results with padj<0.05
res.sig<-res[which(res$padj<0.05),]
res.sig
# write.csv(res.sig,file.path(getwd(),"output/",paste0("Deseq2_balEosPerc_",Sys.Date(),".csv")))

par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Batch") 

#### analyze based on log transformed BAL Eos Count
dds2<-DESeqDataSetFromMatrix(countData = x,colData=p.eos, design=~bal_Eos_ct_log+Batch)
dds2<-DESeq(dds2)
res2<-results(dds2, name=c("bal_Eos_ct_log"))
res2 <- res2[order(res2$padj),]
head(res2)

# signficant results with padj<0.05
res.sig2<-res2[which(res2$padj<0.05),]
res.sig2

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_balEos_ct_",Sys.Date(),".csv")))


par(mfrow=c(1,1))
# Make a basic volcano plot
with(res2, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res2, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res2, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata2 <- vst(dds2, blind=FALSE)
plotPCA(vsdata2, intgroup="Batch") 

#### analyze based on log transformed BAL WBC Count
dds3<-DESeqDataSetFromMatrix(countData = x,colData=p.eos, design=~BAL_WBC_log+Batch)
dds3<-DESeq(dds3)
res3<-results(dds3, name=c("BAL_WBC_log"))
res3 <- res3[order(res3$padj),]
head(res3)

# signficant results with padj<0.05
res.sig3<-res3[which(res3$padj<0.05),]
res.sig3

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_bal_WBC_ct_",Sys.Date(),".csv")))


par(mfrow=c(1,2))
# Make a basic volcano plot
with(res3, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res3, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res3, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata3 <- vst(dds3, blind=FALSE)
plotPCA(vsdata3, intgroup="Batch") 


#### analyze based on serum Eos Count
p.eos<-p.eos%>%mutate(serum_Eos_log=log(serum_Eos+0.001))
hasna<-which(p.eos$serum_Eos_log%>%is.na)
nona<-p.eos[-hasna,"SampleID"]
p.eos.rmna<-p.eos[-hasna,]
x.rmna<-x[,nona]

dds4<-DESeqDataSetFromMatrix(countData = x.rmna,colData=p.eos.rmna, design=~serum_Eos_log+Batch)
dds4<-DESeq(dds4)
res4<-results(dds4, name=c("serum_Eos_log"))
res4 <- res4[order(res4$padj),]
head(res4)

# signficant results with padj<0.05
res.sig4<-res4[which(res4$padj<0.05),]
res.sig4
res.sig4%>%nrow

# write.csv(res.sig2,file.path(getwd(),"output/",paste0("Deseq2_serEos_ct_",Sys.Date(),".csv")))


par(mfrow=c(1,2))

# Make a basic volcano plot
with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata4 <- vst(dds4, blind=FALSE)
plotPCA(vsdata4, intgroup="Batch") 


