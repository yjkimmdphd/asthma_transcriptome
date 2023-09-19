##
# RNA seq analysis of nasal/bronchial data with edgeR and limma
# YJK local machine version
## 
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)

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

# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/asthma-phenotype-filtered-revised.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

# get batch information
filename3<-file.path(getwd(),"input/original_data/MS_asthma/MS_asthma_phenotype.batch1234.txt")
file.exists(filename3)
batch.info<-read.delim(filename3)

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

#### make dataframe of input values that will be used for DEG with DESeq2
# save as 'df.deseq2input'
# x is original count table
# x2 is x filtered for low gene counts
df.deseq2input<-data.frame(count.data="x",col.data="p.counts", design = paste0("~",p.counts.var,"+ Batch"), resoutput = p.counts.var)
df.deseq2input[c(3,4),1:2]<-data.frame(rep("x.BalNeut",2),rep("p.count.BalNeut",2))
df.deseq2input[6:10,1:2]<-data.frame(rep("x.SerCt",5),rep("p.count.SerCt",5))
df.input2<-df.deseq2input%>%mutate(count.data=sub("x","x2",count.data))
df.deseq2input<-rbind(df.deseq2input,df.input2)
print(df.deseq2input)

#### filtering phenotype table based on cell counts
p.BalEos.pos<-p.counts%>%filter(bal_Eos_ct>0)
p.BalNeut.pos<-p.counts%>%filter(BAL_neut_ct>0)
p.serEos.pos<-p.counts%>%filter(serum_Eos>0)
p.serNeut.pos<-p.counts%>%filter(serum_Neut>0)

#### make new count tables with sampleID filtered for cell counts
x2.BalEos.pos<-x2[,p.BalEos.pos$SampleID]
x2.BalNeut.pos<-x2[,p.BalNeut.pos$SampleID]
x2.serEos.pos<-x2[,p.serEos.pos$SampleID]
x2.serNeut.pos<-x2[,p.serNeut.pos$SampleID]

# update 'df.deseq2input'
df.deseq2input[21:28,"count.data"]<-rep(c("x2.BalEos.pos","x2.BalNeut.pos","x2.serEos.pos","x2.serNeut.pos"),each=2)
df.deseq2input[21:28,"col.data"]<-rep(c("p.BalEos.pos","p.BalNeut.pos","p.serEos.pos","p.serNeut.pos"),each=2)
df.deseq2input[21:28,c("design","resoutput")]<-df.deseq2input[c(1:4,6:9),c("design","resoutput")]
print(df.deseq2input)


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
for(i in 1:nrow(df.deseq2input)){
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}

for(i in 13:nrow(df.deseq2input)){
  assign(paste0("res",i),deseq2DEG(df.deseq2input[i,1],df.deseq2input[i,2],df.deseq2input[i,3],df.deseq2input[i,4]))
}

# summarize number of signficant genes in the 'res.table'
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





#######################
#Make a basic volcano plot

with(res4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(res4, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

