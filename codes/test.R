#RNAseq analysis of nasal/bronchial data with edgeR and limma
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
# if working on HPC, set WD as YJ's asthma_transcriptome folder
# if working on local PC, set WD as YJ's asthma-allergy-bunyavanich folder on windows
ifelse(getwd()=="/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome",
       setwd("/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome"),
       setwd("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich"))

ifelse(getwd()=="/sc/arion/projects/asthma-allergy/MS_asthma/Young-Jin/asthma_transcriptome",WDlocation<-"Minerva",WDlocation<-"localPC")

# setting directory for count data when working on local PC
file.pc<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/batch1234_readcount_matrix_allsamples.afterQC.txt")

file.b5.pc<-c("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/input/original_data/MS_asthma/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")

# setting directory for count data when working on HPC minerva
file.minerva<-c("/sc/arion/projects/asthma-allergy/MS_asthma/outputs/batch1234_combined/read_count/batch1234_readcount_matrix_allsamples.afterQC.txt")
file.b5.minerva<-c("/sc/arion/projects/asthma-allergy/MS_asthma/outputs/batch5/read_count/AsthmaNasal_RNAseq.batch5.geneID_readcount_rmscaffold_highqualitysamples.txt")

# reading count data for batch1,2,3,4 and saving as 'counts'
file<-ifelse(WDlocation=="Minerva",file.minerva,file.pc)
file.exists(file)
counts<-read.delim(file[1])
counts.ID<-colnames(counts)

# reading count data for batch 5 as 'counts.b5'
file.b5<-ifelse(WDlocation=="Minerva",file.b5.minerva,file.b5.pc)
file.exists(file.b5)
counts.b5<-read.delim(file.b5[1])
counds.b5ID<-colnames(counts.b5)


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


