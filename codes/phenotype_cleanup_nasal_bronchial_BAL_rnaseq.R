#########################
# cleaning up patient phenotype data table for nasal_bronchial_BAL_rnaseq.R
#########################

library(dplyr)

# calculate absolute Eos from total BAL WBC and BAL Eos %
phenotype<-phenotype%>%mutate(
  bal_Eos_ct=BAL_Eos_perc*0.01*BAL_WBC,
  BAL_neut_ct=BAL_neut_perc*0.01*BAL_WBC)


# Subject 279's serum Eos % was NA while serum WBC and serum Eos absolute counts were available. A code is to calculate the  missing serum Eos % for subject 279.  
r.279<-which(phenotype$Study.ID=="279")
p.279<-phenotype%>%filter(Study.ID=="279")%>%mutate(serum_Eos_perc=round(serum_Eos/serum_WBC,4)*100)
phenotype[r.279,]<-p.279

# log transform serum and BAL cell counts
phenotype<-mutate(phenotype, 
       BAL_Eos_ct_log=log(bal_Eos_ct+0.001)%>%round(2),
       BAL_Eos_perc_log = log(BAL_Eos_perc+0.001)%>%round(2),
       BAL_WBC_log=log(BAL_WBC+0.00001)%>%round(2),
       BAL_neut_perc_log=log(BAL_neut_perc+0.001)%>%round(2),
       BAL_neut_ct_log10=log10(BAL_neut_ct+0.00001)%>%round(2),
       serum_Eos_log10=log10(serum_Eos+0.001)%>%round(2),
       serum_Eos_perc_log=log(serum_Eos_perc+0.001)%>%round(2),
       serum_Neut_log10=log10(serum_Neut+0.001)%>%round(2),
       serum_Neut_perc_log=log(serum_neut_perc+0.001)%>%round(2),
       serum_WBC_log10=log10(serum_WBC)%>%round(2)
       )

# subset  'phenotype' with BAL and blood cell count information into 'phenotype.counts'
p.col<-colnames(phenotype)
count.col<-sapply(c("Eos","WBC","serum","BAL","bal"), function(keyword){grep(keyword,p.col)})%>%unlist%>%unique()
phenotype.counts<-phenotype[,c("subject_assgn",p.col[count.col])]
phenotype.counts<-phenotype.counts[,c("subject_assgn",colnames(phenotype.counts)[-1]%>%sort)]
colnames(phenotype.counts)[1]<-"SampleID"




# find subject assignment ID with nasal and bronchial cell RNAseq data  
nasal.ID<-phenotype[,"subject_assgn"] # nasal study ID that have been collected
bronch.ID<-sub("N","B",nasal.ID) # bronchial subject ID 
nb.ID<-c(nasal.ID,bronch.ID) # nasal and bronchial subject ID
a<-counts.ID%in%nb.ID  
nb.exist.ID<-counts.ID[a] # subject ID nasal and bronchial samples that are in the count matrix

# subset 'phenotype.counts' for which RNA counts data exist to 'p.counts'
sample.ID<-counts.ID[counts.ID%in%nasal.ID]
p.counts<-phenotype.counts[nasal.ID%in%sample.ID,]

# subset batch info filtered based on nb.exist.ID.
batch.info.ID<-which(batch.info$SampleID%in%nb.exist.ID)
batch.info.BAL<-batch.info[batch.info.ID,]
(batch.info.BAL$Type=="Nasal")%>%table

# left join batch info with 'p.counts'
p.counts<-left_join(p.counts,batch.info.BAL, by="SampleID")
p.counts$Batch<-factor(p.counts$Batch, levels=unique(batch.info.BAL$Batch))
p.counts$Race_corrected<-factor(p.counts$Race_corrected, levels=p.counts$Race_corrected%>%unique)

par(mfrow=c(4,6))

for(i in 3:22){
  p.counts[,i]%>%hist(main=colnames(p.counts)[i])
}

# some serum or Eos counts are NA of missing values
# make phenotype and count table removing NA
# could conveniently use na.omit(), but that would remove too many data points and reduce power 
a<-sapply(p.counts,is.na)
a1<- unique(which(a==TRUE)%%45)
p.counts[a1,]
p.count.BalNeut<-p.counts[-c(41,42),]
p.count.SerCt<-p.counts[-c(13,22,28,37,38),]

sapply(p.counts,is.na)%>%colSums() 
sapply(p.count.BalNeut,is.na)%>%colSums() # check if BAL cols with Neut of the new dataframe has 'NA'
sapply(p.count.SerCt,is.na)%>%colSums() # check if cols with serum of the new dataframe has 'NA'

# should use 'p.count.BalNeut' for DEG using BAL Neut Ct or %
# should use 'p.count.SerCt' for DEG using Ser counts 


