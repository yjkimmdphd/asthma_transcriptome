#########################
# cleaning up patient phenotype data table for nasal_bronchial_BAL_rnaseq.R
#########################

library(dplyr)


# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-file.path(getwd(),"input/asthma-phenotype-filtered-revised-2023-09-25.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

# get batch information
filename3<-file.path(getwd(),"input/original_data/MS_asthma/MS_asthma_phenotype.batch1234.txt")
file.exists(filename3)
batch.info<-read.delim(filename3)

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
# find which cell count data is zero and nonzero 
phenotype<-mutate(phenotype,
       zero.BALEos=ifelse(BAL_Eos_perc>0,"nonZero","zero"),
       zero.BALWbc=ifelse(BAL_WBC>0,"nonZero","zero"),
       zero.BALNeut=ifelse(BAL_neut_perc>0,"nonZero","zero"),
       zero.serEos=ifelse(serum_Eos_perc>0,"nonZero","zero"),
       zero.serNeut=ifelse(serum_neut_perc>0,"nonZero","zero"))
z<-lapply(phenotype[,c("zero.BALEos","zero.BALWbc", "zero.BALNeut", "zero.serEos", "zero.serNeut")],
       function(d){factor(d,levels=c("zero","nonZero"))})%>%data.frame()
phenotype[,c("zero.BALEos","zero.BALWbc", "zero.BALNeut", "zero.serEos", "zero.serNeut")]<-z

# subset  'phenotype' with BAL and blood cell count information into 'phenotype.counts'
p.col<-colnames(phenotype)
count.col<-sapply(c("Eos","WBC","serum","BAL","bal","Neut"), function(keyword){grep(keyword,p.col)})%>%unlist%>%unique()
phenotype.counts<-phenotype[,c("subject_assgn",p.col[count.col],"Blood_draw_date")]
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

# calculate time difference between BAL and CBC
p.counts<-p.counts%>%mutate(BAL_Date=BAL_Date%>%as.Date(format="%m/%d/%Y"), 
                            Blood_draw_date=Blood_draw_date%>%as.Date(format="%m/%d/%Y"),
                            BAL_CBC_delay=BAL_Date%>%as.Date(format="%m/%d/%Y")-Blood_draw_date%>%as.Date(format="%m/%d/%Y"))


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

# summary statistics of each cell count
pc.cn<-colnames(p.counts)
pc.cn<-pc.cn[-grep("log",pc.cn)]
pc.cn<-pc.cn[c(3:7,9:13)]
pcrange=sapply(p.counts[,pc.cn],function(d){range(d,na.rm=TRUE)})
pc.cn.df<-data.frame(mean=sapply(p.counts[,pc.cn],function(d){mean(d,na.rm=TRUE)%>%round(2)}),
                  sd=sapply(p.counts[,pc.cn],function(d){sd(d,na.rm=TRUE)%>%round(2)}),
                  range=apply(pcrange,2,function(d){paste(d[1],d[2],sep="-")}))

delay.input<-p.counts
delay.input<-delay.input%>%filter(BAL_CBC_delay<3000)
cbc.bal.delay=data.frame(
  mean.delay=delay.input$BAL_CBC_delay%>%mean(na.rm=TRUE)%>%round(2),
  sd.delay=delay.input$BAL_CBC_delay%>%sd(na.rm=TRUE)%>%round(2),
  range=paste(range(delay.input$BAL_CBC_delay,na.rm=TRUE)[1],
              range(delay.input$BAL_CBC_delay,na.rm=TRUE)[2],
                           sep="-"))
cbc.bal.delay.abs=data.frame(
  mean.delay=delay.input$BAL_CBC_delay%>%abs%>%mean(na.rm=TRUE)%>%round(2),
  sd.delay=delay.input$BAL_CBC_delay%>%abs%>%sd(na.rm=TRUE)%>%round(2),
  range=paste(range(delay.input$BAL_CBC_delay%>%abs,na.rm=TRUE)[1],
              range(delay.input$BAL_CBC_delay%>%abs,na.rm=TRUE)[2],
              sep="-"))

print(pc.cn.df)    
print(cbc.bal.delay)
print(p.counts%>%select(SampleID,BAL_CBC_delay)%>%arrange(desc(BAL_CBC_delay)))
hist(delay.input$BAL_CBC_delay%>%as.numeric, 
     breaks=seq(-1000,4000,100), 
     xlab="BAL date - CBC date (days)", 
     ylab="number of samples")
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

input<-df.deseq2input
z.counts<-c("zero.BALEosnonZero","zero.BALNeutnonZero","zero.serEosnonZero","zero.serNeutnonZero")
z.input<-data.frame(count.data=input[c(11,13,16,17),"count.data"],
           col.data=input[c(11,13,16,17),"col.data"],
           design = paste0("~ ",z.counts," + Batch"),
           resoutput = z.counts)
df.deseq2input<-rbind(input,z.input)
print(df.deseq2input)
