################ Exploratory Data Analysis Part 2 - simplified version #############################


# In part 2 of EDA, we are looking at mapping statistics and how the reads compare before and after removing
# PCR duplicates
# Inputs for this script are the All_head_flagstat_reformat.txt file generated immediately after the mapping step, and
# All_head_filt_flagstat_reformat.txt generated after duplicate removal step
# From this script, you should generate a csv file with statistics reporting percent reads mapped before and after removing PCR dups
# You should also produce histograms of mapping reads and proportion of duplicates in reads

#setup

library(here) #this should set your working directory to your project folder where this script is located
list.files()

library(ggplot2)
library(tidyr)
library(data.table)
library(reshape2)


##### Read in files and reformat

list.files(path="./",pattern="cc_All_head") #list files with this pattern in the title
mapped.raw<-read.delim("./cc_All_head_flagstat_reformat.txt",header=T) #read in file generated after mapping
mapped.fltr<-read.delim("./cc_All_head_filt_flagstat_reformat.txt",header=T) #read in file generated after removing PCR dups
mapped.raw$ID<-gsub("_sort_flagstat.txt", "", mapped.raw$sampleID)
mapped.fltr$ID<-gsub("_sortfltr_flagstat.txt", "", mapped.fltr$sampleID)
mapped.comb<-merge(mapped.raw,mapped.fltr, by="ID") #make one dataframe with both mapped read types
mapped.comb<-mapped.comb[,c(1,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24,13,25)] #reorganize dataframe, this can be changed if needed, but this column format generally works
vars<-colsplit(mapped.comb$ID, "_", c("LABID","well","turtle","mpatch")) #the options after c() will be specifiy to your naming scheme for your samples, CHANGE FOR YOUR DATA
mapped.comb<-cbind(mapped.comb,vars)
mapped.comb<-mapped.comb[,c(1,24:27,2:23)] #tidy up


#get proportions of filtered to unfiltered mapped reads
mapped.comb$prop.NDreads<-mapped.comb$mapped.reads.QC.passed.y/mapped.comb$total.reads.QC.passed.x
summary(mapped.comb$prop.NDreads)#see percentage of raw reads left after filtering
sd(mapped.comb$prop.NDreads) #get sd of how many reads are dropped to filtering
hist(mapped.comb$prop.NDreads) #plot histogram of prop of nonduplicated reads
mapped.comb$prop.dupreads<-(1-mapped.comb$prop.NDreads) #get proportion of reads that were duplicates
summary(mapped.comb$prop.dupreads)
sd(mapped.comb$prop.dupreads)
hist(mapped.comb$prop.dupreads)

mapped.comb.short<-mapped.comb[,c(1:11,29)] #create shorter version of he above dataframe
str(mapped.comb.short)
names(mapped.comb.short) <- c("ID","LABID","well", "turtle", "mpatch", "UF_total_reads_QC_passed"	,"F_total_reads_QC_passed","UF_mapped_reads_QC_passed","F_mapped_reads_QC_passed",	"UF_percent_reads_mapped","F_percent_reads_mapped","prop") #you should change some of these if the variables aren't in your dataset, like mpatch
hist(mapped.comb.short$F_percent_reads_mapped)
write.csv(mapped.comb.short,"../short_mapping_stats.csv") #change path to where you want your mapping stat file written to

######################add in QC category info from section I- OPTIONAL####################
# you can do this step if you included your failed sampels during your mapping step, if not this can be skipped. Basically this just allows us to look at EDA with only qc-passed samples


#SRA1<-sample.reads.binRAshort[,c(1,5:9)] #this dataframe was generated in EDA part 1 script, which is upstream of this one. Can copy from that to get qc info.
#SRA1$platewell<-paste(SRA1$plate,SRA1$well)
#mapped.comb.short$platewell<-paste(mapped.comb.short$plate,mapped.comb.short$well)
#CFS<-merge(mapped.comb.short,SRA1,by="platewell")
#CFS.ok<-subset(CFS,QC_cat!="failed")



