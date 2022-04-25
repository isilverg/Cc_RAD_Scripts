
############## Exploratory Data Analyis Part 1- Raw, demultiplexed sequence data ###################

# The point of this analysis is to identify samples that failed to sequence well, any possible patterns related to metadata that are correlated with failed samples, 
# and look at overall distribution of sequence data for your plate and by metadata factors
# Note that you should adjust metadata here for your own dataset. In this example, I look at sample location and RAD plate well, but you could also look at input DNA amount for libraries, year/time of sample
# collection, sex, anything you'd like.
# A few red flags would be noticing patterns of sequence data distribution by RAD plate well- like if whole plate columns or the outer edges of the palte contained failed samples, RA and RB reads should be EQUAL, 
# and if all of your failed samples came from one metadata component- this would be something to consider in subsequent library preps. Essentially, you shouldn't see really ANY patterns in your data here and that
# everything is roughly evenly distributed. It is very common to have at least a few samples fail, but these can still have important data and could be included in your downstream analysis, but important to know.




########### Setup #################
library(here)
list.files()

library(ggplot2)
library(tidyr)
library(data.table)
library(reshape2)

#Multiplot Functions####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##### I. Initial QC checks of raw read data ########

###Data read in, basic clean up:####
datarenamed<-read.table("cc_linecount_bysample.txt", header=F, sep="") #read in the line counts file you generated in Step 3 of workflow, using the raw data
colnames(datarenamed)<-c("No_seq","sample") #rename columns to be number of sequences and sample name
datarenamed$No_seq<-datarenamed$No_seq/4 #line counts in file need to be divided by 4 since each sequence has 4 lines in fastq files- this gives true number of reads
datarenamed$log_No_seqs<-log(datarenamed$No_seq) #get log transformed number of seqs
vars<-colsplit(datarenamed$sample,"_",c("species","well","turtle", "loc","seq_direction")) #separate sample name by metadata factors that are in the name of your sample- CHANGE THIS FOR YOUR DATA
datarenamed<-cbind(vars,datarenamed)#combine seq info with metadata info

samplesonly<-subset(datarenamed,species != "total") # remove line at the end of your file that has total count of all sequences, if you have blank wells you also need to remove these to look at them separately here, so could add !="BLANK" in addition
levels(samplesonly$species<-factor(samplesonly$species))#double check for typos in categories
levels(samplesonly$sample.id<-factor(samplesonly$turtle))#double check for typos in categories
levels(samplesonly$mpatch<-factor(samplesonly$loc))#double check for typos in categories
levels(samplesonly$well<-factor(samplesonly$well))#double check for typos in categories
levels(samplesonly$seq_direction<-factor(samplesonly$seq_direction))#double check for typos in categories



#Summaries by plate:
mean<-tapply(samplesonly$No_seq, samplesonly$species,mean) #if you did not subset to remove empty wells , then you can do this for the object "datarenamed" instead of "samplesonly". I do this by species since this plate only has ONE species so this will include ALL samples in plate
min<-tapply(samplesonly$No_seq, samplesonly$species,min)
max<-tapply(samplesonly$No_seq, samplesonly$species,max)
median<-tapply(samplesonly$No_seq, samplesonly$species,median)
summary_by_plate<-data.frame(mean,median,min,max)

#####Histograms of raw reads by plate, with binning for thresholds####
sample.reads.bin<-samplesonly #from above
sample.reads.binRA<-subset(sample.reads.bin,seq_direction=="RA.fastq") #subset out only forward reads
sample.reads.binRA$Bin_seq <- ifelse(sample.reads.binRA$No_seq <= 10000, 10000, ifelse((sample.reads.binRA$No_seq >10000) & (sample.reads.binRA$No_seq < 100000), sample.reads.binRA$No_seq, 100000)) #bin samples 
ggplot(sample.reads.binRA, aes(x=Bin_seq)) +
  geom_histogram(colour="black", fill="blue",binwidth = 10000) +theme_bw()+
  facet_wrap(~species,scales="free")

sample.reads.binRA$QC_cat <- ifelse(sample.reads.binRA$No_seq <= 10000, "failed", ifelse((sample.reads.binRA$No_seq >10000) & (sample.reads.binRA$No_seq <= 100000), "good", "high")) #add in column for qc metric
sample.reads.binRAshort<-sample.reads.binRA
sample.reads.binRAshort$Bin_seqlow<-ifelse(sample.reads.binRAshort$No_seq <= 20000, sample.reads.binRAshort$No_seq, 20000)#bin everything high to be able to look at low distribution and confirm we think 10K is reasonable cutoff for pass/fail
ggplot(sample.reads.binRAshort, aes(x=Bin_seqlow)) +
  geom_histogram(colour="black", fill="blue",binwidth = 1000) +theme_bw()+
  facet_wrap(~species,scales="free")

failed<-subset(sample.reads.binRAshort, QC_cat=="failed")
table(failed$species)
qc_pass<-subset(sample.reads.binRAshort, QC_cat!="failed")

#Patterns with metadata for failed samples####
#The point here is to look at patterns pertaining to metadata, so these might include collection date, extraction method, etc.
str(failed)
failed$loc<-as.factor(failed$loc)#change text categories to factors for plotting
plot(failed$loc,failed$No_seq,las=2) #plot by metadata component 1 (e.g. location, year)
plot(failed$well,failed$No_seq,las=2) #plot by RAD plate well


#Summaries by metadata components:

#by location:
mean<-tapply(samplesonly$No_seq, samplesonly$loc,mean) #if you did not subset to remove empty wells in line 50, then do this for the object "datarenamed" instead of "samplesonly"
min<-tapply(samplesonly$No_seq, samplesonly$loc,min)
max<-tapply(samplesonly$No_seq, samplesonly$loc,max)
median<-tapply(samplesonly$No_seq, samplesonly$loc,median)
summary_by_loc<-data.frame(mean,median,min,max) #get summary by metadata factor, here using location
samplesonly$loc<-as.factor(samplesonly$loc)#change character strings to factors for plotting
plot(samplesonly$loc,samplesonly$No_seq,las=2) #plot number sequences by location
qc_pass$loc<-as.factor(qc_pass$loc)
plot(qc_pass$loc,qc_pass$No_seq) #plotting number sequences by location with failed samples removed

mean<-tapply(qc_pass$No_seq, qc_pass$loc,mean) #if you did not subset to remove empty wells in line 50, then do this for the object "datarenamed" instead of "samplesonly"
min<-tapply(qc_pass$No_seq, qc_pass$loc,min)
max<-tapply(qc_pass$No_seq, qc_pass$loc,max)
median<-tapply(qc_pass$No_seq, qc_pass$loc,median)
summary_by_loc<-data.frame(mean,median,min,max)

plot(samplesonly$loc,samplesonly$log_No_seqs,las=2,ylab="log # sequences")

#by RAD plate well
plot(samplesonly$well,samplesonly$No_seq,las=2)

