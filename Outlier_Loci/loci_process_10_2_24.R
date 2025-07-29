#Script to process loci output
#The loci output files are MASSIVE, so best to run this on a cluster

library(tidyverse)
library(dplyr)

#Your filename here
loci <- read.csv("MU_loo_assignments_parallel_2_502_loci.csv")

#Reformat columns
loci$locus <- as.factor(loci$locus)
loci$iteration <- as.factor(loci$iteration)

#Remove unecessary column
loci <- loci[-1]

#Generate summary dataframe
#loci_summary<- loci %>%
#  group_by(`locus`) %>%
#  summarize(Freq=length(`locus`),
#            MeanDiff=mean(`Frequency Difference`, na.rm=TRUE),
#            SDDiff=sd(`Frequency Difference`, na.rm=TRUE),
#            MinDiff=min(`Frequency Difference`),
#            MaxDiff=max(`Frequency Difference`))

#Create empty list
summary_list <- list()

# Loop through each column except the first and last
for (i in 2:(ncol(loci) - 1)) {
  # Extract the column name
  
  col_name <- colnames(loci)[i]
  temp <- loci[c("locus", col_name)]%>%
    na.omit(loci[2])
  
  # Create a summary dataframe for the current column
  summary_df <- temp %>%
    group_by(locus) %>%
    summarize(Freq = length(locus),
              MeanDiff = mean(.data[[col_name]], na.rm = TRUE),
              SDDiff = sd(.data[[col_name]], na.rm = TRUE),
              MinDiff = min(.data[[col_name]]),
              MaxDiff = max(.data[[col_name]])) %>%
    # Rename the summary columns with the current column name
    rename_with(~ paste0(col_name, "_", .), -locus)
  
  # Append the summary dataframe to the list
  summary_list[[col_name]] <- summary_df
}

# Merge all the summary dataframes by 'locus'
loci_summary <- reduce(summary_list, full_join, by = "locus")

#Output CSV
write.csv(loci_summary, "MU_loci_2-502_summary.csv")
