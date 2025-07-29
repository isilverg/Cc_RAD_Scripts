##LOAD PACKAGES

library(tidyverse)
library(dplyr)
library(poppr)
library(dartR) #this calls SNPrelate as a dependency, which seems to have been taken off CRAN
library(vcfR)
library(adegenet)
library(Matrix)

## LOAD DATA

#Read VCF
in_vcf <- read.vcfR("Combined_filtered_normed_norels_locname_newnames_sorted.vcf", convertNA = T)

#Convert vcf to genlight. Save twice: genlight.orig will be preserved if samples are subset and the genlight object will be replaced; if no subsetting is needed, then the file named genlight will be used throughout the rest of the script. 
genlight.orig = genlight <- vcfR2genlight(in_vcf)

#Read in the population map. This should be a two-column tab-delimited text file where the first column contains the name of each individual sample and the second column contains the population. This can be the same as the population map used in Stacks.
popmap.orig <- read_tsv("Combined_sort_popmap.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

#Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

#Run compliance check to ensure genlight is functional
genlight_2 <- gl.compliance.check(genlight)

#Make genind
genind <- gl2gi(genlight_2) # need genind format for some functions

#Check populations
summary(genlight_2@pop)
summary(genind@pop)

# Create a list of individuals for each population
pop_names <- unique(genind$pop) 
indiv_names <- genlight_2@other[["ind.metrics"]][["id"]]
populations <- lapply(pop_names, function(pop) indiv_names[pop == genind$pop])

##LOO GSI

#First, identify all combinations of individuals

#Initialize an empty list to store combinations
combinations_list <- list()

#Nested loops to iterate through each combination
for (indiv_pop1 in populations[[1]]) {
  for (indiv_pop2 in populations[[2]]) {
    for (indiv_pop3 in populations[[3]]) {
      for (indiv_pop4 in populations[[4]]) {
        for (indiv_pop5 in populations[[5]]) {
          for (indiv_pop6 in populations[[6]]) {
            for (indiv_pop7 in populations[[7]]) {
              # Create a list for the current combination
              current_combination <- list(
                population1 = indiv_pop1,
                population2 = indiv_pop2,
                population3 = indiv_pop3,
                population4 = indiv_pop4,
                population5 = indiv_pop5,
                population6 = indiv_pop6,
                population7 = indiv_pop7
              )
              # Append the current combination to the list of combinations
              combinations_list <- c(combinations_list, list(current_combination))
            }
          }
        }
      }
    }
  }
}

save(combinations_list, file = "MU_combinations_list.RData")
