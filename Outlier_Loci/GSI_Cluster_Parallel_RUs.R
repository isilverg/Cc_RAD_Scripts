#Outlier Loci GSI Optimization

#Load Packages
library(tidyverse)
library(dplyr)
library(dartR)
library(poppr)
library(vcfR)
library(adegenet)
library(ggplot2)
library(parallel)

#Bring in filtered data

#Read in VCF file
in_vcf <- read.vcfR("../Combined_filtered_normed_norels_locname_newnames_sorted.vcf", convertNA = T)

#Convert vcf to genlight. Save twice: genlight.orig will be preserved if samples are subset and the genlight object will be replaced; if no subsetting is needed, then the file named genlight will be used throughout the rest of the script. 
genlight.orig = genlight <- vcfR2genlight(in_vcf)

#Read in the population map. This should be a two-column tab-delimited text file where the first column contains the name of each individual sample and the second column contains the population. This can be the same as the population map used in Stacks.
popmap.orig <- read_tsv("../Combined_sort_RUs_popmap.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

# Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

genlight_2 <- gl.compliance.check(genlight)

#genind <- gl2gi(genlight_2) # need genind format for some functions
#summary(genlight_2@pop)
#summary(genind@pop)

# Create a list of individuals for each population
#pop_names <- unique(genind$pop)
#indiv_names <- genlight_2@other[["ind.metrics"]][["id"]]
#populations <- lapply(pop_names, function(pop) indiv_names[pop == genind$pop])

# Find all combinations in serial, then run LOO & loci optimization in parallel
# Initialize an empty list to store combinations
#combinations_list <- list()

#Nested loops to iterate through each combination
# Nested loops to iterate through each combination
#for (indiv_pop1 in populations[[1]]) {
#   for (indiv_pop2 in populations[[2]]) {
#    for (indiv_pop3 in populations[[3]]) {
#      for (indiv_pop4 in populations[[4]]) {
# Create a list for the current combination
#        current_combination <- list(
#          population1 = indiv_pop1,
#          population2 = indiv_pop2,
#          population3 = indiv_pop3,
#          population4 = indiv_pop4
#        )
# Append the current combination to the list of combinations
#        combinations_list <- c(combinations_list, list(current_combination))
#      }
#    }
#  }
#}

#Create fx
process_combination <- function(combination, genind, top_seq) {
  individual_names <- unlist(combination)
  
  # Subset geninds
  training_genind <- genind[!row.names(genind@tab) %in% individual_names]
  validation_genind <- genind[individual_names]
  
  ## ID outlier loci
  genpop <- genind2genpop(training_genind, quiet = TRUE) 
  # Calculate allele frequencies
  freq <- makefreq(genpop, missing = 0, quiet = TRUE)
  
  # Calculate allele frequency differences for all pairs of populations
  freqdiffs <- list()
  for (j in 1:(nrow(freq) - 1)) {
    for (k in (j + 1):nrow(freq)) {
      freqdiffs[[paste(row.names(freq)[j], row.names(freq)[k], sep = '-')]] <- abs(freq[j, ] - freq[k, ])
    }
  }
  
  # Convert to data frame
  as.data.frame(do.call(rbind, freqdiffs))
  
  #Initialize empty list to store results
  datalist<-list()
  
  #Loop through loci
  for (l in seq_along(top_seq)){
    
    n.top<-top_seq[[l]]
    
    # Combination 1: PFL - NRU
    NRU_PFL <- as.data.frame(freqdiffs[["NRU-PFL"]])
    colnames(NRU_PFL) <- "NRUPFL"
    NRU_PFL$locus <- gsub("\\..*", "", rownames(NRU_PFL))
    NRU_PFL <- NRU_PFL %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_PFL) <- NRU_PFL$locus
    NRUPFL_top <- NRU_PFL %>% top_n(n.top, NRUPFL)
    
    # Combination 2: NRU - DRTO
    NRU_DRTO <- as.data.frame(freqdiffs[["DRTO-NRU"]])
    colnames(NRU_DRTO) <- "NRUDRTO"
    NRU_DRTO$locus <- gsub("\\..*", "", rownames(NRU_DRTO))
    NRU_DRTO <- NRU_DRTO %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_DRTO) <- NRU_DRTO$locus
    NRUDRTO_top <- NRU_DRTO %>% top_n(n.top, NRUDRTO)
    
    # Combination 3: NRU - NGM
    NRU_NGM <- as.data.frame(freqdiffs[["NGM-NRU"]])
    colnames(NRU_NGM) <- "NRUNGM"
    NRU_NGM$locus <- gsub("\\..*", "", rownames(NRU_NGM))
    NRU_NGM <- NRU_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_NGM) <- NRU_NGM$locus
    NRUNGM_top <- NRU_NGM %>% top_n(n.top, NRUNGM)
    
    # Combination 4: DRTO - NGM
    DRTO_NGM <- as.data.frame(freqdiffs[["DRTO-NGM"]])
    colnames(DRTO_NGM) <- "DRTONGM"
    DRTO_NGM$locus <- gsub("\\..*", "", rownames(DRTO_NGM))
    DRTO_NGM <- DRTO_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(DRTO_NGM) <- DRTO_NGM$locus
    DRTONGM_top <- DRTO_NGM %>% top_n(n.top, DRTONGM)
    
    # Combination 5: DRTO - PFL
    DRTO_PFL <- as.data.frame(freqdiffs[["DRTO-PFL"]])
    colnames(DRTO_PFL) <- "DRTOPFL"
    DRTO_PFL$locus <- gsub("\\..*", "", rownames(DRTO_PFL))
    DRTO_PFL <- DRTO_PFL %>% distinct(locus, .keep_all = TRUE)
    rownames(DRTO_PFL) <- DRTO_PFL$locus
    DRTOPFL_top <- DRTO_PFL %>% top_n(n.top, DRTOPFL)
    
    # Combination 6: PFL - NGM
    PFL_NGM <- as.data.frame(freqdiffs[["NGM-PFL"]])
    colnames(PFL_NGM) <- "PFLNGM"
    PFL_NGM$locus <- gsub("\\..*", "", rownames(PFL_NGM))
    PFL_NGM <- PFL_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(PFL_NGM) <- PFL_NGM$locus
    PFLNGM_top <- PFL_NGM %>% top_n(n.top, PFLNGM)
    
    diff_loci <- NRUPFL_top %>%
      full_join(NRUDRTO_top, by = "locus") %>%
      full_join(NRUNGM_top, by = "locus") %>%
      full_join(DRTONGM_top, by = "locus") %>%
      full_join(DRTOPFL_top, by = "locus") %>%
      full_join(PFLNGM_top, by = "locus")
    diff_loci<-diff_loci[c(2,1,3:7)]
    
    #Isolate/format outlier locus names for data subsetting
    outlier_loci <- diff_loci %>% dplyr::select(locus) %>% 
      distinct(locus)
    outlier_loci<- as.matrix(outlier_loci)
    
    genlight_out <-gl.keep.loc(genlight_2, loc.list = outlier_loci, verbose=0)
    genind_out <- gl2gi(genlight_out, verbose=0)
    
    training_genind<- genind_out[!row.names(genind@tab) %in% individual_names]
    validation_genind<-genind_out[individual_names]
    
    #Run DAPC
    n.da <- 6 #Set number of DAs
    dapc <- dapc(training_genind, training_genind@pop, n.pca=60, n.da = n.da, var.contrib = TRUE)
    
    set.seed(5)
    temp_score <- optim.a.score(dapc, plot=FALSE) #Determine best number of PCs to retain for DAPC analysis
    n.pc <- temp_score$best
    dapc2 <-dapc(training_genind,training_genind@pop, n.pca = n.pc, n.da = n.da)
    
    dapc_pred<-predict.dapc(dapc2, newdata=validation_genind)
    
    ##Make DF for validation assignments
    val_assignments<-as.data.frame(dapc_pred$posterior)
    val_assignments$ind<-rownames(val_assignments)
    
    for (m in 1:length(val_assignments$ind)){
      if (str_sub(val_assignments$ind[m], 1,2) == "NR"){    
        val_assignments$pop[m] <- "NR"
      } 
      else {
        if (str_sub(val_assignments$ind[m], 1,2) == "CE" || str_sub(val_assignments$ind[m], 1,2) == "SE"|| str_sub(val_assignments$ind[m], 1,2) == "SW" || str_sub(val_assignments$ind[m], 1,2) == "CW"){
          val_assignments$pop[m] <- "PF"
        }
        else {
          if (str_sub(val_assignments$ind[m], 1,2) == "DR"){
            val_assignments$pop[m] <- "DR"
          }
          else {
            if (str_sub(val_assignments$ind[m], 1,2) == "NG"){
              val_assignments$pop[m] <- "NG"
            }
          }
        }
      }
    }
    
    val_assignments$post<-str_sub(colnames(val_assignments[,1:4])[max.col(val_assignments[,1:4],ties.method="first")], 1,2)
    val_assignments$loci<-top_seq[[l]]
    
    datalist[[l]]<-val_assignments
  }
  return(datalist)
}

#Read in combinations list
load("../Combos/RU_combinations_list.RData")
tmp <- bind_rows(combinations_list)

NRUinds <- unique(tmp$population1)
PFLinds <- unique(tmp$population2)
DRTOinds <- unique(tmp$population3)
NGMinds <- unique(tmp$population4)

final_inds <- c(NRUinds, PFLinds, DRTOinds, NGMinds)

genlight_2 <- gl.keep.ind(genlight_2, ind.list = final_inds)

#Make genind
genind <- gl2gi(genlight_2)

#Check populations
summary(genlight_2@pop)
summary(genind@pop)

#Create a list of individuals for each population
pop_names <- unique(genind$pop)
indiv_names <- genlight_2@other[["ind.metrics"]][["id"]]
populations <- lapply(pop_names, function(pop) indiv_names[pop == genind$pop])

#Set parameters
n<-1000
rand_combs <- sample(combinations_list, n, replace = FALSE)
datalist <- vector("list", length = length(rand_combs))
top_seq <- seq(from = 10, to = 10000, by = 100)
loci_list <- vector("list", length = length(top_seq))
cores <- 8

# Apply the process_combination function to each element of rand_combs using mclapply
results <- mclapply(rand_combs, process_combination, genind=genind, top_seq=top_seq, mc.cores=cores)

# Format dataframe, visualize results, and export
loo_assignments<-dplyr::bind_rows(results)
loo_assignments$loci<-as.factor(loo_assignments$loci)
loo_assignments$correct<-ifelse(loo_assignments$pop==loo_assignments$post, "1", "0")
write.csv(loo_assignments, "RU_loo_assignments.csv")
