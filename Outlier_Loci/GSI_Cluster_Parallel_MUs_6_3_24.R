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
popmap.orig <- read_tsv("../Combined_sort_popmap.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

# Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

genlight_2 <- gl.compliance.check(genlight)

genind <- gl2gi(genlight_2) # need genind format for some functions
summary(genlight_2@pop)
summary(genind@pop)

# Create a list of individuals for each population
pop_names <- unique(genind$pop)
indiv_names <- genlight_2@other[["ind.metrics"]][["id"]]
populations <- lapply(pop_names, function(pop) indiv_names[pop == genind$pop])

# Find all combinations in serial, then run LOO & loci optimization in parallel
# Initialize an empty list to store combinations
#combinations_list <- list()

# Nested loops to iterate through each combination
#for (indiv_pop1 in populations[[1]]) {
#  for (indiv_pop2 in populations[[2]]) {
#    for (indiv_pop3 in populations[[3]]) {
#      for (indiv_pop4 in populations[[4]]) {
#        for (indiv_pop5 in populations[[5]]) {
#          for (indiv_pop6 in populations[[6]]) {
#            for (indiv_pop7 in populations[[7]]) {
              # Create a list for the current combination
#              current_combination <- list(
#                population1 = indiv_pop1,
#                population2 = indiv_pop2,
#                population3 = indiv_pop3,
#                population4 = indiv_pop4,
#                population4 = indiv_pop5,
#                population4 = indiv_pop6,
#                population4 = indiv_pop7
#              )
              # Append the current combination to the list of combinations
#              combinations_list <- c(combinations_list, list(current_combination))
#            }
#          }
#        }
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
  locilist<-list()
  
  #Loop through loci
  for (l in seq_along(top_seq)){
    
    n.top<-top_seq[[l]]
    
    # Combination 1: CEFL - NRU
    NRU_CEFL <- as.data.frame(freqdiffs[["CEFL-NRU"]])
    colnames(NRU_CEFL) <- "NRUCEFL"
    NRU_CEFL$locus <- gsub("\\..*", "", rownames(NRU_CEFL))
    NRU_CEFL <- NRU_CEFL %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_CEFL) <- NRU_CEFL$locus
    NRUCEFL_top <- NRU_CEFL %>% top_n(n.top, NRUCEFL)
    
    # Combination 2: NRU - SEFL
    NRU_SEFL <- as.data.frame(freqdiffs[["NRU-SEFL"]])
    colnames(NRU_SEFL) <- "NRUSEFL"
    NRU_SEFL$locus <- gsub("\\..*", "", rownames(NRU_SEFL))
    NRU_SEFL <- NRU_SEFL %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_SEFL) <- NRU_SEFL$locus
    NRUSEFL_top <- NRU_SEFL %>% top_n(n.top, NRUSEFL)
    
    # Combination 3: NRU - DRTO
    NRU_DRTO <- as.data.frame(freqdiffs[["DRTO-NRU"]])
    colnames(NRU_DRTO) <- "NRUDRTO"
    NRU_DRTO$locus <- gsub("\\..*", "", rownames(NRU_DRTO))
    NRU_DRTO <- NRU_DRTO %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_DRTO) <- NRU_DRTO$locus
    NRUDRTO_top <- NRU_DRTO %>% top_n(n.top, NRUDRTO)
    
    # Combination 4: NRU - SWFL
    NRU_SWFL <- as.data.frame(freqdiffs[["NRU-SWFL"]])
    colnames(NRU_SWFL) <- "NRUSWFL"
    NRU_SWFL$locus <- gsub("\\..*", "", rownames(NRU_SWFL))
    NRU_SWFL <- NRU_SWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_SWFL) <- NRU_SWFL$locus
    NRUSWFL_top <- NRU_SWFL %>% top_n(n.top, NRUSWFL)
    
    # Combination 5: NRU - CWFL
    NRU_CWFL <- as.data.frame(freqdiffs[["CWFL-NRU"]])
    colnames(NRU_CWFL) <- "NRUCWFL"
    NRU_CWFL$locus <- gsub("\\..*", "", rownames(NRU_CWFL))
    NRU_CWFL <- NRU_CWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_CWFL) <- NRU_CWFL$locus
    NRUCWFL_top <- NRU_CWFL %>% top_n(n.top, NRUCWFL)
    
    # Combination 6: NRU - NGM
    NRU_NGM <- as.data.frame(freqdiffs[["NGM-NRU"]])
    colnames(NRU_NGM) <- "NRUNGM"
    NRU_NGM$locus <- gsub("\\..*", "", rownames(NRU_NGM))
    NRU_NGM <- NRU_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(NRU_NGM) <- NRU_NGM$locus
    NRUNGM_top <- NRU_NGM %>% top_n(n.top, NRUNGM)
    
    # Combination 7: CEFL - SEFL
    CEFL_SEFL <- as.data.frame(freqdiffs[["CEFL-SEFL"]])
    colnames(CEFL_SEFL) <- "CEFLSEFL"
    CEFL_SEFL$locus <- gsub("\\..*", "", rownames(CEFL_SEFL))
    CEFL_SEFL <- CEFL_SEFL %>% distinct(locus, .keep_all = TRUE)
    rownames(CEFL_SEFL) <- CEFL_SEFL$locus
    CEFLSEFL_top <- CEFL_SEFL %>% top_n(n.top, CEFLSEFL)
    
    # Combination 8: CEFL - DRTO
    CEFL_DRTO <- as.data.frame(freqdiffs[["CEFL-DRTO"]])
    colnames(CEFL_DRTO) <- "CEFLDRTO"
    CEFL_DRTO$locus <- gsub("\\..*", "", rownames(CEFL_DRTO))
    CEFL_DRTO <- CEFL_DRTO %>% distinct(locus, .keep_all = TRUE)
    rownames(CEFL_DRTO) <- CEFL_DRTO$locus
    CEFLDRTO_top <- CEFL_DRTO %>% top_n(n.top, CEFLDRTO)
    
    # Combination 9: CEFL - SWFL
    CEFL_SWFL <- as.data.frame(freqdiffs[["CEFL-SWFL"]])
    colnames(CEFL_SWFL) <- "CEFLSWFL"
    CEFL_SWFL$locus <- gsub("\\..*", "", rownames(CEFL_SWFL))
    CEFL_SWFL <- CEFL_SWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(CEFL_SWFL) <- CEFL_SWFL$locus
    CEFLSWFL_top <- CEFL_SWFL %>% top_n(n.top, CEFLSWFL)
    
    # Combination 10: CEFL - CWFL
    CEFL_CWFL <- as.data.frame(freqdiffs[["CEFL-CWFL"]])
    colnames(CEFL_CWFL) <- "CEFLCWFL"
    CEFL_CWFL$locus <- gsub("\\..*", "", rownames(CEFL_CWFL))
    CEFL_CWFL <- CEFL_CWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(CEFL_CWFL) <- CEFL_CWFL$locus
    CEFLCWFL_top <- CEFL_CWFL %>% top_n(n.top, CEFLCWFL)
    
    # Combination 11: CEFL - NGM
    CEFL_NGM <- as.data.frame(freqdiffs[["CEFL-NGM"]])
    colnames(CEFL_NGM) <- "CEFLNGM"
    CEFL_NGM$locus <- gsub("\\..*", "", rownames(CEFL_NGM))
    CEFL_NGM <- CEFL_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(CEFL_NGM) <- CEFL_NGM$locus
    CEFLNGM_top <- CEFL_NGM %>% top_n(n.top, CEFLNGM)
    
    # Combination 12: SEFL - DRTO
    SEFL_DRTO <- as.data.frame(freqdiffs[["DRTO-SEFL"]])
    colnames(SEFL_DRTO) <- "SEFLDRTO"
    SEFL_DRTO$locus <- gsub("\\..*", "", rownames(SEFL_DRTO))
    SEFL_DRTO <- SEFL_DRTO %>% distinct(locus, .keep_all = TRUE)
    rownames(SEFL_DRTO) <- SEFL_DRTO$locus
    SEFLDRTO_top <- SEFL_DRTO %>% top_n(n.top, SEFLDRTO)
    
    # Combination 13: SEFL - SWFL
    SEFL_SWFL <- as.data.frame(freqdiffs[["SEFL-SWFL"]])
    colnames(SEFL_SWFL) <- "SEFLSWFL"
    SEFL_SWFL$locus <- gsub("\\..*", "", rownames(SEFL_SWFL))
    SEFL_SWFL <- SEFL_SWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(SEFL_SWFL) <- SEFL_SWFL$locus
    SEFLSWFL_top <- SEFL_SWFL %>% top_n(n.top, SEFLSWFL)
    
    # Combination 14: SEFL - CWFL
    SEFL_CWFL <- as.data.frame(freqdiffs[["CWFL-SEFL"]])
    colnames(SEFL_CWFL) <- "SEFLCWFL"
    SEFL_CWFL$locus <- gsub("\\..*", "", rownames(SEFL_CWFL))
    SEFL_CWFL <- SEFL_CWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(SEFL_CWFL) <- SEFL_CWFL$locus
    SEFLCWFL_top <- SEFL_CWFL %>% top_n(n.top, SEFLCWFL)
    
    # Combination 15: SEFL - NGM
    SEFL_NGM <- as.data.frame(freqdiffs[["NGM-SEFL"]])
    colnames(SEFL_NGM) <- "SEFLNGM"
    SEFL_NGM$locus <- gsub("\\..*", "", rownames(SEFL_NGM))
    SEFL_NGM <- SEFL_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(SEFL_NGM) <- SEFL_NGM$locus
    SEFLNGM_top <- SEFL_NGM %>% top_n(n.top, SEFLNGM)
    
    # Combination 16: DRTO - SWFL
    DRTO_SWFL <- as.data.frame(freqdiffs[["DRTO-SWFL"]])
    colnames(DRTO_SWFL) <- "DRTOSWFL"
    DRTO_SWFL$locus <- gsub("\\..*", "", rownames(DRTO_SWFL))
    DRTO_SWFL <- DRTO_SWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(DRTO_SWFL) <- DRTO_SWFL$locus
    DRTOSWFL_top <- DRTO_SWFL %>% top_n(n.top, DRTOSWFL)
    
    # Combination 17: DRTO - CWFL
    DRTO_CWFL <- as.data.frame(freqdiffs[["CWFL-DRTO"]])
    colnames(DRTO_CWFL) <- "DRTOCWFL"
    DRTO_CWFL$locus <- gsub("\\..*", "", rownames(DRTO_CWFL))
    DRTO_CWFL <- DRTO_CWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(DRTO_CWFL) <- DRTO_CWFL$locus
    DRTOCWFL_top <- DRTO_CWFL %>% top_n(n.top, DRTOCWFL)
    
    # Combination 18: DRTO - NGM
    DRTO_NGM <- as.data.frame(freqdiffs[["DRTO-NGM"]])
    colnames(DRTO_NGM) <- "DRTONGM"
    DRTO_NGM$locus <- gsub("\\..*", "", rownames(DRTO_NGM))
    DRTO_NGM <- DRTO_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(DRTO_NGM) <- DRTO_NGM$locus
    DRTONGM_top <- DRTO_NGM %>% top_n(n.top, DRTONGM)
    
    # Combination 19: SWFL - CWFL
    SWFL_CWFL <- as.data.frame(freqdiffs[["CWFL-SWFL"]])
    colnames(SWFL_CWFL) <- "SWFLCWFL"
    SWFL_CWFL$locus <- gsub("\\..*", "", rownames(SWFL_CWFL))
    SWFL_CWFL <- SWFL_CWFL %>% distinct(locus, .keep_all = TRUE)
    rownames(SWFL_CWFL) <- SWFL_CWFL$locus
    SWFLCWFL_top <- SWFL_CWFL %>% top_n(n.top, SWFLCWFL)
    
    # Combination 20: SWFL - NGM
    SWFL_NGM <- as.data.frame(freqdiffs[["NGM-SWFL"]])
    colnames(SWFL_NGM) <- "SWFLNGM"
    SWFL_NGM$locus <- gsub("\\..*", "", rownames(SWFL_NGM))
    SWFL_NGM <- SWFL_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(SWFL_NGM) <- SWFL_NGM$locus
    SWFLNGM_top <- SWFL_NGM %>% top_n(n.top, SWFLNGM)
    
    # Combination 21: CWFL - NGM
    CWFL_NGM <- as.data.frame(freqdiffs[["CWFL-NGM"]])
    colnames(CWFL_NGM) <- "CWFLNGM"
    CWFL_NGM$locus <- gsub("\\..*", "", rownames(CWFL_NGM))
    CWFL_NGM <- CWFL_NGM %>% distinct(locus, .keep_all = TRUE)
    rownames(CWFL_NGM) <- CWFL_NGM$locus
    CWFLNGM_top <- CWFL_NGM %>% top_n(n.top, CWFLNGM)
    
    diff_loci <- NRUCEFL_top%>%
      full_join(NRUSEFL_top, by = "locus") %>%
      full_join(NRUDRTO_top, by = "locus") %>%
      full_join(NRUSWFL_top, by = "locus") %>%
      full_join(NRUCWFL_top, by = "locus") %>%
      full_join(NRUNGM_top, by = "locus") %>%
      full_join(CEFLSEFL_top, by = "locus") %>%
      full_join(CEFLDRTO_top, by = "locus") %>%
      full_join(CEFLSWFL_top, by = "locus") %>%
      full_join(CEFLCWFL_top, by = "locus") %>%
      full_join(CEFLNGM_top, by = "locus") %>%
      full_join(SEFLDRTO_top, by = "locus") %>%
      full_join(SEFLSWFL_top, by = "locus") %>%
      full_join(SEFLCWFL_top, by = "locus") %>%
      full_join(SEFLNGM_top, by = "locus") %>%
      full_join(DRTOSWFL_top, by = "locus") %>%
      full_join(DRTOCWFL_top, by = "locus") %>%
      full_join(DRTONGM_top, by = "locus") %>%
      full_join(SWFLCWFL_top, by = "locus") %>%
      full_join(SWFLNGM_top, by = "locus") %>%
      full_join(CWFLNGM_top, by = "locus")
    diff_loci<-diff_loci[c(2,1,3:22)]
    
    #Add iteration column
    diff_loci$iteration<-top_seq[[l]]
    
    #Store diff_loci in locilist
    locilist[[l]] <- diff_loci

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
    val_assignments$pop<-str_sub(row.names(dapc_pred$ind.scores), 1,2) #assignments$pop<-gsub ("_", "",   assignments$pop)
    val_assignments$ind<-rownames(val_assignments)
    val_assignments$post<-str_sub(colnames(val_assignments[,1:4])[max.col(val_assignments[,1:4],ties.method="first")], 1,2)
    val_assignments$loci<-top_seq[[l]]
    
    datalist[[l]]<-val_assignments
    
  }
  return(list(datalist = datalist, locilist = locilist))
}

#Read in combinations list
load("../Combos/MU_combinations_list.RData")

#Set parameters
n<-1000
rand_combs <- sample(combinations_list, n, replace = FALSE)
datalist <- vector("list", length = length(rand_combs))
top_seq <- seq(from = 2, to = 502, by = 10)
loci_list <- vector("list", length = length(top_seq))
cores <- 56

# Apply the process_combination function to each element of rand_combs using mclapply
results <- mclapply(rand_combs, process_combination, genind=genind, top_seq=top_seq, mc.cores=cores)

# Format dataframe, visualize results, and export
datalist<-list()
for (i in 1:length(results)){
	datalist[[i]]<-results[[i]][["datalist"]]
}

loo_assignments<-dplyr::bind_rows(datalist)
loo_assignments$loci<-as.factor(loo_assignments$loci)
loo_assignments$correct<-ifelse(loo_assignments$pop==loo_assignments$post, "1", "0")
write.csv(loo_assignments, "MU_loo_assignments_parallel_2_502.csv")

locilist<-list()
for (i in 1:length(results)){
	locilist[[i]]<-results[[i]][["locilist"]]
}

loci_list<-dplyr::bind_rows(locilist)
write.csv(loci_list, "MU_loo_assignments_parallel_2_502_loci.csv")
