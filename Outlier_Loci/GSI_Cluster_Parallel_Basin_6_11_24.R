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
in_vcf <- read.vcfR("./data/Combined_filtered_normed_norels_locname_newnames_sorted.vcf", convertNA = T)

#Convert vcf to genlight. Save twice: genlight.orig will be preserved if samples are subset and the genlight object will be replaced; if no subsetting is needed, then the file named genlight will be used throughout the rest of the script. 
genlight.orig = genlight <- vcfR2genlight(in_vcf)

#Read in the population map. This should be a two-column tab-delimited text file where the first column contains the name of each individual sample and the second column contains the population. This can be the same as the population map used in Stacks.
popmap.orig <- read_tsv("./data/Combined_sort_basins_popmap.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

# Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

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

# Nested loops to iterate through each combination
for (indiv_pop1 in populations[[1]]) {
  for (indiv_pop2 in populations[[2]]) {
# Create a list for the current combination
      current_combination <- list(
        population1 = indiv_pop1,
        population2 = indiv_pop2
      )
# Append the current combination to the list of combinations
     combinations_list <- c(combinations_list, list(current_combination))
  }
}

#Create fx
process_combination <- function(combination, genind, top_seq) {
  individual_names <-  unlist(combination)
  
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
  locilist <- list()
  
  #Loop through loci
  for (l in seq_along(top_seq)){

    n.top<-top_seq[[l]]
    
    # Combination 1: ATL - GOM
    ATL_GOM <- as.data.frame(freqdiffs[["ATL-GOM"]])
    colnames(ATL_GOM) <- "ATLGOM"
    ATL_GOM$locus <- gsub("\\..*", "", rownames(ATL_GOM))
    ATL_GOM <- ATL_GOM %>% distinct(locus, .keep_all = TRUE)
    rownames(ATL_GOM) <- ATL_GOM$locus
    ATLGOM_top <- ATL_GOM %>% top_n(n.top, ATLGOM)
    
    diff_loci <- ATLGOM_top 
    diff_loci<-diff_loci[c(2,1)]
    diff_loci$iteration<-top_seq[[l]]
    
    # Store diff_loci in locilist
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
    val_assignments$ind<-rownames(val_assignments)
    
    for (m in 1:length(val_assignments$ind)){
      if (str_sub(val_assignments$ind[m], 1,2) == "NR" || str_sub(val_assignments$ind[m], 1,2) == "CE" || str_sub(val_assignments$ind[m], 1,2) == "SE"){
        val_assignments$pop[m] <- "ATL"
      } 
      else { 
        if (str_sub(val_assignments$ind[m], 1,2) == "DR" || str_sub(val_assignments$ind[m], 1,2) == "SW" || str_sub(val_assignments$ind[m], 1,2) == "CW" || str_sub(val_assignments$ind[m], 1,2) == "NG"){
          val_assignments$pop[m] <- "GOM"
        }
      }
    }
    #Has been giving issues, but can fix in R locally
    val_assignments$post<-str_sub(colnames(val_assignments[,2:3])[max.col(val_assignments[,2:3],ties.method="first")], 1,2)
    val_assignments$loci<-top_seq[[l]]
    
    datalist[[l]]<-val_assignments
  }
  return(list(datalist = datalist, locilist = locilist))
}

#Read in the combination list
load("../Combos/Basin_combinations_list.RData")

#Set parameters
n<-5382
rand_combs <- sample(combinations_list, n, replace = FALSE)
datalist <- vector("list", length = length(rand_combs))
top_seq <- seq(from = 2, to = 1002, by = 25)
locilist <- vector("list", length = length(top_seq))
cores <- 8

# Apply the process_combination function to each element of rand_combs using mclapply
results <- mclapply(rand_combs, process_combination, genind=genind, top_seq=top_seq, mc.cores=cores)

# Format dataframe, visualize results, and export
datalist<-list()
for (i in 1:length(results)){
  datalist[[i]]<-results[[1]][["datalist"]]
}

loo_assignments<-dplyr::bind_rows(datalist)
loo_assignments$loci<-as.factor(loo_assignments$loci)
loo_assignments$correct<-ifelse(loo_assignments$pop==loo_assignments$post, "1", "0")
write.csv(loo_assignments, "Basin_loo_assignments_parallel_1-1001.csv")

locilist<-list()
for (i in 1:length(results)){
  locilist[[i]]<-results[[1]][["locilist"]]
}

loci_list<-dplyr::bind_rows(locilist)
write.csv(loci_list, "Basin_loo_assignments_1-1001_loci.csv")

