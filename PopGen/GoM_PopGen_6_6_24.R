##GoM PopGen Script##

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

#devtools::install_github("thierrygosselin/radiator")
##Load Packages##
library(devtools)
library(tidyverse)
library(dartR) #this calls SNPrelate as a dependency, which seems to have been taken off CRAN
library(poppr)
library(vcfR)
library(adegenet)
library(hierfstat)
#library(pophelper)#not on CRAN
library(gridExtra)
library(ggpubr)
#library(strataG)#not on CRAN
library(mapdata)
library(viridis)
library(wesanderson)
library(ape)
library(reshape2)
#library(stackoverflow)
library(Matrix)
library(sf)
library(ggplot2)
library(ggspatial)
library(plotly)
library(MetBrewer)
library(patchwork)
library(matrixStats)
library(cowplot)
library(grid)
library(gridExtra)
#library(rgdal)
library(circlize)
library(MASS)
#library(OutFLANK)
library(giscoR)
library(radiator)
library(rubias)
library(svglite)
library(gridExtra)
library(pophelper)

#setwd("C:/Users/Ian/Desktop/Dissertation/Chapter 4 Genetics/Genomics Metadata and Data/VCFTools_Out")
##Load Raw Data and Filter Using Whitelists from Previous Steps##
#Read in VCF file
in_vcf <- read.vcfR("./data/Combined_filtered_normed_norels_locname_newnames_sorted.vcf", convertNA = T)

#Convert vcf to genlight. Save twice: genlight.orig will be preserved if samples are subset and the genlight object will be replaced; if no subsetting is needed, then the file named genlight will be used throughout the rest of the script. 
genlight.orig = genlight <- vcfR2genlight(in_vcf)

#Read in the population map. This should be a two-column tab-delimited text file where the first column contains the name of each individual sample and the second column contains the population. This can be the same as the population map used in Stacks.
popmap.orig <- read_tsv("./data/Combined_sort_popmap.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

# Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

# Read in a "regions" file. This file should have at a minimum a column for Population, latitude, and longitude. We assume this is a csv (comma-separate) file.
region <- read.csv("ccatl_regions.csv", header = T)

#Tibbles will only show a couple sig digits for lat and long, but the full coordinates are stored.
latlon <- as.data.frame(cbind(region$Population,
                              region$lat,
                              region$lon)) %>% 
  dplyr::rename(STRATA = V1, lat = V2, lon = V3) %>% 
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))

#Replicate the coordinates for each sample
all_latlon <- popmap.filtered %>% left_join(latlon, by = "STRATA") %>% 
  dplyr::select(-Indiv)

latlon <- all_latlon[,2:3]
genlight@other$latlon<-latlon

genlight_2 <- gl.compliance.check(genlight)

#-------------------Subset the full set of samples for analysis------------------------
#Only run this section if you want to run the analysis on a subset of individuals included in the filtered VCF file. If you want to run on all the samples in the vcf file that was loaded, then skip this section and move on to "Finish creating genlight and genind objects" a few lines below.
#Read in white list of individuals/locs. This should be a text file with a single column containing the names of the individuals that were retained after filtering.
whitelist.indvs <- read.delim("./Whitelist_indivs_CCATL.txt", header = F, strip.white = T) %>% pull()
genlight_3 <- gl.keep.ind(genlight_2, ind.list = whitelist.indvs)
genlight_3@ind.names

whitelist.locs <- read.delim("./Whitelist_loci_CCATL.txt", header = F, strip.white = T) %>%
  dplyr::rename(Chrom = V1, Pos = V2) %>%
  arrange(Chrom, Pos)
whitelist.locs_2 <- str_c(whitelist.locs$Chrom, "_", whitelist.locs$Pos) #ISG: WHYYYYYY
genlight_4 <-gl.keep.loc(genlight_3, loc.list = whitelist.locs_2, verbose=0)
rm(genlight, genlight_2, genlight_3)

#-------------------Finish creating genlight and genind objects--------------------------
genind <- gl2gi(genlight_2) # need genind format for some functions
summary(genlight_2@pop)
summary(genind@pop)

##Some Basic Stats##
dartR_het<-gl.report.heterozygosity(genlight_2)
dartR_div<-gl.report.diversity(genlight_2)
#gl2genepop(genlight_2, outpath = "C:/Users/Ian/Desktop/Dissertation/Chapter 4 Genetics/Genomics Metadata and Data/VCFTools_Out")

fst_boot<-gl.fst.pop(genlight_2, nboots=10000, percent=95)
fst_boot_fsts<-data.matrix(round(fst_boot[["Fsts"]], 3))
fst_boot_ps<-data.matrix(round(fst_boot[["Pvalues"]], 10))
fst_boot_fsts
fst_boot_ps

##PCAs and Plots##
pca1 <- glPca(genlight_2,center = T, scale = T, nf = 5)
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

# Extract PC scores to color by location/yr:
# Adapted from https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$Population <- pop(genlight_2)
pca1.scores$reg <- pca1.scores$pop

#Make PCA plots
set.seed(9)
num_pops <- length(levels(factor(popmap.filtered$STRATA)))
pop_cols <- c("#000000","#E69F00" , "#56B4E9","#009E73")
pop_shapes<-c(15,16,17,18)

# plot PC 1 and 2
PCA1_2<-ggscatter(pca1.scores, x = "PC2", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)")) +
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  #scale_y_continuous(limits=c(-110,80), breaks=c(-75, -50, -25, 0, 25, 50),
                     #labels=c("-75", "-50", "-25", "0", "25", "50"))+
  theme(panel.grid = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_2

# #plot  PC 1 and 3
PCA1_3<-ggscatter(pca1.scores, x = "PC3", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"))+
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_3

pcas<-grid.arrange(PCA1_2, PCA1_3, ncol=2, nrow=1)
ggsave(file="Figure1_PCAs.svg", plot=pcas, width=10, height=5)

##Extract Legend
PCA_Legend<-ggscatter(pca1.scores, x = "PC2", y = "PC1", shape = "Population", color = "Population",
                      palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                      xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                      ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)")) +
  scale_shape_manual(values=pop_shapes) +
  theme(legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20)) +
  #scale_y_continuous(limits=c(-110,80), breaks=c(-75, -50, -25, 0, 25, 50),
  #labels=c("-75", "-50", "-25", "0", "25", "50"))+
  theme(panel.grid = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
dev.off()
legend<-ggpubr::get_legend(PCA_Legend)
legend_plot<-as_ggplot(legend)
ggsave(file="Figure1_PCAs_legend.svg", plot=legend_plot, width=7, height=0.5) 


##Initial DAPC
n_individuals <- nrow(popmap.filtered)
n_pops <- length(levels(factor(popmap.filtered$STRATA)))

grp_all <- find.clusters(genind, max.n.clust=10, n.pca=200,
                         choose.n.clust = FALSE)

BIC<-as.data.frame(cbind(seq(1,10,1), grp_all$Kstat))

ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K)") + font("xlab", face = "bold") + font("ylab", face = "bold")
grp_all$Kstat

k <- 7
grp_all <- find.clusters(genind, max.n.clust=n_pops, n.pca=200, n.clust = k)
grp_all$size

n.da <- 6 #Set number of DAs
dapc <- dapc(genind, grp_all$grp, n.pca=60, n.da = n.da, var.contrib = TRUE)
pal <- get_palette("Accent", k = k)
dpca_result <- scatter(dapc, col=pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)
set.seed(4)
#contrib<-loadingplot(dapc$var.contr,axis=2, thres=.002,lab.jitter=1)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",c(1:k)),
          ncol=1, xlab="Individuals", col=pal, lab=genind@pop, show.lab = T)

dapc$IND <- row.names(dapc$ind.coord)

dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))

colnames(dapc_info) <- c("IND","DPC1", "DPC2", "DPC3", "K")
colnames(dapc_info) <- c("IND", c(paste0("DPC", 1:(n.da-3))), "K")
dapc_info$SITE <- str_sub(dapc_info$IND, 1,2) #Takes first two letters of the Indv name and adds to the columns "SITE"
dapc_info

table.value(table(pop(genind), grp_all$grp), col.lab=paste("inf", 1:4),
            row.lab=paste("ori", 1:4))

#Second DAPC
set.seed(5)
dapc_a_score <- dapc(genind, n.pca=100, n.da=n.da) 
temp_score <- optim.a.score(dapc_a_score) #Determine best number of PCs to retain for DAPC analysis
n.pc <- temp_score$best #Save best number of PCs
dapc2 <-dapc(genind,genind@pop, n.pca = n.pc, n.da = n.da)
dapc2$IND <- row.names(dapc2$ind.coord)
dapc2_info <- as.data.frame(cbind(dapc2$IND, dapc2$ind.coord, grp_all$grp))
dapc2_info$IND <- row.names(dapc2_info)

colnames(dapc2_info) <- c("IND", "LD1", "LD2", "LD3", "K") 
dapc2_info$SITE <- gsub("[^a-zA-Z]","",dapc2_info$IND)
dapc2_info$SITE<-as.factor(dapc2_info$SITE)
dapc2_info

jpeg(filename="dapc_scatter.jpeg", width=10, height=8, units="in", res = 300)
scatter(dapc2,
        scree.pca=TRUE,
        pch = pop_shapes,
        cstar = 0,
        cellipse = 1,
        solid = 0.6,
        cex = 2,
        clab = 1,
        ratio.da = 0.15,
        ratio.pca = 0.15,
        col = pop_cols)
dev.off()

ggscatter(dapc2_info, x = "DPC1", y = "DPC2", shape = 21, fill = "SITE", size = 3,
          xlab = "DA1", ylab = "DA2", palette = pop_cols)


dapc2_info%>% 
  ggplot(aes(x=LD2, y=LD1, color=SITE), group=SITE) +
  geom_point(size = 4, alpha=0.5)+
  stat_ellipse()+
  scale_color_manual(values = pop_cols)+
  scale_shape_manual(values=pop_shapes)
  theme_classic()+
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_blank())


load_dpca2 <- as.data.frame(dapc2$var.contr)
write.table(load_dpca2, "Loadings_SNPs.txt", sep="\t", row.names=FALSE, quote=FALSE)
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, main = "Percent of genetic variance explained by eigenvectors", ylim = c(0, max(percent) + 10), names.arg = paste0("E", c(1:length(percent))), ylab = "Percent", xlab = "Eigenvector")
dapc_prior=as.data.frame(dapc2$ind.coord)
write.table(dapc_prior, "DAPC_prior_results.txt", quote=F, sep="\t", row.names=TRUE)

### Add information to the tab results of the DPCA
dapc_prior$IND <- row.names(dapc_prior)

### Add site info
dapc_prior$SITE <- str_sub(dapc_prior$IND, 1,2) #Again, assumes the site is the first two characters of the individual names
#dapc_prior$SITE <- gsub ('_', '', dapc_prior$SITE)

### Add 'region' info:
#dapc_prior <- merge(dapc_prior, region, by.x = "SITE", by.y = "Population")

### Make a ggplot graph representing the DAPC for the first and second axes for the regions
ggscatter(dapc_prior, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 3,
          xlab = "DA1", ylab = "DA2", palette = pop_cols)
summary(dapc2)

compoplot(dapc2)

assignments<-as.data.frame(dapc2$posterior)
assignments$ind<-row.names(assignments)
assignments$pop<-str_sub(assignments$ind, 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
assignments$pop <- "BLANK"
for (i in 1:length(assignments$ind)){
  if (str_sub(assignments$ind[i], 1,2) == "NR" || str_sub(assignments$ind[i], 1,2) == "CE" || str_sub(assignments$ind[i], 1,2) == "SE"){
    assignments$pop[i] <- "ATL"
  } 
  else { 
    if (str_sub(assignments$ind[i], 1,2) == "DR" || str_sub(assignments$ind[i], 1,2) == "SW" || str_sub(assignments$ind[i], 1,2) == "CW" || str_sub(assignments$ind[i], 1,2) == "NG"){
      assignments$pop[i] <- "GOM"
    }
  }
}

for (i in 1:length(assignments$ind)){
  if (str_sub(assignments$ind[i], 1,2) == "NR"){    
    assignments$pop[i] <- "NRU"
  } 
  else {
    if (str_sub(assignments$ind[i], 1,2) == "CE" || str_sub(assignments$ind[i], 1,2) == "SE"|| str_sub(assignments$ind[i], 1,2) == "SW" || str_sub(assignments$ind[i], 1,2) == "CW"){
    assignments$pop[i] <- "PFL"
    }
    else {
      if (str_sub(assignments$ind[i], 1,2) == "DR"){
      assignments$pop[i] <- "DRTO"
      }
      else {
        if (str_sub(assignments$ind[i], 1,2) == "NG"){
        assignments$pop[i] <- "NGM"
        }
      }
  }
  }
}

assignments$post<-colnames(assignments[,1:7])[max.col(assignments[,1:7],ties.method="first")] #change based on number pops
ass_melt<-melt(assignments)

ass_melt<-ass_melt%>%
  arrange(pop) %>%
  arrange(ind)
ass_ind_table <- table(ass_melt$ind)
ind_levels <- names(ass_ind_table)[order(ass_ind_table)]
ass_melt$ind <- factor(ass_melt$ind, levels = ind_levels)

ass_melt$ind <- factor(ass_melt$ind, levels = popmap.filtered$Indiv) #order based on sorted ind list

compoplot_dapc2<-ggplot(ass_melt, 
       aes(fill=variable, y=value, x=ind),
       color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols)+
  theme(legend.position="none") +
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank())

ggsave(file="MU_Barplot_DAPC.jpeg", plot=compoplot_dapc2, width=10, height=5) 

compoplot_dapc2_legend<-ggplot(ass_melt, 
                               aes(fill=variable, y=value, x=ind),
                               color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols, name= "Population")+
  theme(legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20))+
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank())

dev.off()
legend_compoplot<-ggpubr::get_legend(compoplot_dapc2_legend)
legend_comp_plot<-as_ggplot(legend_compoplot)
ggsave(file="Figure2_DAPC_legend.svg", plot=legend_comp_plot) 


##fast Structure BarPlot##
sfiles<-list.files(path="./faststructure_q", full.names = TRUE)
slist<-readQ(files=sfiles)
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)
plotQ(slist[5],
               barbordercolour="black",
               showsp = FALSE,
                splab = "",
               showindlab=F,
               indlabangle=45,
                clustercol = c("lightgrey", "darkgrey", "brown", "red"),
               returnplot=F,
               exportplot=T,
               outputfilename="plotq_4",
               height = 5,
               width = 10,
               dpi = 500,
               units = "in",
               imgtype="tiff",
               exportpath=getwd())





#################################################################################################
#####-----GSI------##############################################################################
#################################################################################################
genpop<-genind2genpop(genind)
freq<-makefreq(genpop, missing=0)

freqdiffs<-list()
for (i in 1:(nrow(freq) - 1)) {
  for (j in (i + 1):nrow(freq)) {
    freqdiffs[[paste(row.names(freq)[i], row.names(freq)[j], sep = '-')]] <- abs(freq[i, ] - freq[j, ])
  }
}
as.data.frame(do.call(rbind, freqdiffs))

CWFL_DRTO <- as.data.frame(freqdiffs[["CWFL-DRTO"]])
colnames(CWFL_DRTO)<-"CWFLDRTO"
CWFL_DRTO$locus <- gsub("\\..*", "", rownames(CWFL_DRTO))
CWFL_DRTO<-CWFL_DRTO %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_DRTO)<-CWFL_DRTO$locus
CWFLDRTO_top<-CWFL_DRTO %>% top_n(40, CWFLDRTO)
#CWFLDRTO_top<- CWFLDRTO_top %>% dplyr::select("CWFLDRTO")

CWFL_NGM <- as.data.frame(freqdiffs[["CWFL-NGM"]])
colnames(CWFL_NGM)<-"CWFLNGM"
CWFL_NGM$locus <- gsub("\\..*", "", rownames(CWFL_NGM))
CWFL_NGM<-CWFL_NGM %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_NGM)<-CWFL_NGM$locus
CWFLNGM_top<-CWFL_NGM %>% top_n(40, CWFLNGM)
#CWFLNGM_top<- CWFLNGM_top %>% dplyr::select("CWFLNGM")

CWFL_SWFL <- as.data.frame(freqdiffs[["CWFL-SWFL"]])
colnames(CWFL_SWFL)<-"CWFLSWFL"
CWFL_SWFL$locus <- gsub("\\..*", "", rownames(CWFL_SWFL))
CWFL_SWFL<-CWFL_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_SWFL)<-CWFL_SWFL$locus
CWFLSWFL_top<-CWFL_SWFL %>% top_n(40, CWFLSWFL)
#CWFLSWFL_top<- CWFLSWFL_top %>% dplyr::select("CWFLSWFL")

DRTO_NGM <- as.data.frame(freqdiffs[["DRTO-NGM"]])
colnames(DRTO_NGM)<-"DRTONGM"
DRTO_NGM$locus <- gsub("\\..*", "", rownames(DRTO_NGM))
DRTO_NGM<-DRTO_NGM %>% distinct(locus, .keep_all = TRUE)
rownames(DRTO_NGM)<-DRTO_NGM$locus
DRTONGM_top<-DRTO_NGM %>% top_n(40, DRTONGM)
#DRTONGM_top<- DRTONGM_top %>% dplyr::select("DRTONGM")

DRTO_SWFL <- as.data.frame(freqdiffs[["DRTO-SWFL"]])
colnames(DRTO_SWFL)<-"DRTOSWFL"
DRTO_SWFL$locus <- gsub("\\..*", "", rownames(DRTO_SWFL))
DRTO_SWFL<-DRTO_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(DRTO_SWFL)<-DRTO_SWFL$locus
DRTOSWFL_top<-DRTO_SWFL %>% top_n(40, DRTOSWFL)
#DRTOSWFL_top<- DRTOSWFL_top %>% dplyr::select("DRTOSWFL")

NGM_SWFL <- as.data.frame(freqdiffs[["NGM-SWFL"]])
colnames(NGM_SWFL)<-"NGMSWFL"
NGM_SWFL$locus <- gsub("\\..*", "", rownames(NGM_SWFL))
NGM_SWFL<-NGM_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(NGM_SWFL)<-NGM_SWFL$locus
NGMSWFL_top<-NGM_SWFL %>% top_n(40, NGMSWFL)
#NGMSWFL_top<- NGMSWFL_top %>% dplyr::select("NGMSWFL")

diff_loci<- CWFLDRTO_top %>% full_join(CWFLNGM_top, by="locus") %>% 
  full_join(CWFLSWFL_top, by="locus") %>%
  full_join(DRTONGM_top, by="locus") %>%
  full_join(DRTOSWFL_top, by="locus") %>%
  full_join(NGMSWFL_top, by="locus")

diff_loci<-diff_loci[c(2,1,3,4,5,6,7)]

outlier_loci <- diff_loci %>% dplyr::select(locus) %>% 
  distinct(locus)
outlier_loci<- as.matrix(outlier_loci)
genlight_out <-gl.keep.loc(genlight_2, loc.list = outlier_loci, verbose=0)
genind_out <- gl2gi(genlight_out)

genpop_out<-genind2genpop(genind_out)
freq_out<-as.data.frame(t(makefreq(genpop_out, missing=0)))

##Check for loci out of HWE : No loci out of HWE
gl.filter.hwe(genlight_out, subset="each", n.pop.threshold = 3, alpha_val = 0.1)

##PCAs and Plots##
pca1 <- glPca(genlight_out,center = T, scale = T, nf = 5)
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

# Extract PC scores to color by location/yr:
# Adapted from https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$Population <- pop(genlight_out)
pca1.scores$reg <- pca1.scores$pop

#Make PCA plots
set.seed(9)
num_pops <- length(levels(factor(popmap.filtered$STRATA)))
pop_cols <- c("#000000","#E69F00" , "#56B4E9"  ,"#009E73")
pop_shapes<-c(15,16,17,18)

# plot PC 1 and 2
PCA1_2<-ggscatter(pca1.scores, x = "PC2", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)")) +
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  #scale_y_continuous(limits=c(-110,80), breaks=c(-75, -50, -25, 0, 25, 50),
  #labels=c("-75", "-50", "-25", "0", "25", "50"))+
  theme(panel.grid = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_2

# #plot  PC 1 and 3
PCA1_3<-ggscatter(pca1.scores, x = "PC3", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"))+
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_3

pcas<-grid.arrange(PCA1_2, PCA1_3, ncol=2, nrow=1)
ggsave(file="Figure##__out_PCAs.svg", plot=pcas, width=10, height=5)

##Extract Legend
PCA_Legend<-ggscatter(pca1.scores, x = "PC2", y = "PC1", shape = "Population", color = "Population",
                      palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                      xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                      ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)")) +
  scale_shape_manual(values=pop_shapes) +
  theme(legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20)) +
  #scale_y_continuous(limits=c(-110,80), breaks=c(-75, -50, -25, 0, 25, 50),
  #labels=c("-75", "-50", "-25", "0", "25", "50"))+
  theme(panel.grid = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
dev.off()
legend<-ggpubr::get_legend(PCA_Legend)
legend_plot<-as_ggplot(legend)
ggsave(file="Figure1_PCAs_legend.svg", plot=legend_plot, width=7, height=0.5) 

##LDA
#Create random list of training and validation individuals, using half of each population
ctemp<-sample(rep(c("T", "V"), round((0.5*length(which(popmap.filtered$STRATA == "CWFL"))), 0)))
dtemp<-sample(rep(c("T", "V"), round((0.5*length(which(popmap.filtered$STRATA == "DRTO"))), 0)))
ntemp<-sample(rep(c("T", "V"), round((0.5*length(which(popmap.filtered$STRATA == "NGM"))), 0)))
stemp<-sample(rep(c("T", "V"), (0.5*length(which(popmap.filtered$STRATA == "SWFL"))))) #no round because odd number
extra<-sample (c("T","V"), 1) #add because SWFL odd
tv_list<-c(ctemp, dtemp, ntemp,stemp, extra)

pca1.scores$lda_status<-c(tv_list)
pca1.scores$lda_status<-as.factor(pca1.scores$lda_status)
training<-which(pca1.scores$lda_status == "T")
validation<-which(pca1.scores$lda_status == "V")
pca1.scores.lda<-pca1.scores[,-7]

z <- lda(Population~.,pca1.scores.lda, subset=training)
z1<-predict(z,pca1.scores.lda[training,])
z1$class==pca1.scores.lda[training,]$Population
z1_lds<-as.data.frame(z1$x)
z1_lds$lda_status<-c("Training")
z1_lds$orig<-gsub('[[:digit:]]+',"",rownames(z1_lds))
z1_lds$post<-z1$class

out <- predict(z,pca1.scores.lda[validation,])
out$class==pca1.scores.lda[validation,]$Population
out_lds<-as.data.frame(out$x)
out_lds$lda_status<-c("Validation")
out_lds$orig<-gsub('[[:digit:]]+',"",rownames(out_lds))
out_lds$post<-out$class

LD1_2<-ggplot() +
  geom_point(data=z1_lds, aes(x=LD2, y=LD1, color=orig), shape = 15, size = 4, alpha=0.5)+
  stat_ellipse(aes(x=z1_lds$LD2, y=z1_lds$LD1, color=z1_lds$orig), type="norm", linetype=2)+
  geom_point(data=out_lds, aes(x=LD2, y=LD1, color=orig), shape = 16, size = 4)+
  stat_ellipse(aes(x=out_lds$LD2, y=out_lds$LD1, color=out_lds$orig), type="norm")+
  scale_color_manual(values = pop_cols)+
  theme_classic()+
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
LD1_2

LD1_3<-ggplot() +
  geom_point(data=z1_lds, aes(x=LD3, y=LD1, color=orig), shape = 15, size = 4, alpha=0.5)+
  stat_ellipse(aes(x=z1_lds$LD3, y=z1_lds$LD1, color=z1_lds$orig), type="norm", linetype=2)+
  geom_point(data=out_lds, aes(x=LD3, y=LD1, color=orig), shape = 16, size = 4)+
  stat_ellipse(aes(x=out_lds$LD3, y=out_lds$LD1, color=out_lds$orig), type="norm")+
  scale_color_manual(values = pop_cols)+
  theme_classic()+
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
LD1_3

LD_legend<-ggplot() +
  geom_point(data=z1_lds, aes(x=LD2, y=LD1, color=orig), shape = 15, size = 4, alpha=0.5)+
  stat_ellipse(aes(x=z1_lds$LD2, y=z1_lds$LD1, color=z1_lds$orig), type="norm", linetype=2)+
  geom_point(data=out_lds, aes(x=LD2, y=LD1, color=orig), shape = 16, size = 4)+
  stat_ellipse(aes(x=out_lds$LD2, y=out_lds$LD1, color=out_lds$orig), type="norm")+
  scale_color_manual(name= "Population", values = pop_cols)+
  theme_classic()+
  theme(legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))


ldas<-grid.arrange(LD1_2, LD1_3, ncol=2, nrow=1)
ggsave(file="Figure##_LDAs.svg", plot=ldas, width=10, height=5)

legend<-ggpubr::get_legend(LD_legend)
legend_plot<-as_ggplot(legend)
ggsave(file="Figure##_LDAs_legend.svg", plot=legend_plot, width=7, height=0.5) 


ld_df<-rbind(z1_lds, out_lds)


#out2 <- out$posterior
#out3 <- data.frame(out2)
#out3$reg <- gsub('[[:digit:]]+',"",rownames(out3))
#out3$reg <- rownames(out3)
#md <- melt(out3)
#md2 <- dcast(reg~variable,data=md,sum)
#md2 <- data.frame(md2)
#rownames(md2) <- md2[,1]
#md2 <- md2[,-1]
#md3 <- data.frame(t(md2))
#md3 <- cbind(md3,data.frame(matrix(rep(0,4*10),ncol=10)))
#tmp <- data.frame(matrix(rep(0,4*14),nrow=4));rownames(tmp) <- c("CWFL","DRTO","NGM", "SWFL");colnames(tmp) <- colnames(md3)
#md3 <- rbind(md3,tmp)
#md3<-md3[1:8,1:4]
#md3 <- as.matrix(md3)
#par(mar = rep(0, 4), cex=0.9)
#circos.par(start.degree = 90, gap.degree = 4)
#cols.to.use <- pop_cols
#chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
            # grid.col=cols.to.use,
            # annotationTrack = "grid", 
             #transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
             #diffHeight  = -0.04,
            # direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
            # link.arr.length =  0.15)


#new Diversity MEtrics
dartR_div<-gl.report.heterozygosity(genlight_out)

#Calculate pairwise FST
pop.fst <- genet.dist(genind_out, method = "WC84") # Run time ~20 mins
pop.fst.df <- as.data.frame(as.matrix(round(pop.fst, 3)))
pop.fst.df #Look at matrix of Fst values
write.csv(pop.fst.df,"fst_df.csv")
#pop.fst.df<-read.csv(file="fileName.csv", header = T)

#Format Fst dataframe
pop.fst.tri <- pop.fst.df
pop.fst.tri[lower.tri(pop.fst.df, diag=TRUE)] <- NA
fst.mat = data.matrix(pop.fst.tri)
melted <- melt(fst.mat, na.rm =TRUE)

#Plot heatmap of Fst among populations
ggplot(data = melted, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "blue", high = "red", name="FST")  + ggtitle(expression(atop("Pairwise FST, WC (1984)", atop(italic("N = 75, L = 28,710"), ""))))+labs( x = "Sampling Site", y = "Sampling Site") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed()

##BootStrap + Pvalues
#All Data
fst_boot<-gl.fst.pop(genlight_out, nboots=10000, percent=95)
fst_boot_fsts<-data.matrix(round(fst_boot[["Fsts"]], 3))
fst_boot_ps<-data.matrix(round(fst_boot[["Pvalues"]], 10))
fst_boot_fsts
fst_boot_ps

#####DAPC/GSI######
#Split genind by randomly assigned training and validation groups
train.id<-row.names(pca1.scores[pca1.scores$lda_status=="T",])
val.id<-row.names(pca1.scores[pca1.scores$lda_status=="V",])
training_genind<-genind_out[train.id]
validation_genind<-genind_out[val.id]
  
##Initial DAPC
n_individuals <- nInd(training_genind)
n_pops <- length(levels(factor(popmap.filtered$STRATA)))

grp_all <- find.clusters(training_genind, max.n.clust=10, n.pca=200,
                         choose.n.clust = FALSE)

BIC<-as.data.frame(cbind(seq(1,10,1), grp_all$Kstat))

ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K)") + font("xlab", face = "bold") + font("ylab", face = "bold")
grp_all$Kstat

k <- 4
grp_all <- find.clusters(training_genind, max.n.clust=n_pops, n.pca=200, n.clust = k)
grp_all$size

n.da <- 6 #Set number of DAs
dapc <- dapc(training_genind, grp_all$grp, n.pca=60, n.da = n.da, var.contrib = TRUE)
pal <- get_palette("Accent", k = k)
dpca_result <- scatter(dapc, col=pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)
set.seed(4)
#contrib<-loadingplot(dapc$var.contr,axis=2, thres=.002,lab.jitter=1)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",c(1:k)),
          ncol=1, xlab="Individuals", col=pal, lab=training_genind@pop, show.lab = T)

dapc$IND <- row.names(dapc$ind.coord)

dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))

colnames(dapc_info) <- c("IND","DPC1", "DPC2", "DPC3", "K")
colnames(dapc_info) <- c("IND", c(paste0("DPC", 1:(n.da-3))), "K")
dapc_info$SITE <- str_sub(dapc_info$IND, 1,2) #Takes first two letters of the Indv name and adds to the columns "SITE"
dapc_info

table.value(table(pop(training_genind), grp_all$grp), col.lab=paste("inf", 1:4),
            row.lab=paste("ori", 1:4))

#Second DAPC
set.seed(5); dapc_a_score <- dapc(training_genind, n.pca=100, n.da=n.da) 
temp_score <- optim.a.score(dapc_a_score) #Determine best number of PCs to retain for DAPC analysis
n.pc <- temp_score$best #Save best number of PCs
dapc2 <-dapc(training_genind,training_genind@pop, n.pca = n.pc, n.da = n.da)
dapc2$IND <- row.names(dapc2$ind.coord)
dapc2_info <- as.data.frame(cbind(dapc2$IND, dapc2$ind.coord, grp_all$grp))

dpca_result <- scatter(dapc2, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)

jpeg(filename="out_dapc_scatter.jpeg", width=10, height=8, units="in", res = 400)
scatter(dapc2,
        scree.pca = TRUE,
        pch = pop_shapes,
        cell = 0,
        cstar = 0,
        cellipse = 1.5,
        solid = 0.8,
        cex = 3,
        clab = 0,
        ratio.da = 0.15,
        ratio.pca = 0.15,
        col = pop_cols)
dev.off()


contrib1 <- loadingplot(dapc2$var.contr, axis=1,
                       thres=.007, lab.jitter=1)
contrib2 <- loadingplot(dapc2$var.contr, axis=2,
                        thres=.007, lab.jitter=1)
contrib3 <- loadingplot(dapc2$var.contr, axis=3,
                        thres=.007, lab.jitter=1)

da_loci <- c(contrib1$var.names, contrib2$var.names, contrib3$var.names)
da_loci <- as.data.frame(da_loci)
da_loci$loci <- gsub("\\...*","",da_loci$da_loci)
da_loci <- da_loci %>% distinct(loci, .keep_all = TRUE)
da_loci<-da_loci[,2]

genlight_dapc <-gl.keep.loc(genlight_2, loc.list = da_loci, verbose=0)
genind_dapc <- gl2gi(genlight_dapc)

load_dpca2 <- as.data.frame(dapc2$var.contr)
write.table(load_dpca2, "Loadings_SNPs.txt", sep="\t", row.names=FALSE, quote=FALSE)
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, main = "Percent of genetic variance explained by eigenvectors", ylim = c(0, max(percent) + 10), names.arg = paste0("E", c(1:length(percent))), ylab = "Percent", xlab = "Eigenvector")
dapc_prior=as.data.frame(dapc2$ind.coord)
write.table(dapc_prior, "DAPC_prior_results.txt", quote=F, sep="\t", row.names=TRUE)

### Add information to the tab results of the DPCA
dapc_prior$IND <- row.names(dapc_prior)

### Add site info
dapc_prior$SITE <- str_sub(dapc_prior$IND, 1,2) #Again, assumes the site is the first two characters of the individual names
#dapc_prior$SITE <- gsub ('_', '', dapc_prior$SITE)

### Add 'region' info:
#dapc_prior <- merge(dapc_prior, region, by.x = "SITE", by.y = "Population")

### Make a ggplot graph representing the DAPC for the first and second axes for the regions
ggscatter(dapc_prior, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 3,
          xlab = "DA1", ylab = "DA2", palette = pop_cols)
summary(dapc2)
compoplot(dapc2)

##Make DF for training assignments
assignments<-as.data.frame(dapc2$posterior)
assignments$pop<-str_sub(dapc_prior$IND, 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
assignments$ind<-rownames(assignments)
assignments$post<-colnames(assignments[,1:4])[max.col(assignments[,1:4],ties.method="first")]
assignments$train_val<-"train"

##Predict
dapc_pred<-predict.dapc(dapc2, newdata=validation_genind)
dapc_pred$posterior

##Make DF for validation assignments
val_assignments<-as.data.frame(dapc_pred$posterior)
val_assignments$pop<-str_sub(row.names(dapc_pred$ind.scores), 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
val_assignments$ind<-rownames(val_assignments)
val_assignments$post<-colnames(val_assignments[,1:4])[max.col(val_assignments[,1:4],ties.method="first")]
val_assignments$train_val<-"val"

whole_ass<-rbind(assignments, val_assignments)
ass_melt<-melt(whole_ass)

ass_melt<-ass_melt%>%
  arrange(pop) %>%
  arrange(ind)%>%
  arrange(train_val)
ass_ind_table <- table(ass_melt$ind)
ind_levels <- names(ass_ind_table)[order(ass_ind_table)]
ass_melt$ind <- factor(ass_melt$ind, levels = ind_levels)

compoplot_dapc2<-ggplot(ass_melt, 
                        aes(fill=variable, y=value, x=ind),
                        color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols)+
  theme(legend.position="none") +
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank())

ggsave(file="Figure##_DAPC_out.svg", plot=compoplot_dapc2, width=10, height=5) 

compoplot_dapc2_legend<-ggplot(ass_melt, 
                               aes(fill=variable, y=value, x=ind),
                               color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols, name= "Population")+
  theme(legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20))+
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank())

dev.off()
legend_compoplot<-ggpubr::get_legend(compoplot_dapc2_legend)
legend_comp_plot<-as_ggplot(legend_compoplot)
ggsave(file="Figure2_DAPC_legend.svg", plot=legend_comp_plot) 

###########################################################
##second outlier method:OutFLANK
#########################################################
outflank<-gl.outflank(genind)
sum(outflank[["outflank"]][["results"]][["qvalues"]]<0.8, na.rm=TRUE)

plot(outflank$outflank$results$He, outflank$outflank$results$FST, pch=20, col="grey")
points(outflank$outflank$results$He[outflank$outflank$results$qvalues<0.8], outflank$outflank$results$FST[outflank$outflank$results$qvalues<0.8], pch=21, col="blue")
top_candidates <- outflank$outflank$results$qvalues<0.8 & outflank$outflank$results$He>0.1
topcan <- outflank$outflank$results[top_candidates,]
topcan[order(topcan$LocusName),]      

flank_loci <- gsub("\\...*","",topcan$LocusName)
flank_loci <- as.data.frame(flank_loci)
flank_loci <- flank_loci %>% distinct(flank_loci, .keep_all = TRUE)
flank_loci<- as.matrix(flank_loci)
genlight_flank <-gl.keep.loc(genlight_2, loc.list = flank_loci, verbose=0)
genind_flank <- gl2gi(genlight_flank)


##PCA
##PCAs and Plots##
pca1 <- glPca(genlight_flank,center = T, scale = T, nf = 5)
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

# Extract PC scores to color by location/yr:
# Adapted from https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$Population <- pop(genlight_flank)
pca1.scores$reg <- pca1.scores$pop

#Make PCA plots
set.seed(9)
num_pops <- length(levels(factor(popmap.filtered$STRATA)))
pop_cols <- c("#000000","#E69F00" , "#56B4E9"  ,"#009E73")
pop_shapes<-c(15,16,17,18)

# plot PC 1 and 2
PCA1_2<-ggscatter(pca1.scores, x = "PC2", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)")) +
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  #scale_y_continuous(limits=c(-110,80), breaks=c(-75, -50, -25, 0, 25, 50),
  #labels=c("-75", "-50", "-25", "0", "25", "50"))+
  theme(panel.grid = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_2

# #plot  PC 1 and 3
PCA1_3<-ggscatter(pca1.scores, x = "PC3", y = "PC1", shape = "Population", color = "Population",
                  palette = pop_cols, size=4, ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"),
                  ylab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"))+
  scale_shape_manual(values=pop_shapes) +
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
PCA1_3

#LDA
#Create random list of training and validation individuals, using half of each population
ctemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "CWFL"))), 0)))
dtemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "DRTO"))), 0)))
ntemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "NGM"))), 0)))
stemp<-sample(rep(c("T", "V"), (0.5*length(which(pca1.scores$Population == "SWFL"))))) #no round because odd number
extra<-sample (c("T","V"), 1) #add because SWFL odd
tv_list<-c(ctemp, dtemp, ntemp,stemp, extra)

pca1.scores$lda_status<-c(tv_list)
pca1.scores$lda_status<-as.factor(pca1.scores$lda_status)
training<-which(pca1.scores$lda_status == "T")
validation<-which(pca1.scores$lda_status == "V")
pca1.scores.lda<-pca1.scores[,-7]

z <- lda(Population~.,pca1.scores.lda, subset=training)
z1<-predict(z,pca1.scores.lda[training,])
z1$class==pca1.scores.lda[training,]$Population
z1_lds<-as.data.frame(z1$x)
z1_lds$lda_status<-c("Training")
z1_lds$orig<-gsub('[[:digit:]]+',"",rownames(z1_lds))
z1_lds$post<-z1$class

out <- predict(z,pca1.scores.lda[validation,])
out$class==pca1.scores.lda[validation,]$Population
out_lds<-as.data.frame(out$x)
out_lds$lda_status<-c("Validation")
out_lds$orig<-gsub('[[:digit:]]+',"",rownames(out_lds))
out_lds$post<-out$class

LD1_2<-ggplot() +
  geom_point(data=z1_lds, aes(x=LD2, y=LD1, color=orig), shape = 15, size = 4, alpha=0.5)+
  stat_ellipse(aes(x=z1_lds$LD2, y=z1_lds$LD1, color=z1_lds$orig), type="norm", linetype=2)+
  geom_point(data=out_lds, aes(x=LD2, y=LD1, color=orig), shape = 16, size = 4)+
  stat_ellipse(aes(x=out_lds$LD2, y=out_lds$LD1, color=out_lds$orig), type="norm")+
  scale_color_manual(values = pop_cols)+
  theme_classic()+
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
LD1_2

LD1_3<-ggplot() +
  geom_point(data=z1_lds, aes(x=LD3, y=LD1, color=orig), shape = 15, size = 4, alpha=0.5)+
  stat_ellipse(aes(x=z1_lds$LD3, y=z1_lds$LD1, color=z1_lds$orig), type="norm", linetype=2)+
  geom_point(data=out_lds, aes(x=LD3, y=LD1, color=orig), shape = 16, size = 4)+
  stat_ellipse(aes(x=out_lds$LD3, y=out_lds$LD1, color=out_lds$orig), type="norm")+
  scale_color_manual(values = pop_cols)+
  theme_classic()+
  theme(legend.position="none") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=25, face="bold"), 
        axis.text = element_text(size=25))
LD1_3

##DAPC
#Split genind by randomly assigned training and validation groups
train.id<-row.names(pca1.scores[pca1.scores$lda_status=="T",])
val.id<-row.names(pca1.scores[pca1.scores$lda_status=="V",])
training_genind<-genind_flank[train.id]
validation_genind<-genind_flank[val.id]

##Initial DAPC
n_individuals <- nInd(training_genind)
n_pops <- length(levels(factor(popmap.filtered$STRATA)))

grp_all <- find.clusters(training_genind, max.n.clust=10, n.pca=200,
                         choose.n.clust = FALSE)

BIC<-as.data.frame(cbind(seq(1,10,1), grp_all$Kstat))

ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K)") + font("xlab", face = "bold") + font("ylab", face = "bold")
grp_all$Kstat

k <- 4
grp_all <- find.clusters(training_genind, max.n.clust=n_pops, n.pca=200, n.clust = k)
grp_all$size

n.da <- 6 #Set number of DAs
dapc <- dapc(training_genind, grp_all$grp, n.pca=60, n.da = n.da, var.contrib = TRUE)
pal <- get_palette("Accent", k = k)
dpca_result <- scatter(dapc, col=pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)
set.seed(4)
#contrib<-loadingplot(dapc$var.contr,axis=2, thres=.002,lab.jitter=1)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",c(1:k)),
          ncol=1, xlab="Individuals", col=pal, lab=training_genind@pop, show.lab = T)

dapc$IND <- row.names(dapc$ind.coord)

dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))

colnames(dapc_info) <- c("IND","DPC1", "DPC2", "DPC3", "K")
colnames(dapc_info) <- c("IND", c(paste0("DPC", 1:(n.da-3))), "K")
dapc_info$SITE <- str_sub(dapc_info$IND, 1,2) #Takes first two letters of the Indv name and adds to the columns "SITE"
dapc_info

table.value(table(pop(training_genind), grp_all$grp), col.lab=paste("inf", 1:4),
            row.lab=paste("ori", 1:4))

#Second DAPC
set.seed(5); dapc_a_score <- dapc(training_genind, n.pca=100, n.da=n.da) 
temp_score <- optim.a.score(dapc_a_score) #Determine best number of PCs to retain for DAPC analysis
n.pc <- temp_score$best #Save best number of PCs
dapc2 <-dapc(training_genind,training_genind@pop, n.pca = n.pc, n.da = n.da)
dapc2$IND <- row.names(dapc2$ind.coord)
dapc2_info <- as.data.frame(cbind(dapc2$IND, dapc2$ind.coord, grp_all$grp))

dpca_result <- scatter(dapc2, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)

contrib1 <- loadingplot(dapc2$var.contr, axis=1,
                        thres=.007, lab.jitter=1)
contrib2 <- loadingplot(dapc2$var.contr, axis=2,
                        thres=.007, lab.jitter=1)
contrib3 <- loadingplot(dapc2$var.contr, axis=3,
                        thres=.007, lab.jitter=1)

da_loci <- c(contrib1$var.names, contrib2$var.names, contrib3$var.names)
da_loci <- as.data.frame(da_loci)
da_loci$loci <- gsub("\\...*","",da_loci$da_loci)
da_loci <- da_loci %>% distinct(loci, .keep_all = TRUE)
da_loci<-da_loci[,2]

genlight_dapc <-gl.keep.loc(genlight_2, loc.list = da_loci, verbose=0)
genind_dapc <- gl2gi(genlight_dapc)

load_dpca2 <- as.data.frame(dapc2$var.contr)
write.table(load_dpca2, "Loadings_SNPs.txt", sep="\t", row.names=FALSE, quote=FALSE)
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, main = "Percent of genetic variance explained by eigenvectors", ylim = c(0, max(percent) + 10), names.arg = paste0("E", c(1:length(percent))), ylab = "Percent", xlab = "Eigenvector")
dapc_prior=as.data.frame(dapc2$ind.coord)
write.table(dapc_prior, "DAPC_prior_results.txt", quote=F, sep="\t", row.names=TRUE)

### Add information to the tab results of the DPCA
dapc_prior$IND <- row.names(dapc_prior)

### Add site info
dapc_prior$SITE <- str_sub(dapc_prior$IND, 1,2) #Again, assumes the site is the first two characters of the individual names
#dapc_prior$SITE <- gsub ('_', '', dapc_prior$SITE)

### Add 'region' info:
#dapc_prior <- merge(dapc_prior, region, by.x = "SITE", by.y = "Population")

### Make a ggplot graph representing the DAPC for the first and second axes for the regions
ggscatter(dapc_prior, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 3,
          xlab = "DA1", ylab = "DA2", palette = pop_cols)
summary(dapc2)
compoplot(dapc2)

##Make DF for training assignments
assignments<-as.data.frame(dapc2$posterior)
assignments$pop<-str_sub(dapc_prior$IND, 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
assignments$ind<-rownames(assignments)
assignments$post<-colnames(assignments[,1:4])[max.col(assignments[,1:4],ties.method="first")]
assignments$train_val<-"train"

##Predict
dapc_pred<-predict.dapc(dapc2, newdata=validation_genind)
dapc_pred$posterior

##Make DF for validation assignments
val_assignments<-as.data.frame(dapc_pred$posterior)
val_assignments$pop<-str_sub(row.names(dapc_pred$ind.scores), 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
val_assignments$ind<-rownames(val_assignments)
val_assignments$post<-colnames(val_assignments[,1:4])[max.col(val_assignments[,1:4],ties.method="first")]
val_assignments$train_val<-"val"

whole_ass<-rbind(assignments, val_assignments)
ass_melt<-melt(whole_ass)

ass_melt<-ass_melt%>%
  arrange(pop) %>%
  arrange(ind)%>%
  arrange(train_val)
ass_ind_table <- table(ass_melt$ind)
ind_levels <- names(ass_ind_table)[order(ass_ind_table)]
ass_melt$ind <- factor(ass_melt$ind, levels = ind_levels)

compoplot_dapc2<-ggplot(ass_melt, 
                        aes(fill=variable, y=value, x=ind),
                        color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols)+
  theme(legend.position="none") +
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank())
compoplot_dapc2


###########################################################
###Rubias##
###########################################################
#Prior to This, gl2genalex, get rid of indiv and headers in excel
loc_names<-as.character(unique(genind_out$loc.fac))
locnames.1<-paste0(loc_names,".1")
locnames_rdy<-c(rbind(loc_names,locnames.1))

rubias_data<-read.csv("gl.csv", header = F)
colnames(rubias_data)<-locnames_rdy #add locus names

#get sample names into proper column
rubias_data$indiv<-row.names(genind_out$tab)

#Need to create list of reference and mixture individuals
ctemp<-sample(rep(c("reference", "mixture"), round((0.5*length(which(pca1.scores$Population == "CWFL"))), 0)))
dtemp<-sample(rep(c("reference", "mixture"), round((0.5*length(which(pca1.scores$Population == "DRTO"))), 0)))
ntemp<-sample(rep(c("reference", "mixture"), round((0.5*length(which(pca1.scores$Population == "NGM"))), 0)))
stemp<-sample(rep(c("reference", "mixture"), (0.5*length(which(pca1.scores$Population == "SWFL"))))) #no round because odd number
extra<-sample (c("reference","mixture"), 1) #add because SWFL odd
samptype_list<-c(ctemp, dtemp, ntemp,stemp, extra)

rubias_data$sample_type<-samptype_list

#RepUnit: closely related populations that might aggregate. Will try with GoM. Only for mixture samples
refrows<-as.numeric(which(rubias_data$sample_type=="reference")) #for downstream row selections
mixrows<-as.numeric(which(rubias_data$sample_type=="mixture")) #for downstream row selections

rubias_data$repunit<- ifelse(rubias_data$sample_type == "reference", "GoM", NA)

#Collection is the population
rubias_data$collection<-ifelse(rubias_data$sample_type == "reference", gsub("[^a-zA-Z]", "",rubias_data$indiv), "GoM")

#Order columns
rubias_data_rdy<- rubias_data %>%
  relocate(indiv) %>%
  relocate(collection) %>%
  relocate(repunit) %>%
  relocate(sample_type)

rubias_data_rdy<-as_tibble(rubias_data_rdy)

#Make mixture and reference tibbles
rubias_mixture<-rubias_data_rdy %>% dplyr::filter(sample_type=="mixture")
rubias_ref<-rubias_data_rdy %>% dplyr::filter(sample_type=="reference")

##Check self assignment of references
sa_gom <- self_assign(reference = rubias_ref, gen_start_col = 5)

##Conduct mixture analysis
mix_est <- infer_mixture(reference = rubias_ref, 
                         mixture = rubias_mixture, 
                         gen_start_col = 5,
                         reps = 100000,
                         burn_in = 10000,
                         pb_iter = 1000,
                         method= "PB")
rubias_ass<-mix_est$indiv_posteriors

ggplot(rubias_ass, 
       aes(fill=collection, y=PofZ, x=indiv),
       color = "variable")+ 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=pop_cols)+
  theme(legend.position="none") +
  labs(y="Membership Probability", x="Individuals")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=20, face="bold"), 
        axis.text = element_text(size=20),
        axis.text.x = element_blank())

nsweeps <- max(mix_est$mix_prop_traces$sweep)
trace_subset <- mix_est$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, collection) %>%
  summarise(repprop = sum(pi))

ggplot(trace_subset, aes(x = repprop, colour = collection)) +
  geom_density()


###POPHELPER##

#setwd("./CCATL_harvester/CCATL_structureselector/")
#setwd("./COMBINED_harvester/Combined_cut_Harvester_Upload/harvester_input")

setwd("C:/Users/Ian/Desktop/Dissertation/Chapter 4 Genetics/Genomics Metadata and Data/VCFTools_Out/")
setwd("./CCGOM_harvester/CCGOM_cut_Harvester_Upload/harvester_input/")

setwd("./CCGOM_fastStructure")
#setwd("./COMBINED_fastStructure")
setwd("./CCATL_fastStructure")

sfile<-list.files(path="./") #, pattern="\\.meanQ$")

#Different inputs for STRUCTURE vs fastStructure

slist<-readQ(files=sfile, filetype="structure", indlabfromfile = TRUE)

slist<-readQ(files=sfile)

#For non-structure files, load individual names list and add to slist

inds<-read.delim("C:/Users/Ian/Desktop/Dissertation/Chapter 4 Genetics/Genomics Metadata and Data/VCFTools_Out/data/CCGOM_filt_indnames_sorted.txt", header=F)
#inds<-read.delim("Combined_sortednames.txt", header=F)
inds<-read.delim("CCATL_sortednames.txt", header=F)

if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$V1)

tab_slist<-tabulateQ(slist)
sum_slist<-summariseQ(tab_slist)

#Set WD for plot outputs
#setwd("C:/Users/Ian/Desktop/Dissertation/Chapter 4 Genetics/Genomics Metadata and Data/VCFTools_Out/COMBINED_harvester")

#Evanno Method only works for STRUCTURE runs
em <- evannoMethodStructure(data=sum_slist,exportplot=T,writetable=T,na.rm=T,exportpath=getwd())
p <- evannoMethodStructure(data=sum_slist,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

slist1 <- alignK(slist)

p1<- plotQ(slist1[c(1,11,21,31,41,51,61,71,81,91)],
           imgoutput="join",
           returnplot=T,
           exportplot=F,
           basesize=11,
           showindlab = FALSE,
           useindlab = FALSE,
           ordergrp = TRUE,
           #sortind = "all",
           sharedindlab = TRUE)

jpeg("CCGOM_Structure_k1-10.jpg", width = 15, height = 10, units = "in", res = 500) 
grid.arrange(p1$plot[[1]])
dev.off()

###Relatedness###
#For this analysis, using KING estimation as well as kinship from Combined dataset

rel_cols_2<-c("#F8766DFF", "#7CAE00FF")
rel_cols_4<-c("#F8766DFF", "#7CAE00FF", "#00BFC4FF", "#C77CFFFF")

data<-read.csv("./COMBINED_Relatedness/Combined_ngsrelate.csv")
data<-data %>%
  dplyr::select(c("ida", "idb", "R0", "R1", "KING"))
data$label <- paste(data$ida, data$idb)

plink<-read.csv("./COMBINED_Relatedness/Combined_plink.csv", header=T)
plink<-plink %>%
  dplyr::select(c("IID1", "IID2", "Z0", "Z1", "Z2", "PI_HAT"))
plink$label<-paste(plink$IID1, plink$IID2)

alldata<-merge(data, plink)
alldata$kinship<-((alldata$Z1/4)+(alldata$Z2/2))
#Cutoffs = Unrelated < 0.044 < 3rd  < 0.088 < 2nd < 0.177 < FS/PO < 0.354 < Twins
alldata$Relationship<-cut(alldata$kinship,
                          breaks=c(-Inf, (1/(2^(9/2))), (1/(2^(7/2))), (1/(2^(5/2))), (1/(2^(3/2))), Inf),
                          labels=c("Unrelated", "3rd Degree", "2nd Degree", "FS/PO", "Twins"))
alldata$kplus0_01 <- alldata$kinship + 0.01
alldata$kminus0_01 <- alldata$kinship - 0.01

# Create a function to extract population abbreviations
extract_population <- function(individual) {
  # Extract population abbreviation by splitting at the numeric part
  Population <- gsub("\\d", "", individual)
  return(Population)
}

# Create a new column 'population'
alldata <- alldata %>%
  mutate(Population = ifelse(extract_population(ida) == extract_population(idb), 
                             extract_population(ida), 
                             paste(extract_population(ida), "-", extract_population(idb), sep = "")))


#Subset for higher KING values (may not be necessary if you have already filtered relatives)
#alldata_pos<- alldata %>% subset(KING>-0.4)

#Subset for internal population comparisons
alldata_internal_pop <- alldata %>%
  filter(!grepl("-", Population))

alldata_internal_pop$Population <- factor(alldata_internal_pop$Population, levels = c("NRU", "CEFL", "SEFL", "DRTO", "SWFL", "CWFL", "NGM"))

jpeg("Combined_KINGvsR1_Relatedness.jpg", width = 6, height = 6, units = "in", res = 500)
alldata_internal_pop %>% ggplot(aes(R1, KING, label=label)) +
  geom_point(aes(color=Population), size=3, alpha=0.6) +
  #scale_color_manual(values = rel_cols_4)+
  theme_classic()+
  theme(axis.title = element_text(size=14, face="bold"), 
        axis.text = element_text(size=14, face="bold"))
dev.off()

#        ,
#        legend.text = element_blank(),
#        legend.title = element_blank(),
#        legend.position = "none")

jpeg("Combined_KING__Violin_Relatedness.jpg", width = 6, height = 6, units = "in", res = 500)
alldata_internal_pop %>% 
  ggplot(aes(x=Population, y=KING)) + 
  geom_violin(trim=FALSE)+
  theme_classic()+
  geom_boxplot(width=0.1, color="black")
dev.off()

#Summary Statistics

population_stats <- alldata_internal_pop %>%
  group_by(Population) %>%
  summarize(
    MeanKING = mean(KING),
    MedianKING = median(KING),
    StdDevKING = sd(KING),
    MeanR1 = mean(R1),
    MedianR1 = median(R1),
    StdDevR1 = sd(R1)
  )

write.csv(population_stats, "Pop_Relatedness.csv")

#Subset for between population comparisons
alldata_between_pop <- alldata %>%
  filter(grepl("-", population))


alldata_between_pop %>% ggplot(aes(R1, KING, label=label)) +
  geom_point(aes(color=population), size=3, alpha=0.6) +
  #scale_color_manual(values = rel_cols_4)+
  theme_classic()+
  theme(axis.title = element_text(size=14, face="bold"), 
        axis.text = element_text(size=14, face="bold"))
        
#        ,
#        legend.text = element_blank(),
#        legend.title = element_blank(),
#        legend.position = "none")
