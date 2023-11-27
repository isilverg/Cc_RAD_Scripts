##GoM Cc PopGen Script##

#This script makes a map of sample sites (provided appropriate shapefiles),
#imports and filters raw data (based on retained individuals and loci from upstream filtering),
#and conducts various pop-gen analyses (PCA, DAPC, Fst, Diversity).
#Most results are exported as csv files or figures (jpeg) to the working directory.

#---------------Load Packages-------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
devtools::install_github("thierrygosselin/radiator")

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
library(stackoverflow)
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
library(rgdal)
library(circlize)
library(MASS)
#library(OutFLANK)
library(giscoR)
library(radiator)
library(rubias)

#---------------Map Study Sites-------------------------------------
#Load data
Excel <- read.csv("regions.csv") #excel containing the GPS points for each population
Excel$Population <- as.factor(Excel$Population)
Florida<-st_read("FL.shp") #shapefile for the state of Florida
Bathy<-st_read("bathy200m.shp") #shapefile for 200m bathymetry contour

#Create the main map
map <- ggplot(data=Bathy)+
  geom_sf() +
  coord_sf(xlim = c(-87, -80), ylim = c(24, 30.5))+
  geom_point(data=Excel, aes(lon, lat, color=Population, shape=Population), size=6.5, alpha=1)+
  scale_color_manual(values=pop_cols) +
  scale_shape_manual(values=pop_shapes)+
  geom_label(data=Excel, aes(lon, lat), label=Excel$Population, size=6, fontface="bold", nudge_x = -0.5, nudge_y = -0.15) +
  labs(x="Longitude", y="Latitude")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.title = element_text(size=22, face="bold"), 
        axis.text = element_text(size=22, face="bold"),
        legend.position = "none") +
  annotation_scale(location="bl", 
                   width_hint=0.4,
                   height = unit(0.15, "in"), 
                   text_cex = 1.5,
                   text_face = "bold")+
  annotation_north_arrow(location="bl", which_north="true",
                         pad_x = unit(0.0, "in"), pad_y = unit(0.23, "in"),
                         style=north_arrow_fancy_orienteering)
#View the map
map 

#Background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = "+proj=ortho +lat_0=25 +lon_0=-80 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

#Country polygons, cut to size
world <- gisco_countries %>%
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = "+proj=ortho +lat_0=25 +lon_0=-80 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

#Create bbox for NWADPS
bbox <- st_sfc(st_point(c(-98.5, 35)), st_point(c(-60, 16.5)), crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(crs = "+proj=ortho +lat_0=25 +lon_0=-80 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

#Create bbox for FL GOM
bbox_fl <- st_sfc(st_point(c(-89, 31.5)), st_point(c(-80, 23.5)), crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(crs = "+proj=ortho +lat_0=25 +lon_0=-80 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

#Create inset map
inset.plot <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "lightblue", color = NA) + # background first
  geom_sf(lwd = 0.1) + # now land over the oceans
  geom_sf(data = bbox, fill = "transparent", color = "black", lty="11", lwd=0.75) +
  geom_sf(data = bbox_fl, fill = "transparent", color = "red", lty=1, lwd=1) +
  coord_sf(expand = FALSE) +
  theme_void()

#View inset plot
inset.plot

#Align map and inset
full_map<-map+inset_element(inset.plot, left = 0.15, bottom = 0.15, right = 0.55, top = 0.55, align_to = "plot")

#Export Figure 1: Map#
jpeg("Figure1_Map.jpg", width = 10, height = 10, units = "in", res = 500) 
full_map
dev.off()

#-----------------Import and Filter SNP Data----------------------
setwd("FilePath/With/RawVCF/And/Whitelist/And/PopMap/And/RegionFile")

#Read in VCF file
in_vcf <- read.vcfR("your_filename.vcf", convertNA = T)

#Convert vcf to genlight. Save twice: genlight.orig will be preserved if samples are subset and the genlight object will be replaced; if no subsetting is needed, then the file named genlight will be used throughout the rest of the script. 
genlight.orig = genlight <- vcfR2genlight(in_vcf)

#Read in the population map. This should be a two-column tab-delimited text file where the first column contains the name of each individual sample and the second column contains the population. This can be the same as the population map used in Stacks.
popmap.orig <- read_tsv("your_filename.txt", col_names = F) %>% 
  dplyr::rename(Indiv = X1, STRATA = X2)

#Filter the population map for whitelisted individuals i.e. those we will keep for the population genetics analysis.
popmap.filtered <- tibble(Indiv = genlight@ind.names) %>% left_join(popmap.orig, by = "Indiv")
genlight@pop <- as.factor(popmap.filtered$STRATA)
summary(genlight@pop)

#Read in a "regions" file. This file should have at a minimum a column for Population, latitude, and longitude. We assume this is a csv (comma-separate) file.
region <- read.csv("your_filename.csv", header = T)

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

#Check compliance so that individuals and loci can be filtered
genlight_2 <- gl.compliance.check(genlight)

#Subset the full set of samples for analysis
#Only run this section if you want to run the analysis on a subset of individuals included in the filtered VCF file. If you want to run on all the samples in the vcf file that was loaded, then skip this section and move on to "Finish creating genlight and genind objects" a few lines below.
#Read in white list of individuals/locs. This should be a text file with a single column containing the names of the individuals that were retained after filtering.
whitelist.indvs <- read.delim("your_filename.txt", header = F, strip.white = T) %>% pull()
genlight_3 <- gl.keep.ind(genlight_2, ind.list = whitelist.indvs)
genlight_3@ind.names

whitelist.locs <- read.delim("your_filename.txt", header = F, strip.white = T) %>%
  dplyr::rename(Chrom = V1, Pos = V2) %>%
  arrange(Chrom, Pos)
whitelist.locs_2 <- str_c(whitelist.locs$Chrom, "_", whitelist.locs$Pos)
genlight_4 <-gl.keep.loc(genlight_3, loc.list = whitelist.locs_2, verbose=0)
rm(genlight, genlight_2, genlight_3)

#Finish creating genlight and genind objects
genind <- gl2gi(genlight_4) # need genind format for some functions
summary(genlight_4@pop)
summary(genind@pop)

#-------------Analysis--------------------------------------------------------------------

#Some Basic Stats##
dartR_het<-gl.report.heterozygosity(genlight_4)
dartR_div<-gl.report.diversity(genlight_4)

##PCAs and Plots##
pca1 <- glPca(genlight_4,center = T, scale = T, nf = 5)
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

# Extract PC scores to color by location/yr:
# Adapted from https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$Population <- pop(genlight_4)
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

#Make Figure of PCAs
dev.off()
legend<-ggpubr::get_legend(PCA_Legend)
legend_plot<-as_ggplot(legend)
ggsave(file="Figure1_PCAs_legend.svg", plot=legend_plot, width=7, height=0.5) 

#Initial DAPC
n_individuals <- nrow(popmap.filtered)
n_pops <- length(levels(factor(popmap.filtered$STRATA)))

#Search for optimal clusters in data
grp_all <- find.clusters(genind, max.n.clust=10, n.pca=200,
                         choose.n.clust = FALSE)

BIC<-as.data.frame(cbind(seq(1,10,1), grp_all$Kstat))

ggline(BIC, x = "V1", y = "V2", plot_type = "b",
       col = "navy",
       xlab = "Number of clusters (K)",
       ylab = "BIC Value",
       title = "Selection of optimum number of clusters (K)") + font("xlab", face = "bold") + font("ylab", face = "bold")
grp_all$Kstat

#Set k to optimal # of clusters as indicated by minimum in prior graph
#If none, set k to number of populations
k <- 4
grp_all <- find.clusters(genind, max.n.clust=n_pops, n.pca=200, n.clust = k)
grp_all$size

n.da <- 6 #Set number of DAs
dapc <- dapc(genind, grp_all$grp, n.pca=60, n.da = n.da, var.contrib = TRUE)
pal <- get_palette("Accent", k = k)

#Visualize first DAPC
dapc_result <- scatter(dapc, col=pal, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)
set.seed(4)
compoplot(dapc, posi="topleft", txt.leg=paste("Cluster",c(1:k)),
          ncol=1, xlab="Individuals", col=pal, lab=genind_out@pop, show.lab = T)

#Add individual names
dapc$IND <- row.names(dapc$ind.coord)

#Create new dataframe for custom plotting
dapc_info <- as.data.frame(cbind(dapc$IND, dapc$ind.coord, grp_all$grp))
colnames(dapc_info) <- c("IND","DPC1", "DPC2", "DPC3", "K")
colnames(dapc_info) <- c("IND", c(paste0("DPC", 1:(n.da-3))), "K")
dapc_info$SITE <- str_sub(dapc_info$IND, 1,2) #Takes first two letters of the Indv name and adds to the columns "SITE"
dapc_info

#See where individuals fall out in DA by populations
table.value(table(pop(genind), grp_all$grp), col.lab=paste("inf", 1:4),
            row.lab=paste("ori", 1:4))

#Second, calibrated DAPC
set.seed(5); dapc_a_score <- dapc(genind, n.pca=100, n.da=n.da) 
temp_score <- optim.a.score(dapc_a_score) #Determine best number of PCs to retain for DAPC analysis
n.pc <- temp_score$best #Save best number of PCs
dapc2 <-dapc(genind,genind@pop, n.pca = n.pc, n.da = n.da)
dapc2$IND <- row.names(dapc2$ind.coord)
dapc2_info <- as.data.frame(cbind(dapc2$IND, dapc2$ind.coord, grp_all$grp))

#Visualize second DAPC
dapc_result <- scatter(dapc2, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)

percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, main = "Percent of genetic variance explained by eigenvectors", ylim = c(0, max(percent) + 10), names.arg = paste0("E", c(1:length(percent))), ylab = "Percent", xlab = "Eigenvector")
dapc_final=as.data.frame(dapc2$ind.coord)
write.table(dapc_final, "DAPC_final_results.txt", quote=F, sep="\t", row.names=TRUE)

#Add information to the tab results of the DAPC
dapc_final$IND <- row.names(dapc_final)

#Add site info
dapc_final$SITE <- str_sub(dapc_final$IND, 1,2) #Assumes the site is the first two characters of the individual names


#Add 'region' info if desired
#dapc_final <- merge(dapc_final, region, by.x = "SITE", by.y = "Population")

#Make a ggplot graph representing the DAPC for the first and second axes for the regions
ggscatter(dapc_final, x = "LD1", y = "LD2", shape = 21, fill = "SITE", size = 3,
          xlab = "DA1", ylab = "DA2", palette = pop_cols)

#Sumarize the dapc
summary(dapc2)

#View assignment barplot
compoplot(dapc2)

#Make custom graph of second DAPC
assignments<-as.data.frame(dapc2$posterior)
assignments$pop<-str_sub(dapc_prior$IND, 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
assignments$ind<-rownames(assignments)
assignments$post<-colnames(assignments[,1:4])[max.col(assignments[,1:4],ties.method="first")]
assign_melt<-melt(assignments)

assign_melt<-assign_melt%>%
  arrange(pop) %>%
  arrange(ind)
assign_ind_table <- table(assign_melt$ind)
ind_levels <- names(assign_ind_table)[order(assign_ind_table)]
assign_melt$ind <- factor(assign_melt$ind, levels = ind_levels)

compoplot_dapc2<-ggplot(assign_melt, 
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

ggsave(file="Figure2_DAPC.svg", plot=compoplot_dapc2, width=10, height=5) 

compoplot_dapc2_legend<-ggplot(assign_melt, 
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

#################################################################################################
#####-----GSI------##############################################################################
#################################################################################################

#Convert genind to genpop for allele frequency difference calculations
genpop<-genind2genpop(genind)

#Calculate allele frequencies
freq<-makefreq(genpop, missing=0)

#Calculate allele frequency differences for all pairs of populations
freqdiffs<-list()
for (i in 1:(nrow(freq) - 1)) {
  for (j in (i + 1):nrow(freq)) {
    freqdiffs[[paste(row.names(freq)[i], row.names(freq)[j], sep = '-')]] <- abs(freq[i, ] - freq[j, ])
  }
}
as.data.frame(do.call(rbind, freqdiffs))

##For each pair of populations, format allele frequency differences by locus, 
##and select the top 40 rows as outliers

CWFL_DRTO <- as.data.frame(freqdiffs[["CWFL-DRTO"]])
colnames(CWFL_DRTO)<-"CWFLDRTO"
CWFL_DRTO$locus <- gsub("\\..*", "", rownames(CWFL_DRTO))
CWFL_DRTO<-CWFL_DRTO %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_DRTO)<-CWFL_DRTO$locus
CWFLDRTO_top<-CWFL_DRTO %>% top_n(40, CWFLDRTO)

CWFL_NGM <- as.data.frame(freqdiffs[["CWFL-NGM"]])
colnames(CWFL_NGM)<-"CWFLNGM"
CWFL_NGM$locus <- gsub("\\..*", "", rownames(CWFL_NGM))
CWFL_NGM<-CWFL_NGM %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_NGM)<-CWFL_NGM$locus
CWFLNGM_top<-CWFL_NGM %>% top_n(40, CWFLNGM)

CWFL_SWFL <- as.data.frame(freqdiffs[["CWFL-SWFL"]])
colnames(CWFL_SWFL)<-"CWFLSWFL"
CWFL_SWFL$locus <- gsub("\\..*", "", rownames(CWFL_SWFL))
CWFL_SWFL<-CWFL_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(CWFL_SWFL)<-CWFL_SWFL$locus
CWFLSWFL_top<-CWFL_SWFL %>% top_n(40, CWFLSWFL)

DRTO_NGM <- as.data.frame(freqdiffs[["DRTO-NGM"]])
colnames(DRTO_NGM)<-"DRTONGM"
DRTO_NGM$locus <- gsub("\\..*", "", rownames(DRTO_NGM))
DRTO_NGM<-DRTO_NGM %>% distinct(locus, .keep_all = TRUE)
rownames(DRTO_NGM)<-DRTO_NGM$locus
DRTONGM_top<-DRTO_NGM %>% top_n(40, DRTONGM)

DRTO_SWFL <- as.data.frame(freqdiffs[["DRTO-SWFL"]])
colnames(DRTO_SWFL)<-"DRTOSWFL"
DRTO_SWFL$locus <- gsub("\\..*", "", rownames(DRTO_SWFL))
DRTO_SWFL<-DRTO_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(DRTO_SWFL)<-DRTO_SWFL$locus
DRTOSWFL_top<-DRTO_SWFL %>% top_n(40, DRTOSWFL)

NGM_SWFL <- as.data.frame(freqdiffs[["NGM-SWFL"]])
colnames(NGM_SWFL)<-"NGMSWFL"
NGM_SWFL$locus <- gsub("\\..*", "", rownames(NGM_SWFL))
NGM_SWFL<-NGM_SWFL %>% distinct(locus, .keep_all = TRUE)
rownames(NGM_SWFL)<-NGM_SWFL$locus
NGMSWFL_top<-NGM_SWFL %>% top_n(40, NGMSWFL)

#Combine outlier loci into one dataframe
diff_loci<- CWFLDRTO_top %>% full_join(CWFLNGM_top, by="locus") %>% 
  full_join(CWFLSWFL_top, by="locus") %>%
  full_join(DRTONGM_top, by="locus") %>%
  full_join(DRTOSWFL_top, by="locus") %>%
  full_join(NGMSWFL_top, by="locus")

#Isolate/format outlier locus names for data subsetting
outlier_loci <- diff_loci %>% dplyr::select(locus) %>% 
  distinct(locus)
outlier_loci<- as.matrix(outlier_loci)

#Create outlier locus genlight and genind objects
genlight_out <-gl.keep.loc(genlight_2, loc.list = outlier_loci, verbose=0)
genind_out <- gl2gi(genlight_out)

#Create genpop of outlier loci, and calculate allele frequencies if desired
genpop_out<-genind2genpop(genind_out)
freq_out<-as.data.frame(t(makefreq(genpop_out, missing=0)))

#Check for loci out of HWE: No loci out of HWE
gl.filter.hwe(genlight_out, subset="each", n.pop.threshold = 3, alpha_val = 0.1)

##PCAs and Plots##
##Need the pca scores for LDA
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
ggsave(file="Figure##_out_PCAs_legend.svg", plot=legend_plot, width=7, height=0.5) 

#Calculate diversity metrics if desired
dartR_div<-gl.report.heterozygosity(genlight_out)

#Calculate pairwise Fst if desired
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

##LDA: Conduct linear discriminant analysis to test GSI

#First create random list of training and validation individuals, using half of each population
#These same delineations will be used in DAPC later
ctemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "CWFL"))), 0)))
dtemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "DRTO"))), 0)))
ntemp<-sample(rep(c("T", "V"), round((0.5*length(which(pca1.scores$Population == "NGM"))), 0)))
stemp<-sample(rep(c("T", "V"), (0.5*length(which(pca1.scores$Population == "SWFL"))))) #no round because odd number
extra<-sample (c("T","V"), 1) #add because SWFL odd
tv_list<-c(ctemp, dtemp, ntemp,stemp, extra)

#Format training and validation lists
pca1.scores$lda_status<-c(tv_list)
pca1.scores$lda_status<-as.factor(pca1.scores$lda_status)
training<-which(pca1.scores$lda_status == "T")
validation<-which(pca1.scores$lda_status == "V")
pca1.scores.lda<-pca1.scores[,-7]

#Conduct LDA
#First use training individuals to define axes
z <- lda(Population~.,pca1.scores.lda, subset=training)

#As a test, see how well training individuals reassign
z1<-predict(z,pca1.scores.lda[training,])

#If any individuals reassigning incorrectly, this will output "FALSE"
z1$class==pca1.scores.lda[training,]$Population

#Format dataframe of training individuals + reassignments
z1_lds<-as.data.frame(z1$x)
z1_lds$lda_status<-c("Training")
z1_lds$orig<-gsub('[[:digit:]]+',"",rownames(z1_lds))
z1_lds$post<-z1$class

#Now conduct assignment test with validation individuals
out <- predict(z,pca1.scores.lda[validation,])

#If any individuals reassigning incorrectly, this will output "FALSE"
out$class==pca1.scores.lda[validation,]$Population

#Format dataframe of validation individuals + reassignments
out_lds<-as.data.frame(out$x)
out_lds$lda_status<-c("Validation")
out_lds$orig<-gsub('[[:digit:]]+',"",rownames(out_lds))
out_lds$post<-out$class

#Make plots of LDA to vissualize assignments. 
#Transparent squares/dashed ellipses = training individuals
#Circles/solid ellipses = validation individuals

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

#Save figure
ldas<-grid.arrange(LD1_2, LD1_3, ncol=2, nrow=1)
ggsave(file="Figure##_LDAs.svg", plot=ldas, width=10, height=5)

#Save legend
legend<-ggpubr::get_legend(LD_legend)
legend_plot<-as_ggplot(legend)
ggsave(file="Figure##_LDAs_legend.svg", plot=legend_plot, width=7, height=0.5) 


## Conduct DAPC to test GSI

#Split genind by randomly assigned training and validation groups from above
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
set.seed(4)

#Calibrated DAPC
set.seed(5) 
dapc_a_score <- dapc(training_genind, n.pca=100, n.da=n.da) 
temp_score <- optim.a.score(dapc_a_score) #Determine best number of PCs to retain for DAPC analysis
n.pc <- temp_score$best #Save best number of PCs
dapc2 <-dapc(training_genind,training_genind@pop, n.pca = n.pc, n.da = n.da)
dapc2$IND <- row.names(dapc2$ind.coord)
dapc2_info <- as.data.frame(cbind(dapc2$IND, dapc2$ind.coord, grp_all$grp))

dapc_result <- scatter(dapc2, scree.pca = TRUE,
                       pch = 20, cell = 0, cstar = 1,
                       solid = 0.8, cex = 3, clab = 1)

##Make dataframe for training assignments
assignments<-as.data.frame(dapc2$posterior)
assignments$pop<-str_sub(dapc2$grp, 1,2)

assignments$ind<-rownames(assignments)
assignments$post<-colnames(assignments[,1:4])[max.col(assignments[,1:4],ties.method="first")]
assignments$train_val<-"train"

##Assign validation individuals
dapc_pred<-predict.dapc(dapc2, newdata=validation_genind)
dapc_pred$posterior

##Make dataframe for validation assignments
val_assignments<-as.data.frame(dapc_pred$posterior)
val_assignments$pop<-str_sub(row.names(dapc_pred$ind.scores), 1,2)
#assignments$pop<-gsub ("_", "", assignments$pop)
val_assignments$ind<-rownames(val_assignments)
val_assignments$post<-colnames(val_assignments[,1:4])[max.col(val_assignments[,1:4],ties.method="first")]
val_assignments$train_val<-"val"

#Format dataframes for bar plot and output barplot figure
whole_assign<-rbind(assignments, val_assignments)
assign_melt<-melt(whole_assign)

assign_melt<-assign_melt%>%
  arrange(pop) %>%
  arrange(ind)%>%
  arrange(train_val)
assign_ind_table <- table(assign_melt$ind)
ind_levels <- names(assign_ind_table)[order(assign_ind_table)]
assign_melt$ind <- factor(assign_melt$ind, levels = ind_levels)

c<-ggplot(assign_melt, 
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

compoplot_dapc2_legend<-ggplot(assign_melt, 
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
rubias_assign<-mix_est$indiv_posteriors

ggplot(rubias_assign, 
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
