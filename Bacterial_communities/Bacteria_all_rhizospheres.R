# # # Fungi analyses # # # 

#### Packages installations ####

# Required packages for taxonomy analysis
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
# install the necessary packages for beta-gamma diversity
# For windows system:
# install.packages("doMC", repos = "http://R-Forge.R-project.org")
# For linux or mac
# install.packages("doMC")
# Then install the following packages
# install.packages("lokern")
# install.packages("monomvn")
# install.packages("pspline")
# devtools::install_github('csb5/beem')
# install.packages("mice")

# Packages for MICROECO analysis
# List may not be complete, please see tutorial MICROECO
# From tutorial Microeco package (https://chiliubio.github.io/microeco_tutorial/)
# install.packages("microeco")
# allow more waiting time to download each package
# options(timeout = 1000)
# If a package is not installed, it will be installed from CRAN
# First select the packages of interest
# tmp <- c("microeco", "mecoturn", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "igraph", "picante", "pheatmap", "rgexf", 
# "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally", "ggraph", "networkD3", "poweRlaw", "ggtern", "SRS", "performance")
# Now check or install
# for(x in tmp){
#  if(!require(x, character.only = TRUE)) {
#   install.packages(x, dependencies = TRUE)
# }
# }
# install.packages("BiocManager")
# install.packages("file2meco", repos = BiocManager::repositories())
# install.packages("MicrobiomeStat", repos = BiocManager::repositories())
# install.packages("WGCNA", repos = BiocManager::repositories())
# BiocManager::install("ggtree")
# BiocManager::install("metagenomeSeq")
# BiocManager::install("ALDEx2")
# BiocManager::install("ANCOMBC")
# download link of the compressed packages archive
# Alternative from Gitee "https://gitee.com/chiliubio/microeco_dependence/releases/download/v0.20.0/microeco_dependence.zip"
# url <- "https://github.com/ChiLiubio/microeco_dependence/releases/download/v0.20.0/microeco_dependence.zip"
# allow more time to download the zip file in R
# options(timeout = 2000)
# Another way is to open the upper url in browser to download the zip file and move it to the current R working directory
# download.file(url = url, destfile = "microeco_dependence.zip")
# uncompress the file in R
# tmp <- "microeco_dependence"
# unzip(paste0(tmp, ".zip"))
# install devtools
# if(!require("devtools", character.only = TRUE)){install.packages("devtools", dependencies = TRUE)}
# run these one by one
# devtools::install_local(paste0(tmp, "/", "SpiecEasi-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "mixedCCA-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "SPRING-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "NetCoMi-main.zip"), repos = BiocManager::repositories())
# devtools::install_local(paste0(tmp, "/", "beem-static-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "chorddiag-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "ggradar-master.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "ggnested-main.zip"), dependencies = TRUE)
# devtools::install_local(paste0(tmp, "/", "ggcor-1-master.zip"), dependencies = TRUE)
#  Either seqinr or Biostrings package should be installed for reading and writing fasta file
# install.packages("seqinr", dependencies = TRUE)
# or install Biostrings from bioconductor https://bioconductor.org/packages/release/bioc/html/Biostrings.html
# see https://github.com/ChiLiubio/file2meco to install file2meco if needed
# Since two of NetCoMi's dependencies are only available on GitHub, 
# it is recommended to install them first:
# devtools::install_github("zdk123/SpiecEasi")
# devtools::install_github("GraceYoon/SPRING")
# Install NetCoMi
# devtools::install_github("stefpeschel/NetCoMi", 
# repos = c("https://cloud.r-project.org/",
#  BiocManager::repositories()))

# install_github("zdk123/SpiecEasi")



#### Load libraries for taxonomy ####

library('remotes')
library("dada2")
library('BiocManager')
library('phyloseq')
library('openssl')
library('theseus')
library('tidyverse')
library('Biostrings')
library("ggpubr")
library("dplyr")
library("dendextend")
library('microViz')
library("tidyr")
library("reshape2")
library("gridExtra")
library("scales")
library("parallel")
library("permute")
library("sunburstR")
library("htmltools")
library("lattice")
library("Rmisc")
library("knitr")
library("kableExtra")
library("data.table")
library("grid")
library("microbiome")
library("ape")
library("adespatial")
library("gtable")
library("Biostrings")
library("BAT")
library("car")
library("tibble")
library("cowplot")
library("RVAideMemoire")
library("ggtern")
library("bestNormalize")
library("readxl")
library("ggsci")
library("ggrepel")
library("ggplot2")
library("biohelper")
library("vegan")  
library("agricolae")  
library("eulerr")  
library("stats")  #  For Kruskal-Wallis test
library("viridis")
library("ggforce")
library("reshape")
library("mice")

#### Load libraries for microeco analysis ####

library('devtools')
library('GUniFrac') #For phylogenetic tree
library('MiscMetabar') #For phylogenetic tree
library('microeco')
library('file2meco')
library('ggalluvial') #For alluvional plots
library('magrittr') #For PCoA visualisation
library('aplot') #For PCoA visualisation
library('mecoturn')
library ('ggcor')
library('glmmTMB')
library('lmerTest')
library('WGCNA')
library('ggraph')
library("networkD3")
library ('chorddiag')
library('circlize')
library('NetCoMi')
library('SpiecEasi')
library("pheatmap")
library ('ggtree')
library('caret')
library('rfPermute')
library("Boruta")
library("parallel")
library("rsample") 
library("randomForest")  
library("gridExtra")
library("multiROC")
library("ranger")

#### Load libraries for other stat ####

library(metagMisc)

#### GitHub Reference Tutorial Microeco ####

# https://chiliubio.github.io/microeco_tutorial/
# http://search.r-project.org/CRAN/refmans/microeco/html/trans_env.html#method-trans_env-cal_mantel

#### Environment and objects preparation ####

# Load Data Set in your objects created when using DADA2
# Select the folder in which the taxonomic assignment is located
# Assign this location to path_results
path_results<-getwd()
# Create the other objects that will be filled later
seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
taxa=readRDS(paste0(path_results,"/taxaUNITE.rds"))
rownames(taxa) <-openssl::md5(rownames(taxa))
# Upload the table with the metadata
Metadata<-read.table("Fungi_Metadata.csv",h=T,sep=",",row.names = 1)
Metadata
# Save the sequences. Create a DNAStringSet from the ASVs
dna <- readDNAStringSet(paste0(path_results,"/uniqueSeqs.fasta"))
# Check the presence of the sequences in the dna file
dna
# Creation of Phyloseq object for the analysis
colnames(seqtab.nochim) <-md5(colnames(seqtab.nochim))
ps_run<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),sample_data(Metadata),tax_table(taxa),refseq(dna))
# Check the object
ps_run
# May emerge error if run too many things altogether
# Save the phyloseq object
saveRDS(ps_run, file = paste0(path_results,"/ps_run.rds"))
# Prepare the ASVs table for output  
otab = pstoveg_otu(ps_run) %>% t() %>% as.data.frame()
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table.txt"),quote=FALSE,sep="\t", row.names = FALSE)
# Save taxonomy table
taxa = cbind(ASVs = rownames(taxa), taxa)
write.table(taxa, paste0(path_results,"/tax_table.txt"),quote=FALSE,sep="\t", row.names = F)

# Optional (requires some time to load):
# Plot the table (before filtering the data). Fill can be changed into Family, Order, Species, ...
# plot_bar(ps_run, fill="Phylum", merge)

# Allocation of objects for analyses
asv_counts <- as.matrix(otu_table(ps_run))
total_sequence_reads <- sum(asv_counts)
# Check
total_sequence_reads
num_asvs <- ntaxa(ps_run)
# Check
num_asvs


#### Addressing contamination in experimental design ####

# The microDecon option is a wrapper of the decon function from microDecon R package. 
# It performs the decon function on a phyloseq object
# with sample data and returns a decontaminated phyloseq object with sample data,
# taxonomy and reference sequences if present.
# The 'max_v' option subtract read associated to putative contaminant ASV by using
# their max read count in blank(s).I use microDecon because more accurate

# Apply decontamination method
ps_DE<-ps_decon(ps_run,method = "microDecon")

# Create dataset only with samples, no blanks
ps_Final <- subset_samples(ps_DE, amplicon_type == "sample")
# Check the object
summarize_phyloseq(ps_Final)
# New counts allocation
FinalASV_counts <- as.matrix(otu_table(ps_Final))
final_sequence_reads <- sum(FinalASV_counts)
final_sequence_reads
final_num_asvs <- ntaxa(ps_Final)
final_num_asvs

#### Subset creation for different projects and studies ####
# Create subsets, including only the sequences related to the samples, no controls or blanks,
# even if retained from the passage above
# OTU table, Taxa table, and others are stored in the phyloseq object
# Total
ps_Final <- subset_samples(ps_DE, amplicon_type == "sample")
# Different landcover types (Main project (MPJ) for the thesis)
ps_MainPJ <- subset_samples(ps_DE, ProjectFocus == "Field")
# Rhizosphere (RZ)
ps_Rhizosphere <- subset_samples(ps_DE, ProjectFocus == "Rhizosphere")

#### Processing sample-rarefying for each of the working dataset #### 
# This step is done to ensure that the samples are comparable, having equal readings-depth
# (done for Main Project)
# orders samples from highest to lowest
sample_sums(ps_Rhizosphere)[order(sample_sums(ps_Rhizosphere))]
# Observe the number of reads, retain the minimum number similar to others (eg. if present 
# 18-490-23897-23786-[], retain 23897).
# We see that the sample 221 and 186 have really low reads. To discard. They are from R1SB (just 18 reads), and 
# from R1SA (481 reads)
# The lowest number after them is 24803, which I will retain as minimum
# To display the plots to recognize better were to cut:
asv_before_MainPJ <- as(otu_table(ps_Rhizosphere), "matrix")
out<-vegan::rarecurve(asv_before_MainPJ, step=100,lwd=2, ylab="ASV Richness", xlab="Sequence Sample Size", main="INSERT_GENE_TARGET rRNA", label=F)
# Set seed to reproduce the data, since the rarefaction will sample
set.seed(100)
# Rarefy the samples without replacement. 
# Rarefaction is used to simulate an even number of reads per sample. 
ps_rarefiedMPJ <- rarefy_even_depth(ps_Rhizosphere, sample.size =   24803 , replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# Allocation of objects after rarefying
asv_counts_rarefiedMPJ <- as.matrix(otu_table(ps_rarefiedMPJ))
After_rar_sequence_readsMPJ <- sum(asv_counts_rarefiedMPJ)
After_rar_sequence_readsMPJ
# Taxa
after_rar_num_asvsMPJ<- ntaxa(ps_rarefiedMPJ)
after_rar_num_asvsMPJ
# Check
summarize_phyloseq(ps_rarefiedMPJ)
ps_rarefiedMPJ
# Check the structure
str(sample_data(ps_rarefiedMPJ))
# Dataset to use from now:
ps_rarefiedMPJ

#### Summary of the list of taxa ####

# (for Main Project)
# List of taxonomic ranks
rank_names(ps_rarefiedMPJ)
taxonomic_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# Initialize a list to store counts for each rank
taxa_counts <- list()
# Loop through the taxonomic ranks and calculate counts
for (rank in taxonomic_ranks) {
  glommed_ps <- ps_rarefiedMPJ %>%
    tax_glom(taxrank = rank)
  taxa_count <- ntaxa(glommed_ps)
  taxa_counts[[rank]] <- taxa_count
}
#  Display the counts for each taxonomic rank
for (rank in taxonomic_ranks) {
  cat("Number of", rank, ":", taxa_counts[[rank]], "\n")
}

#### MICROECO Analysis ####

# Convert the dataset to work in microeco, convert the phyloseq file
meco_datasetMPJ <- phyloseq2meco(ps_rarefiedMPJ)
# check
meco_datasetMPJ
# remove ASVs which are not assigned in the Kingdom of interest (even if not expected), and should not be present
# use R subset function to filter taxa in tax_table
meco_datasetMPJ$tax_table %<>% base::subset(Kingdom == "k__Fungi")
# another way with grepl function
# meco_datasetMPJ$tax_table %<>% .[grepl("Fungi", .$Kingdom), ]
# Check the result
meco_datasetMPJ
# If I want to remove others:
# meco_datasetMPJ$filter_pollution(taxa = c("INSERTNAME", "E.G.-->chloroplast"))
# To make the ASVs and sample information consistent across all files in the object:
meco_datasetMPJ$tidy_dataset()
# check
meco_datasetMPJ
# sample_sums() to check the sequence numbers in each sample
meco_datasetMPJ$sample_sums() %>% range
# function save_table can be performed to save all the basic data in microtable object to local files, including feature abundance, metadata, taxonomic table, phylogenetic tree and representative sequences.
meco_datasetMPJ$save_table(dirpath = "basic_files", sep = ",")

#### Taxa abundance calculation ####
# Calculate the taxa abundance at each taxonomic rank using cal_abund(). 
# This function generate a list called taxa_abund stored in the microtable object. 
# This list contains several data frame of the abundance information at each taxonomic rank.
# default parameter (rel = TRUE) denotes relative abundance
meco_datasetMPJ$cal_abund()
# return taxa_abund list in the object
class(meco_datasetMPJ$taxa_abund)
# show part of the relative abundance at Phylum level
meco_datasetMPJ$taxa_abund$Phylum[1:5, 1:5]
# The function save_abund() can be used to save the taxa abundance file to a local place easily.
meco_datasetMPJ$save_abund(dirpath = "taxa_abund")

#### Dataset manipulation ####

# To group or subset see https://chiliubio.github.io/microeco_tutorial/basic-class.html
# Paragraph 3.1.4 and beyond
# See also "The filter_taxa function can be applied to filter the features with 
# low abundance or occurrence frequency when needed. 
# For other operations on the features, please directly manipulate the otu_table of your microtable object."

#### Alpha diversity calculation ####
# To add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow, and the tree has to be provided
# If not specified, different indexes of diversity are calculated
meco_datasetMPJ$cal_alphadiv(PD = F)
# return alpha_diversity in the object
class(meco_datasetMPJ$alpha_diversity)
# save alpha_diversity to a directory
meco_datasetMPJ$save_alphadiv(dirpath = "alpha_diversity")
# See levels
meco_datasetMPJ$sample_table$Landcover
# Relevel
meco_datasetMPJ$sample_table$Landcover <- factor(meco_datasetMPJ$sample_table$Landcover, levels = c("Grassland", "TwoYearsOld", "TenYearsold",
                                                                                                    "TwentyfourYearsold", "Remnant")) 
# Creating a trans_alpha object can return two data.frame with the prefix ‘data_’: 
# data_alpha and data_stat. data_alpha is used for subsequent differential test 
# and visualization.
t1 <- trans_alpha$new(dataset = meco_datasetMPJ, group = "Landcover")
t1$data_alpha
# Save table alpha diversity with ASV counts, number of ASVs for each sample
write.csv(t1$data_alpha, "alpha_landcover.csv", row.names = TRUE)
# return t1$data_stat
head(t1$data_stat)
t1$data_stat
# Then, we test the differences among groups using Kruskal-Wallis Rank Sum Test (overall test when groups > 2),
# Wilcoxon Rank Sum Tests (for paired groups), 
# Dunn’s Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2) and Anova with multiple comparisons.
t1$cal_diff(method = "KW")
#  return t1$res_diff
head(t1$res_diff)
t1$cal_diff(method = "KW_dunn")
#  return t1$res_diff
head(t1$res_diff)
# The result is stored in object$res_diff
# more options
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
# head(t1$res_diff)
t1$cal_diff(method = "wilcox")
head(t1$res_diff)
# t1$cal_diff(method = "t.test")
# head(t1$res_diff)
# Try anova
t1$cal_diff(method = "anova")
t1$data_stat
t1$plot_alpha(measure = "Observed", shape = "Landcover")
# Further Post hoc analyses to implement if significant results from ANOVA,
# for e.g. duncan.test
# The multi-factor analysis of variance is also supported with the formula parameter, such as two-way anova.
t1 <- trans_alpha$new(dataset = meco_datasetMPJ, group = "Landcover")
t1$cal_diff(method = "anova", formula = "Landcover+Hill__side")
head(t1$res_diff)
# see the help document for the usage of formula

#### Hills numbers calculation ####
phyloseq_inext(ps_rarefiedMPJ)

#### Alpha diversity calculation at higher taxonomic level ####
# Aggregate to Genus level
t1 <- trans_alpha$new(dataset = meco_datasetMPJ$merge_taxa("Genus"), group = "Landcover")
t1$data_alpha
# Try anova
t1$cal_diff(method = "anova")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "Landcover")
# Aggregate to Family level
t1 <- trans_alpha$new(dataset = meco_datasetMPJ$merge_taxa("Family"), group = "Landcover")
t1$data_alpha
# Try anova
t1$cal_diff(method = "anova")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "Landcover")
# Aggregate to Order level
t1 <- trans_alpha$new(dataset = meco_datasetMPJ$merge_taxa("Order"), group = "Landcover")
t1$data_alpha
# Try anova
t1$cal_diff(method = "anova")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "Landcover")
# Aggregate to Class level
t1 <- trans_alpha$new(dataset = meco_datasetMPJ$merge_taxa("Class"), group = "Landcover")
t1$data_alpha
# Try anova
t1$cal_diff(method = "anova")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "Landcover")
# Aggregate to Phylum level
t1 <- trans_alpha$new(dataset = meco_datasetMPJ$merge_taxa("Phylum"), group = "Landcover")
t1$data_alpha
# Try anova
t1$cal_diff(method = "anova")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "Landcover")

#### Alpha diversity within each landcover type ####
# ASVs
# Grassland
# remember first clone the whole dataset
group_Grass <- clone(meco_datasetMPJ)
# select 'Grassland'
group_Grass$sample_table <- subset(group_Grass$sample_table, Landcover == "Grassland")
# trim all the data
group_Grass$tidy_dataset()
# Check
group_Grass
# Extract numeric parts (transect) and create grouping factor
group_Grass$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_Grass$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_Grass$sample_table$Numeric_Group[is.na(group_Grass$sample_table$Numeric_Group)] <- 0
group_Grass$sample_table$Numeric_Group <- factor(group_Grass$sample_table$Numeric_Group)
# Check
group_Grass$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_Grass<- trans_alpha$new(dataset = group_Grass, group = "Numeric_Group")
# See
group_Grass$data_alpha
# Try anova
group_Grass$cal_diff(method = "anova")
group_Grass$res_diff
# Plot
group_Grass$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# Two Years Old
# remember first clone the whole dataset
group_2YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_2YO$sample_table <- subset(group_2YO$sample_table, Landcover == "TwoYearsOld")
# trim all the data
group_2YO$tidy_dataset()
# Check
group_2YO
# Extract numeric parts (transect) and create grouping factor
group_2YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_2YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_2YO$sample_table$Numeric_Group[is.na(group_2YO$sample_table$Numeric_Group)] <- 0
group_2YO$sample_table$Numeric_Group <- factor(group_2YO$sample_table$Numeric_Group)
# Check
group_2YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_2YO<- trans_alpha$new(dataset = group_2YO, group = "Numeric_Group")
# See
group_2YO$data_alpha
# Try anova
group_2YO$cal_diff(method = "anova")
group_2YO$res_diff
# Plot
group_2YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")


# Ten Years Old
# remember first clone the whole dataset
group_10YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_10YO$sample_table <- subset(group_10YO$sample_table, Landcover == "TenYearsold")
# trim all the data
group_10YO$tidy_dataset()
# Check
group_10YO
# Extract numeric parts (transect) and create grouping factor
group_10YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_10YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_10YO$sample_table$Numeric_Group[is.na(group_10YO$sample_table$Numeric_Group)] <- 0
group_10YO$sample_table$Numeric_Group <- factor(group_10YO$sample_table$Numeric_Group)
# Check
group_10YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_10YO<- trans_alpha$new(dataset = group_10YO, group = "Numeric_Group")
# See
group_10YO$data_alpha
# Try anova
group_10YO$cal_diff(method = "anova")
group_10YO$res_diff
# Plot
group_10YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# TwentyfourYearsold
# remember first clone the whole dataset
group_24YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_24YO$sample_table <- subset(group_24YO$sample_table, Landcover == "TwentyfourYearsold")
# trim all the data
group_24YO$tidy_dataset()
# Check
group_24YO
# Extract numeric parts (transect) and create grouping factor
group_24YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_24YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_24YO$sample_table$Numeric_Group[is.na(group_24YO$sample_table$Numeric_Group)] <- 0
group_24YO$sample_table$Numeric_Group <- factor(group_24YO$sample_table$Numeric_Group)
# Check
group_24YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_24YO<- trans_alpha$new(dataset = group_24YO, group = "Numeric_Group")
# See
group_24YO$data_alpha
# Try anova
group_24YO$cal_diff(method = "anova")
group_24YO$res_diff
# Plot
group_24YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# Remnant
# remember first clone the whole dataset
group_Remnant <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_Remnant$sample_table <- subset(group_Remnant$sample_table, Landcover == "Remnant")
# trim all the data
group_Remnant$tidy_dataset()
# Check
group_Remnant
# Extract numeric parts (transect) and create grouping factor
group_Remnant$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_Remnant$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_Remnant$sample_table$Numeric_Group[is.na(group_Remnant$sample_table$Numeric_Group)] <- 0
group_Remnant$sample_table$Numeric_Group <- factor(group_Remnant$sample_table$Numeric_Group)
# Check
group_Remnant$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_Remnant<- trans_alpha$new(dataset = group_Remnant, group = "Numeric_Group")
# See
group_Remnant$data_alpha
# Try anova
group_Remnant$cal_diff(method = "anova")
group_Remnant$res_diff
# Plot
group_Remnant$plot_alpha(measure = "Observed", shape = "Numeric_Group")

#### Alpha diversity Phylum level for each landcover type ####
# Grassland
# remember first clone the whole dataset
group_Grass <- clone(meco_datasetMPJ)
# select 'Grassland'
group_Grass$sample_table <- subset(group_Grass$sample_table, Landcover == "Grassland")
# trim all the data
group_Grass$tidy_dataset()
# Check
group_Grass
# Extract numeric parts (transect) and create grouping factor
group_Grass$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_Grass$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_Grass$sample_table$Numeric_Group[is.na(group_Grass$sample_table$Numeric_Group)] <- 0
group_Grass$sample_table$Numeric_Group <- factor(group_Grass$sample_table$Numeric_Group)
# Check
group_Grass$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_Grass<- trans_alpha$new(dataset = group_Grass$merge_taxa('Phylum'), group = "Numeric_Group")
# See
group_Grass$data_alpha
# Try anova
group_Grass$cal_diff(method = "anova")
head(group_Grass$res_diff)
# Plot
group_Grass$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# Two Years Old
# remember first clone the whole dataset
group_2YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_2YO$sample_table <- subset(group_2YO$sample_table, Landcover == "TwoYearsOld")
# trim all the data
group_2YO$tidy_dataset()
# Check
group_2YO
# Extract numeric parts (transect) and create grouping factor
group_2YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_2YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_2YO$sample_table$Numeric_Group[is.na(group_2YO$sample_table$Numeric_Group)] <- 0
group_2YO$sample_table$Numeric_Group <- factor(group_2YO$sample_table$Numeric_Group)
# Check
group_2YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_2YO<- trans_alpha$new(dataset = group_2YO$merge_taxa('Phylum'), group = "Numeric_Group")
# See
group_2YO$data_alpha
# Try anova
group_2YO$cal_diff(method = "anova")
head(group_2YO$res_diff)
# Plot
group_2YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# Ten Years Old
# remember first clone the whole dataset
group_10YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_10YO$sample_table <- subset(group_10YO$sample_table, Landcover == "TenYearsold")
# trim all the data
group_10YO$tidy_dataset()
# Check
group_10YO
# Extract numeric parts (transect) and create grouping factor
group_10YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_10YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_10YO$sample_table$Numeric_Group[is.na(group_10YO$sample_table$Numeric_Group)] <- 0
group_10YO$sample_table$Numeric_Group <- factor(group_10YO$sample_table$Numeric_Group)
# Check
group_10YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_10YO<- trans_alpha$new(dataset = group_10YO$merge_taxa('Phylum'), group = "Numeric_Group")
# See
group_10YO$data_alpha
# Try anova
group_10YO$cal_diff(method = "anova")
head(group_10YO$res_diff)
# Plot
group_10YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# TwentyfourYearsold
# remember first clone the whole dataset
group_24YO <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_24YO$sample_table <- subset(group_24YO$sample_table, Landcover == "TwentyfourYearsold")
# trim all the data
group_24YO$tidy_dataset()
# Check
group_24YO
# Extract numeric parts (transect) and create grouping factor
group_24YO$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_24YO$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_24YO$sample_table$Numeric_Group[is.na(group_24YO$sample_table$Numeric_Group)] <- 0
group_24YO$sample_table$Numeric_Group <- factor(group_24YO$sample_table$Numeric_Group)
# Check
group_24YO$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_24YO<- trans_alpha$new(dataset = group_24YO$merge_taxa('Phylum'), group = "Numeric_Group")
# See
group_24YO$data_alpha
# Try anova
group_24YO$cal_diff(method = "anova")
head(group_24YO$res_diff)
# Plot
group_24YO$plot_alpha(measure = "Observed", shape = "Numeric_Group")

# Remnant
# remember first clone the whole dataset
group_Remnant <- clone(meco_datasetMPJ)
# select 'Two Years Old'
group_Remnant$sample_table <- subset(group_Remnant$sample_table, Landcover == "Remnant")
# trim all the data
group_Remnant$tidy_dataset()
# Check
group_Remnant
# Extract numeric parts (transect) and create grouping factor
group_Remnant$sample_table$Numeric_Group <- as.numeric(stringr::str_extract(group_Remnant$sample_table$original_sample_id, "(?<=\\D)(\\d+)"))
group_Remnant$sample_table$Numeric_Group[is.na(group_Remnant$sample_table$Numeric_Group)] <- 0
group_Remnant$sample_table$Numeric_Group <- factor(group_Remnant$sample_table$Numeric_Group)
# Check
group_Remnant$sample_table$Numeric_Group
# Calculate alpha diversity at ASVs level
group_Remnant<- trans_alpha$new(dataset = group_Remnant$merge_taxa('Phylum'), group = "Numeric_Group")
# See
group_Remnant$data_alpha
# Try anova
group_Remnant$cal_diff(method = "anova")
head(group_Remnant$res_diff)
# Plot
group_Remnant$plot_alpha(measure = "Observed", shape = "Numeric_Group")

#### Plot Alpha diversity ####
# The plot_alpha function add the significance label by searching the results in object$res_diff 
# instead of recalculating the significance. 
# Plot the alpha diversity for each group and include the anova result:
t1$cal_diff(method = "anova")
#  y_increase can adjust the distance from the letters to the highest point, see for example:
t1$plot_alpha(measure = "Chao1", y_increase = 0.3)
t1$plot_alpha(measure = "Chao1", y_increase = 0.1)
#  add_sig_text_size: letter size adjustment
t1$plot_alpha(measure = "Chao1", add_sig_text_size = 6, add = "jitter", order_x_mean = TRUE)
# Statistical difference from another statistical analysis:
t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Observed", shape = "Landcover")
#  y_start: starting height for the first label
#  y_increase: increased height for each label
t1$plot_alpha(measure = "Chao1", shape = "Landcover", add = "jitter", y_start = 0.1, y_increase = 0.1)
# Possibility to remove the 'ns' with:
t1$res_diff %<>% base::subset(Significance != "ns")
t1$plot_alpha(measure = "Chao1", add = "dotplot", xtext_size = 15)
# The trans_alpha class supports the differential test of groups within each group 
# using the by_group parameter.
t1 <- trans_alpha$new(dataset = meco_datasetMPJ, group = "Hill__side", by_group = "Landcover")
t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Shannon")
# Note: Scheirer Ray Hare test is a nonparametric test that is suitable for a two-way factorial experiment.
# t1 <- trans_alpha$new(dataset = meco_datasetMPJ)
# require rcompanion package to be installed
# t1$cal_diff(method = "scheirerRayHare", formula = "Landcover+Hill__side")

# Additional code for taxonomic plot
# Abundance for each landcover type
# try<-trans_abund$new(dataset=tmp, ntaxa = 5)
# try$plot_bar(facet = "Landcover")
# try$plot_heatmap(facet = "Landcover")

#### Alpha diversity and lmer ####
# Linear mixed-effects model can be selected with the method = "lme". 
# This model is implemented based on the lmerTest package. For more parameters, 
# please see lmerTest::lmer function. Please use parameter passing when more parameters are needed. 
# For the formula usage, please follow this (https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html). 
# In the return table, conditional R2 is the total variance explained by fixed and random effects,
# and marginal R2 is the variance explained by fixed effects.
# if(!require("lmerTest")) install.packages("lmerTest")
t1 <- trans_alpha$new(dataset = meco_datasetMPJ)
t1$cal_diff(method = "lme", formula = "Landcover + (1|Hill__side)")
View(t1$res_diff)
# return_model = TRUE can return original models, i.e. object$res_model
t1$cal_diff(method = "lme", formula = "Landcover + (1|Hill__side)", return_model = TRUE)
# Note that from v1.9.0, the parameter plot_type control which type of plot is employed. 
# All the options starting with “gg” (e.g., “ggboxplot”, “ggdotplot”, “ggviolin”, “ggstripchart”, “ggerrorplot”) 
# means they are the functions coming from the ggpubr package.
# t1 <- trans_alpha$new(dataset = mt, group = "Type")
# t1$cal_diff(method = "KW_dunn", KW_dunn_letter = TRUE)
# t1$plot_alpha(plot_type = "ggboxplot", add = "none")

# t1$plot_alpha(plot_type = "ggdotplot")
# t1$plot_alpha(plot_type = "ggdotplot", fill = "Type", alpha = 0.3)
# t1$plot_alpha(plot_type = "ggdotplot", add = "mean_se")
# t1$plot_alpha(plot_type = "ggdotplot", add = c("mean_se", "violin"))
# t1$plot_alpha(plot_type = "ggdotplot", add = c("mean_se", "violin"), fill = "Type", alpha = 0.2)

# t1$plot_alpha(plot_type = "ggviolin")
# t1$plot_alpha(plot_type = "ggviolin", y_increase = 0.4, add = "mean_se")
# t1$plot_alpha(plot_type = "ggviolin", fill = "Type", alpha = 0.2, y_increase = 0.4, add = "mean_se", add_sig_text_size = 6)

# t1$plot_alpha(plot_type = "ggstripchart", add = "mean_se")

# t1$plot_alpha(plot_type = "ggerrorplot")

# The option "errorbar" or "barerrorbar" in plot_alpha will invoke the data_stat instead of data_alpha for the Mean±SE (or SD) plot (based on the ggplot2). 
# The line is optional to be added between points (Mean) for the case with a gradient.
# t1 <- trans_alpha$new(dataset = mt, group = "Group")
# t1$cal_diff(method = "KW_dunn", measure = "Chao1", KW_dunn_letter = TRUE)
# t1$plot_alpha(measure = "Chao1")
# t1$plot_alpha(plot_type = "errorbar", measure = "Chao1")
# t1$plot_alpha(plot_type = "errorbar", measure = "Chao1", y_increase = -0.2)
# t1$plot_alpha(plot_type = "errorbar", measure = "Chao1", y_increase = -0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1)
# t1$plot_alpha(plot_type = "errorbar", plot_SE = FALSE, measure = "Chao1", y_increase = 0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1)
# t1$plot_alpha(plot_type = "barerrorbar", measure = "Chao1")
# t1$plot_alpha(plot_type = "barerrorbar", measure = "Chao1", y_increase = -0.3)
# t1$plot_alpha(plot_type = "barerrorbar", measure = "Chao1", bar_width = 0.6, errorbar_width = 0.2, errorbar_size = 1, errorbar_addpoint = FALSE)
# by_group example
# t1 <- trans_alpha$new(dataset = mt, group = "Landcover", by_group = "Landcover")
# t1$cal_diff(method = "KW_dunn", measure = "Shannon", KW_dunn_letter = TRUE)
# View(t1$res_diff)
# t1$plot_alpha(plot_type = "errorbar", measure = "Shannon")
# t1$plot_alpha(plot_type = "errorbar", measure = "Shannon", add_line = TRUE, line_type = 2)
# t1$plot_alpha(plot_type = "errorbar", plot_SE = FALSE, measure = "Shannon", add_line = TRUE, line_type = 2)
# From v1.4.0, the heatmap can be used to visualize the significances for the case with multiple factors in the formula.
# t1 <- trans_alpha$new(dataset = meco_datasetMPJ, group = "Landcover")
# t1$cal_diff(method = "anova", formula = "Group+Type+Group:Type")
# t1$plot_alpha(color_palette = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), trans = "log10")
# t1$plot_alpha(color_palette = c("#053061", "white", "#A50026"), trans = "log10")
# t1$plot_alpha(color_values = c("#053061", "white", "#A50026"), trans = "log10")
# t1$plot_alpha(color_values = c("#053061", "white", "#A50026"), trans = "log10", filter_feature = "", text_y_position = "left")
# t1$plot_alpha(color_values = c("#053061", "white", "#A50026"), trans = "log10", filter_feature = "", text_y_position = "left", cluster_ggplot = "row")

#### Beta diversity calculation ####
# If method parameter is not provided, the function automatically calculates Bray-curtis, Jaccard, 
# weighted Unifrac and unweighted unifrac matrixes (Lozupone and Knight 2005).
# unifrac = FALSE means do not calculate unifrac metric (which is based on phylogeny)
# Below the function to build the phylogenetic tree diversity for the core community, excluding the taxa
# that have less than 100 sequences overall. Thus analyses of the core community.
# Phylogenetic tree construction. See https://adrientaudiere.github.io/MiscMetabar/reference/build_phytree_pq.html
# if (requireNamespace("phangorn")) {
# set.seed(22)
#   df <- subset_taxa_pq(ps_rarefiedMPJ, taxa_sums(ps_rarefiedMPJ) > 100)
#   df_tree <- build_phytree_pq(df, nb_bootstrap = 10)
#   plot(df_tree$UPGMA)
#   phangorn::plotBS(df_tree$UPGMA, df_tree$UPGMA_bs, main = "UPGMA")
#   plot(df_tree$NJ, "unrooted")
#   plot(df_tree$ML)
#   phangorn::plotBS(df_tree$ML$tree, df_tree$ML_bs, p = 20, frame = "circle")
#   phangorn::plotBS(
#     df_tree$ML$tree,
#     df_tree$ML_bs,
#     p = 20,
#     frame = "circle",
#     method = "TBE"
#   )
#   plot(phangorn::consensusNet(df_tree$ML_bs))
#  plot(phangorn::consensusNet(df_tree$NJ_bs))
#   ps_tree <- merge_phyloseq(df, df_tree$ML$tree)
# }
# Calculation without Phylogenetic tree, if unifrac=T it will use
# the phylogenetic tree, but note that the tree is built for the core community (taxa with>100 sequences)
meco_datasetMPJ$cal_betadiv(unifrac = F)
#  return beta_diversity list in the object
class(meco_datasetMPJ$beta_diversity)
#  save beta_diversity to a directory
meco_datasetMPJ$save_betadiv(dirpath = "beta_diversity")
# create trans_beta object
# For PCoA and NMDS, measure parameter must be provided.
# measure parameter should be either one of names(mt$beta_diversity) or a customized symmetric matrix
t1 <- trans_beta$new(dataset = meco_datasetMPJ, group = "Landcover", measure = "bray")
# Visualisation of the beta diversity
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
t1$res_ordination
# plot the PCoA result with confidence ellipse
t1$sample_table$Hill__side<-as.factor(t1$sample_table$Hill__side)
# t1$plot_ordination(plot_color = "Landcover", plot_type = c("point","ellipse")  , plot_shape = "Hill__side")
plot<-t1$plot_ordination(plot_color = "Landcover", plot_type = c("point","ellipse"))
plot + labs(title="Fungi communities")
# The warning appears because it says that the ellipses are draw using all the points, 
# not taking into account the N and S side of the hill, and thus the shape of them.
# Other examples and options
t1$plot_ordination(plot_color = "Landcover", plot_type = "point")
t1$plot_ordination(plot_color = "Landcover", point_size = 5, point_alpha = .2, plot_type = c("point", "ellipse"), ellipse_chull_fill = FALSE)
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("point", "centroid"))
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("point", "ellipse", "centroid"))
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("point", "chull"))
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("point", "chull", "centroid"))
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("chull", "centroid"))
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = c("point", "chull", "centroid"), add_sample_label = "SampleID")
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = "centroid")
t1$plot_ordination(plot_color = "Landcover", plot_shape = "Hill__side", plot_type = "centroid", centroid_segment_alpha = 0.9, centroid_segment_size = 1, centroid_segment_linetype = 1)
t1$plot_ordination(plot_type = c("point", "centroid"), plot_color = "Landcover", centroid_segment_linetype = 1)
t1$plot_ordination(plot_color = "Landcover", point_size = 5, point_alpha = 2, plot_type = c("point", "chull"), ellipse_chull_fill = FALSE, ellipse_chull_alpha = 0.1)
t1$plot_ordination(plot_color = "Landcover") + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)
# PCA with Genus data
tmp <- meco_datasetMPJ$merge_taxa(taxa = "Genus")
tmp$tax_table %<>% .[.$Genus != "g__", ]
tmp$tidy_dataset()
rownames(tmp$otu_table) <- tmp$tax_table[rownames(tmp$otu_table), "Genus"]
rownames(tmp$tax_table) <- tmp$tax_table[, "Genus"]
# Plot and compare between distances
t1 <- trans_beta$new(dataset = tmp)
t1$cal_ordination(method = "PCA")
t1$plot_ordination(plot_color = "Landcover", loading_arrow = TRUE, loading_text_italic = TRUE)
# Return exploring the previous general dataset
t1 <- trans_beta$new(dataset = meco_datasetMPJ, group = "Landcover", measure = "bray")
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(add = "mean")
# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "wilcox")
# parameters in plot_group_distance function will be passed to the plot_alpha function of trans_alpha class
t1$plot_group_distance(plot_type = "ggviolin", add = "mean_se")
t1$plot_group_distance(add = "mean")
# Clustering
tmp <- clone(meco_datasetMPJ)
# extract a part of data if needed
# tmp$sample_table %<>% subset(Landcover %in% c("Grassland", "Remnant"))
# tmp$tidy_dataset()
t1 <- trans_beta$new(dataset = tmp, group = "Landcover")
# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Landcover", replace_name = c("Landcover"))
# Permanova analysis to test beta diversity
t1 <- trans_beta$new(dataset = meco_datasetMPJ, group = "Landcover", measure = "bray")
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# manova for each paired groups, using FALSE (the parameter "by_group=" can be used to restrict the comparison between certain levels)
# see other specific parameters that can be used in this analysis in microeco tutorial
t1$cal_manova(manova_all = FALSE)
t1$res_manova
# Anosim analysis
t1$cal_anosim(group = "Landcover")
t1$res_anosim
t1$cal_anosim(group = "Landcover", paired = TRUE)
t1$res_anosim
# Extra analysis:
# PERMDISP(Anderson et al. 2011) is implemented to test multivariate homogeneity of groups dispersions (variances) based on the betadisper function of vegan package.
# for the whole comparison and for each paired groups
# t1$cal_betadisper()

#### Composition-based class exploration with plots ####
# Composition-based class. These analyses are to visualise the taxonomic abundance, considering both the different ASVs,
# and the number of reads for each of them.
# The trans_abund class and trans_venn class are organised into the section ‘Composition-based class’, 
# since they are mainly used to show the composition information of communities.
# create trans_abund object
# Phyla level (ntaxa=): select top 8 abundant Phyla.
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Phylum", ntaxa = 8)
# t1 object now include the transformed abundance data t1$abund_data 
# and other elements for the following plotting
# Plots
# Adjusting ladncover factors levels in the right order
t1$plot_bar(others_color = "grey70", facet = "Landcover", xtext_keep = FALSE, legend_text_italic = FALSE)
#  return a ggplot2 object
#  require package ggh4x, first run install.packages("ggh4x") if not installed
t1$plot_bar(others_color = "grey70", facet = c("Landcover", "Hill__side"), xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
# Genus level
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Genus", ntaxa = 8)
t1$plot_bar(others_color = "grey70", facet = "Landcover", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
t1$plot_bar(others_color = "grey70", facet = c("Landcover", "Hill__side"), xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
# use_alluvium = TRUE make the alluvial plot, clustering =TRUE can be used 
# to reorder the samples by clustering. bar_type = FALSE can remove 'others'
t1$plot_bar(bar_full = FALSE, use_alluvium = TRUE, clustering = TRUE, xtext_angle = 30, xtext_size = 3, color_values = RColorBrewer::brewer.pal(8, "Set2"))
# Barplot with mean values
# The bar plot can also be performed with group mean values. Note that, from v0.16.0, the parameter group_morestats = TRUE can be used to add more summary statistics in the return data_abund when groupmean parameter is provided.
# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Phylum", ntaxa = 10, groupmean = "Landcover")
g1 <- t1$plot_bar(order_x = c("Grassland", "TwoYearsOld", "TenYearsold", "TwentyfourYearsold", "Remnant"), legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
# show 15 taxa at Genus level
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Genus", ntaxa = 15)
t1$plot_box(group = "Landcover", xtext_angle = 30)
# Heatmap
# Show the heatmap with the high abundant genera.
# show 40 taxa at Genus level
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Genus", ntaxa = 40)
g1 <- t1$plot_heatmap(facet = "Landcover", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
#Show the plot
g1
g1 + theme(axis.text.y = element_text(face = 'italic'))
# Donut chart
# The donut and radar charts are implemented from v0.17.0. Please install the dependent packages according to the steps (https://chiliubio.github.io/microeco_tutorial/intro.html# dependence).
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Phylum", ntaxa = 8, groupmean = "Landcover")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)
# The coord_flip parameter in plot_bar function can be changed to make the coordinate 
# axis flipped. The clustering plot can also be added in the bar plot. 
# In this case, the coordinate axis will be flipped automatically for better visualization.
t1 <- trans_abund$new(dataset = meco_datasetMPJ, taxrank = "Phylum", ntaxa = 4, groupmean = "Landcover")
g1 <- t1$plot_bar(coord_flip = TRUE)
g1 <- g1 + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
# Show the plot
g1
g1 <- t1$plot_bar(clustering_plot = TRUE)
g1
# In this case, g1 (aplot object) is the combination of different ggplot objects
# to adjust the main plot, please select g1[[1]]
g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
#  save the figure (WRITE THE DISTANCE METRIC USED)
ggsave("ClusteringBCdistance.png", g1, width = 8, height = 5)

#### Venn analysis ####
# Analysis of the shared and unique taxa among sites
# The trans_venn class is developed for venn analysis, i.e. shared and unique taxa across samples/groups.
# First merge samples as one community for each group
tmp <- meco_datasetMPJ$merge_samples("Landcover")
# tmp is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(tmp, ratio = NULL)
t1$data_summary
t1$plot_venn()
# create venn plot with more information
t1 <- trans_venn$new(tmp, ratio = "seqratio")
t1$plot_venn()
# The integer is ASVs number, while the percentage data is the sequence number/total sequence number
# Now, we transform the results of venn plot to the traditional feature-sample table, that is, another object of microtable class
# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
tmp <- meco_datasetMPJ$merge_samples("Landcover")
t1 <- trans_venn$new(tmp)
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)
# Calculate taxa abundance, that is, the frequency
t2$cal_abund()
#  transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Phylum", ntaxa = 8)
t3$plot_bar(bar_full = F, legend_text_italic = T, xtext_angle = 30, 
            color_values = RColorBrewer::brewer.pal(8, "Set2"),ylab("Frequency (%)"))
# Pie chart for the composition at the Phylum level
t3 <- trans_abund$new(dataset = t2, taxrank = "Phylum", ntaxa = 8)
t3$data_abund$Sample %<>% factor(., levels = unique(.))
t3$plot_pie(facet_nrow = 3, color_values = c(RColorBrewer::brewer.pal(8, "Dark2"), "grey50"))

#### Test to see community composition differences among sites ####
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "KW", group = "Landcover", taxa_level = "all", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:20)
# Dunn's Kruskal-Wallis Multiple Comparisons when group number > 2
# Select the level of interest
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "KW_dunn", group = "Landcover", taxa_level = "Genus", filter_thres = 0.0001)
t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F)
# Wilcoxon Rank Sum and Signed Rank Tests for all paired groups
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "wilcox", group = "Landcover", taxa_level = "Genus", filter_thres = 0.001)
# filter something not needed to show
t1$res_diff %<>% subset(Significance %in% "***")
t1$plot_diff_abund()
# y_start and y_increase control the position of labels; for the details, please see the document of plot_alpha function in trans_alpha class
t1$plot_diff_abund(y_start = 0.05, y_increase = 0.1)
# anova
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "anova", group = "Landcover", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(coord_flip = F, plot_type = "barerrorbar", errorbar_addpoint = FALSE)
t1$res_diff
# LEfSe analysis (if used, cite Nicola Segata, Jacques Izard, Levi Walron, Dirk Gevers, Larisa Miropolsky, Wendy Garrett, Curtis Huttenhower.
# "Metagenomic Biomarker Discovery and Explanation" Genome Biology, 2011 Jun 24;12(6):R60)
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "lefse", group = "Landcover", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Genus")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Grassland", "TwoYearsOld", "TenYearsold", "TwentyfourYearsold", "Remnant"))
# show part of the table
t1$res_diff[1:5, c(1, 3, 4, 6)]
# Visualise better
t1$plot_diff_abund(use_number = 1:30)
t1$plot_diff_abund(fill = "Landcover", alpha = 0.5, add_sig = FALSE)
t1$plot_diff_abund(group_order = c("Grassland", "TwoYearsOld", "TenYearsold", "TwentyfourYearsold", "Remnant"))
t1$plot_diff_abund(coord_flip = FALSE)
t1$plot_diff_abund(plot_type = "errorbar")
t1$plot_diff_abund(plot_type = "barerrorbar", coord_flip = FALSE)
t1$plot_diff_abund(plot_type = "barerrorbar", errorbar_addpoint = FALSE, errorbar_color_black = TRUE, plot_SE = TRUE)
#t1$plot_diff_abund(plot_type = "ggviolin", coord_flip = FALSE)
# Random Forest
# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = meco_datasetMPJ, method = "rf", group = "Landcover", taxa_level = "Genus")
# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:20, group_order = c("Grassland", "TwoYearsOld", "TenYearsold", "TwentyfourYearsold", "Remnant"))
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("Grassland", "TwoYearsOld", "TenYearsold", "TwentyfourYearsold", "Remnant"), select_taxa = t1$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = FALSE, errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())
p <- g1 %>% aplot::insert_right(g2)
p
# For other methods see microeco tutorial Chapter 6.1
# Possibility to see the differences considering  lmm

#### General Network analysis ####
# Network analysis has been frequently used to study microbial co-occurrence patterns (Deng et al. 2012; Faust and Raes 2012;
# Coyte, Schluter, and Foster 2015). 
# In this part, we describe part of the implemented methods in the trans_network class.
# Using correlation method
# The parameter cor_method in trans_network is used to select correlation calculation method.
# default pearson or spearman correlation invoke R base cor.test, a little slow
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "spearman", filter_thres = 0.001)
# return t1$res_cor_p list, containing two tables: correlation coefficient table and p value table
# Spearman correlation based on WGCNA package is applied to show all the following operations:
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)
# The parameter COR_cut can be used to select the correlation threshold. Furthermore, COR_optimization = TRUE can be used to find the optimized coefficient threshold (potential transition point of network eigenvalues) 
# instead of the COR_cut based on the RMT theory (Deng et al. 2012).
# construct network; require igraph package
t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t1$res_network
# Do partition modules for the network.
# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")
# See chapter 6.2 in tutorial to save in Gephi format
# require rgexf package to be installed
# t1$save_network(filepath = "network.gexf")
# calculate network attributes
t1$cal_network_attr()
t1$res_network_attr
# The function get_node_table, get_edge_table and get_adjacency_matrix are designed 
# to get node properties table, edge properties table and adjacency matrix from network, respectively.
# get node properties
t1$get_node_table(node_roles = TRUE)
# return t1$res_node_table
t1$res_node_table
# get edge properties
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
# return t1$res_adjacency_matrix
# Plot the node classification in terms of the within-module connectivity and among-module connectivity
# add_label = TRUE can be used to directly add text label for points
t1$plot_taxa_roles(use_type = 1)
# plot node roles with phylum information
t1$plot_taxa_roles(use_type = 2)
# Now, we show the eigengene analysis of modules. 
# The eigengene of a module, i.e. the first principal component of PCA,
# represents the main variance of the abundance in the species of the module.
t1$cal_eigen()
# return t1$res_eigen
# Perform correlation heatmap to show the associations between eigengenes and environmental factors
# create trans_env object
Metadata_Rhizosphere<-Metadata[Metadata$ProjectFocus=="Rhizosphere",]
t2 <- trans_env$new(dataset = meco_datasetMPJ, add_data = Metadata_Rhizosphere[,21:31],complete_na = TRUE)
# calculate correlations
t2$cal_cor(add_abund_table = t1$res_eigen)
# plot the correlation heatmap
t2$plot_cor()
# Other functions are present to subset the network
# Network visualisation
# default parameter represents using igraph plot.igraph function
# use ggraph method; require ggraph package
t2 <- clone(t1)
t2$plot_network(method = "ggraph", node_color = "Phylum")
# use networkD3 package method for the dynamic network visualization in R
t1$plot_network(method = "networkD3", node_color = "module")
t1$plot_network(method = "networkD3", node_color = "Phylum")
#The trans_comm function can be used to convert the node classification to a new microtable object for other analysis.
# use_col is used to select a column of t1$res_node_table
tmp <- t1$trans_comm(use_col = "module", abundance = FALSE)
tmp
tmp$otu_table[tmp$otu_table > 0] <- 1
tmp$tidy_dataset()
tmp$cal_abund()
tmp2 <- trans_abund$new(tmp, taxrank = "Phylum", ntaxa = 10)
tmp2$data_abund$Sample %<>% factor(., levels = rownames(tmp$sample_table))
tmp2$plot_line(xtext_angle = 30, color_values = RColorBrewer::brewer.pal(12, "Paired")) + ylab("ASVs ratio (%)")
# The function cal_sum_links can sum the links (edge) number from one taxa to another or 
# within the same taxa. The function plot_sum_links is used to show the result from the function cal_sum_links. 
# This is very useful to fast see how many nodes are connected between different taxa or within one taxa. 
# In terms of ‘Phylum’ level in the tutorial, the function cal_sum_links() sum the linkages number from one Phylum to 
# another Phylum or the linkages in the same Phylum. 
# So the numbers along the outside of the circular plot represent how many edges or linkages are related with the Phylum. 
t1$cal_sum_links(taxa_level = "Phylum")
# interactive visualization
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
# Calculation options for the correlation network
# use Bray–Curtis index (1-dissimilarity)
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "bray", filter_thres = 0.001)
# Pearson correlation
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "pearson", filter_thres = 0.001)
# Pearson correlation
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "pearson", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.001)
# Pearson correlation
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "pearson", use_NetCoMi_pearson_spearman = TRUE, filter_thres = 0.001)
# Spearman correlation
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.001)
# Spearman correlation
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "spearman", use_NetCoMi_pearson_spearman = TRUE, filter_thres = 0.001)
# SparCC method
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "sparcc", use_sparcc_method = "SpiecEasi", filter_thres = 0.003)
# SparCC method
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "sparcc", use_sparcc_method = "NetCoMi", filter_thres = 0.001)
# CCLasso method
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "cclasso", filter_thres = 0.001)
# CCREPE method
t1 <- trans_network$new(dataset = meco_datasetMPJ, cor_method = "ccrepe", filter_thres = 0.001)

#### Network for each land cover type ####
# first create a list
soil_amp_network <- list()
# select samples of "Grassland" Landcover
# use clone to get a deep copy of soil_amp (R6 object)
tmp <- clone(meco_datasetMPJ)
# change sample_table directly
tmp$sample_table %<>% subset(Landcover == "Grassland")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
# put the network into the list
soil_amp_network$Grassland <- tmp
# select samples of "TwoYearsOld" group
tmp <- clone(meco_datasetMPJ)
tmp$sample_table %<>% subset(Landcover == "TwoYearsOld")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$TwoYearsOld <- tmp
# select samples of "TenYearsold" group
tmp <- clone(meco_datasetMPJ)
tmp$sample_table %<>% subset(Landcover == "TenYearsold")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$TenYearsold <- tmp
# select samples of "TwentyfourYearsold" group
tmp <- clone(meco_datasetMPJ)
tmp$sample_table %<>% subset(Landcover == "TwentyfourYearsold")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$TwentyfourYearsold <- tmp
# select samples of "Remnant" group
tmp <- clone(meco_datasetMPJ)
tmp$sample_table %<>% subset(Landcover == "Remnant")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$Remnant <- tmp
# Now we have the list soil_amp_network
# The function cal_module in meconetcomp package is designed to partition modules for all the networks in the list.
soil_amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
# we extracted all the res_network_attr tables in the networks and merged them into one final table by using cal_network_attr function in meconetcomp package.
tmp <- cal_network_attr(soil_amp_network)
# tmp is a data.frame object
# See the properties of the networks
tmp
# The get_node_table and get_edge_table functions of meconetcomp package can be used to directly extract node and edge properties for all the networks. The return table is stored in each network object.
soil_amp_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table
# The nodes in all the networks can be converted to a new microtable object by using 
# the node_comp function of meconetcomp package. Then, it is easy to analyse the nodes overlap with trans_venn class.
# obtain the node distributions by searching the res_node_table in the object
tmp <- node_comp(soil_amp_network, property = "name")
# obtain nodes intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
tmp1
g1 <- tmp1$plot_venn(fill_color = T)
g1
ggsave("soil_amp_node_overlap.pdf", g1, width = 7, height = 6)
# calculate bray distance to reflect the overall differences of networks
tmp$cal_betadiv(method = "bray")
tmp$beta_diversity$bray
# The pipeline of studying edges overlap is similar with the above operations 
# of nodes comparison. The edge_comp function of meconetcomp package is used to convert edges 
# distribution to a new microtable object.
# get the edge distributions across networks
tmp <- edge_comp(soil_amp_network)
# obtain edges intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
ggsave("soil_amp_edge_overlap.pdf", g1, width = 7, height = 6)
# calculate bray distance
tmp$cal_betadiv(method = "bray")
tmp$beta_diversity$bray
# Then we extracted the subset of edges according to the intersections of edges across networks, 
# which can be accomplished with the subset_network function in meconetcomp package.
# first obtain edges distribution and intersection
tmp <- edge_comp(soil_amp_network)
tmp1 <- trans_venn$new(tmp)
# convert intersection result to a microtable object
tmp2 <- tmp1$trans_comm()
tmp2$sample_table
# extract the intersection of all the five networks
# please use colnames(tmp2$otu_table) to find the required name
Intersec_all <- subset_network(soil_amp_network, venn = tmp2, name = "TwentyfourYearsold&Remnant")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
Intersec_all$save_network("Intersec_all.gexf")
# To know which taxa constitute the nodes in edges is important in understanding species 
# co-occurrence patterns and answering ecological questions. In this part, as an instance, 
# we used edge_tax_comp function of meconetcomp package to get the sums of node sources (at Phylum level) 
# in the positive edges. In other words, how many linked nodes of positive edges come from
# different phyla or the same phyla. Then, to make the results comparable, the ratio was calculated 
# with the positive edge number as denominator.
soil_amp_network_edgetax <- edge_tax_comp(soil_amp_network, taxrank = "Phylum", label = "+", rel = TRUE)
# filter the features with small number
soil_amp_network_edgetax <- soil_amp_network_edgetax[apply(soil_amp_network_edgetax, 1, mean) > 0.01, ]
# visualization
g1 <- pheatmap::pheatmap(soil_amp_network_edgetax, display_numbers = TRUE)
ggsave("soil_amp_edge_tax_comp.pdf", g1, width = 7, height = 7)
# Robustness of the network
tmp <- robustness$new(soil_amp_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
View(tmp$res_table)
View(tmp$res_summary)
tmp$plot(linewidth = 1)
# We can also extract the result for a specified metric and perform the visualization, such as the scatter plot.
tmp1 <- tmp$res_table %>% .[.$remove_strategy == "node_rand" & .$measure == "Eigen", ]
t1 <- trans_env$new(dataset = NULL, add_data = tmp1)
t1$dataset$sample_table <- t1$data_env
t1$plot_scatterfit(x = "remove_ratio", y = "value", type = "cor", group = "Network") + 
  xlab("Ratio of randomly removed nodes") + ylab("Network connectivity") + theme(axis.title = element_text(size = 15))
# another way
t1$plot_scatterfit(x = tmp1$remove_ratio, y = tmp1$value, type = "cor", group = tmp1$Network)
# the vulnerability of nodes can be calculated with the vulnerability function. 
# The vulnerability of one node is defined as the efficiency of network after removing this
# targeted node (Yuan et al. 2021). For the details, please swith to the help document.
vul_table <- vulnerability(soil_amp_network)
View(vul_table)
# The cohesion is a method for quantifying the connectivity of microbial communities (Herren and McMahon 2017). 
t1 <- cohesionclass$new(soil_amp_network)
View(t1$res_list$sample)
View(t1$res_list$feature)
t1$cal_diff(method = "anova")
t1$plot(measure = "r_pos")
#### Analyses with environmental parameters ####
# Analyses between environmental parameters
# First, it is better to clone the dataset
tmp_mt <- clone(meco_datasetMPJ)
# Metadata dataset for Rhizosphere
Metadata_Rhizosphere<-Metadata[Metadata$ProjectFocus=="Rhizosphere",]
Metadata_Rhizosphere[,21:30]
# Create the object and complete_na= True because 0 and NA 
t1 <- trans_env$new(dataset = tmp_mt, add_data = Metadata_Rhizosphere[,21:30],complete_na = TRUE)
# Difference cross groups with Wilcoxon
t1$cal_diff(group = "Landcover", method = "wilcox")
head(t1$res_diff)
# Do with Anova
t1$cal_diff(method = "anova", group = "Landcover")
t1$res_diff
# Observe the correlation among variables
t1$cal_autocor()
# Correlation for different landcover types
t1$cal_autocor(group = "Landcover")

# Do dbRDA ordination
# use bray-curtis distance for dbRDA
t1$cal_ordination(method = "dbRDA", use_measure = "bray", group="Landcover")
# show the orginal results
t1$trans_ordination()
t1$res_ordination_trans
t1$plot_ordination(plot_color = "Landcover")
# the main results of RDA are related with the projection and angles between arrows
# adjust the length of the arrows to show them better
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
P<-t1$plot_ordination(plot_color = "Landcover")
P + labs(title="Fungi db-RDA")
# The function cal_ordination_envfit can be used to get the contribution of each variables to the model.
t1$cal_ordination_anova()
t1$cal_ordination_envfit()
t1$res_ordination_envfit

# RDA at the Genus level
# use Genus
t1$cal_ordination(method = "RDA", taxa_level = "Genus")
# select 10 features and adjust the arrow length
t1$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
t1$plot_ordination(plot_color = "Landcover")
# See other visualisation styles in tutorial
# Do RDA checking for the significance of the differences
t1$cal_ordination(method = "RDA", taxa_level = "Genus", use_measure="bray")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
t1$trans_ordination(adjust_arrow_length = F)
g1 <- t1$plot_ordination(plot_color = "Landcover", plot_shape = "Landcover", type="taxa")
g1
ggplot2::ggsave("RDA.pdf", g1, width = 8, height = 6.5)
# use capture.output to save output
capture.output(t1$res_ordination_R2, file = "RDA_R2.txt")
capture.output(t1$res_ordination_envfit, file = "RDA_envfit.txt")
# save data.frame objects
write.table(t1$res_ordination_terms, "RDA_anova_termsig.txt", sep = "\t")
write.table(t1$res_ordination_axis, "RDA_anova_axissig.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_sites, "RDA_axis_sample.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_arrows, "RDA_axis_term.txt", sep = "\t")
write.table(t1$res_ordination_trans$df_arrows_spe, "RDA_axis_taxa.txt", sep = "\t")
# Check RDA which column to include/exclude

# Mantel test can be used to check whether there is 
# significant correlations between environmental variables and distance matrix.
# Only enviromental variables
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
head(t1$res_mantel)
# mantel test for different groups
t1$cal_mantel(by_group = "Landcover", use_measure = "bray")
# partial mantel test
t1$cal_mantel(partial_mantel = TRUE, by_group = "Landcover", method="spearman")
t1$res_mantel


# Heatmap to see correlation between taxa and environmental variables
t1 <- trans_env$new(dataset = meco_datasetMPJ, add_data = Metadata_Rhizosphere[,21:30],complete_na = TRUE)
# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.
t1$cal_cor(use_data = "Genus", p_adjust_method = "none", by_group="Landcover")
# return the results
t1$res_cor
# default ggplot2 method with clustering
t1$plot_cor()
# filter genera that do not have at least one ***
t1$plot_cor(filter_feature = c("", "*", "**"))
# use other_taxa to select taxa you need
# t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:40])
# t1$plot_cor()
# clustering heatmap; require pheatmap package
# Another color pallete
t1$plot_cor(pheatmap = TRUE, color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))

# To study the correlation between ALPHA diversity and environment
t1 <- trans_env$new(dataset = meco_datasetMPJ, add_data = Metadata_Rhizosphere[,21:30],complete_na = TRUE)
# use add_abund_table parameter to add the extra data table
t1$cal_cor(add_abund_table = t1$alpha_diversity, by_group="Landcover")
# try to use ggplot2 with clustering plot
t1$plot_cor(cluster_ggplot = "both")
# see the results
head(t1$res_cor)

#### ML Classification ####
# Classification with random forest
# Load and prepare the dataset
t1 <- trans_env$new(dataset = meco_datasetMPJ, add_data = Metadata_Rhizosphere[,21:30])
tmp <- t1$data_env %>% t %>% as.data.frame
# create the object
t2 <- trans_classifier$new(dataset = meco_datasetMPJ, x.predictors = "All", y.response = "Landcover")
# Pre-processing
t2$cal_preProcess(method = c("center", "scale", "nzv"))
# All samples are used in training if cal_split function is not performed (SO, IT MUST BE PERFORMED!)
# generate train and test set
t2$cal_split(prop.train = 3/4)
# Before training the model,  run the set_trainControl to invoke the trainControl function of caret package 
# to generate the parameters used for training. 
# Here are used the default parameters in trainControl function.
t2$set_trainControl()
# Train the model, method "rf" is by default
t2$cal_train(method = "rf")
# cal_predict function to predict the testing data set.
t2$cal_predict()
# plot the confusionMatrix to check out the performance
t2$plot_confusionMatrix()
# cal_ROC and plot_ROC to get the ROC (Receiver Operator Characteristic) curve.
t2$cal_ROC()
# select one group to plot ROC
t2$plot_ROC(plot_group = "Remnant")
t2$plot_ROC(plot_group = "Remnant", color_values = "black")
# default all groups
t2$plot_ROC(size = 0.5, alpha = 0.7)
# Feature selection to improve the accuracy and reduce overfit and complexity of the model
t2$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
# Perfom all the analysis with feature selection and show the results.
t2$cal_split(prop.train = 3/4)
t2$set_trainControl()
# Methods available are SVM (t3$cal_train(method = "svmRadial", tuneLength = 15))
t2$cal_train()
t2$cal_predict()
t2$plot_confusionMatrix()
t2$cal_ROC()
t2$plot_ROC(size = 0.5, alpha = 0.7)
# To plot the Precision-Recall curve (PR curve) make plot_type = “PR” in plot_ROC function.
t2$plot_ROC(plot_type = "PR", size = 0.5, alpha = 0.7)
# To show the ROC curve or PR curve of the training result make input = “train” in plot_ROC function.
t2$cal_ROC(input = "train")
t2$plot_ROC(plot_type = "ROC", size = 0.5, alpha = 0.7)
# To obtain the feature importance use cal_feature_imp function.
# default method in caret package without significance
t2$cal_feature_imp()
t2$plot_feature_imp(colour = "red", fill = "red", width = 0.6)
# generate significance with rfPermute package
t2$cal_feature_imp(rf_feature_sig = TRUE, num.rep = 1000)
# add_sig = TRUE: add significance label
t2$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6, add_sig = TRUE)
# show_sig_group = TRUE: show different colors in groups with different significance labels
t2$plot_feature_imp(show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE)
t2$plot_feature_imp(show_sig_group = TRUE, coord_flip = TRUE, width = 0.6, add_sig = TRUE)
# rf_sig_show = "MeanDecreaseGini": switch to MeanDecreaseGini
t2$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE)
# group_aggre = FALSE: donot aggregate features for each group
t2$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)

#### Mentel test, correlations between beta diversity and environmental variables ####
# Show the 4 most abundant phyla (select them from previous analyses)
# Prepare the dataset
d1 <- clone(meco_datasetMPJ)
d1$tax_table <- d1$tax_table[d1$tax_table$Phylum == "p__Ascomycota", ]
d1$tidy_dataset()
d1$cal_betadiv()
d2 <- clone(meco_datasetMPJ)
d2$tax_table <- d2$tax_table[d2$tax_table$Phylum == "p__Basidiomycota", ]
d2$tidy_dataset()
d2$cal_betadiv()
d3 <- clone(meco_datasetMPJ)
d3$tax_table <- d3$tax_table[d3$tax_table$Phylum == "p__Rozellomycota", ]
d3$tidy_dataset()
d3$cal_betadiv()
d4 <- clone(meco_datasetMPJ)
d4$tax_table <- d4$tax_table[d4$tax_table$Phylum == "p__Glomeromycota", ]
d4$tidy_dataset()
d4$cal_betadiv()
# first perform mantel test
t1 <- trans_env$new(dataset = d1, env_cols = 22:34)
t1$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
t2 <- trans_env$new(dataset = d2, env_cols = 22:34)
t2$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
t3 <- trans_env$new(dataset = d3, env_cols = 22:34)
t3$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
t4 <- trans_env$new(dataset = d4, env_cols = 22:34)
t4$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
# extract a part of the results 
x1 <- data.frame(spec = "Ascomycota", t1$res_mantel) %>% .[, c(1, 3, 6, 8)]
x2 <- data.frame(spec = "Basidiomycota", t2$res_mantel) %>% .[, c(1, 3, 6, 8)]
x3 <- data.frame(spec = "Rozellomycota", t3$res_mantel) %>% .[, c(1, 3, 6, 8)]
x4 <- data.frame(spec = "Glomeromycota", t4$res_mantel) %>% .[, c(1, 3, 6, 8)]
# rename columns
colnames(x1) <- colnames(x2) <-colnames(x3)<-colnames(x4)<- c("spec", "env", "r", "p.value")
# generate interval data
x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x3 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x4 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# cobine the tables
plot_table <- rbind(x1, x2,x3,x4)
# Visualise
set_scale()
g1 <- quickcor(t1$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "black", size = 6) +
  anno_link(aes(colour = pd, size = rd), data = plot_table) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

g1

#### Further stats ####
# see:https://chiliubio.github.io/microeco_tutorial/diversity-based-class.html
# The trans_env and trans_func classes are placed into the section 
# ‘Explainable class’, as environmental factors and microbial functions can
# be generally applied to explain microbial community structure and assembly.

#### Functional assignment ####
# Paragraph 7.2 of the tutorial (https://chiliubio.github.io/microeco_tutorial/explainable-class.html#trans_func-class)
# create object of trans_func
t2 <- trans_func$new(meco_datasetMPJ)
# mapping the taxonomy to the database. This can recognize prokaryotes or fungi automatically 
# if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html# fungi-data
# default database for prokaryotes is FAPROTAX database
# mapping
t2$cal_spe_func(prok_database = "FUNGuild")
# return t2$res_spe_func, 1 represent trait exists, 0 represent no or cannot confirmed.
t2$res_spe_func[1:5, 1:2]
# The percentages of the ASVs having the same trait 
# can reflect the functional redundancy of this function in the community.
# calculate the percentages for communities
# here consider the abundance
t2$cal_spe_func_perc(abundance_weighted = T)
# To see part of the result
t2$res_spe_func_perc[1:5, 1:2]
# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list<-list(t2$func_group_list,group=t2$sample_table$Landcover)
t2$func_group_list
# To see
t2$trans_spe_func_perc()
t2$res_spe_func_perc_trans
t2$plot_spe_func_perc()
#### Functions and general network ####
# construct a network to show the percentages of the OTUs for each trait in network modules
network <- trans_network$new(dataset = meco_datasetMPJ, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network$cal_network(p_thres = 0.01, COR_cut = 0.7)
network$cal_module()
soil_amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
# convert module info to microtable object
meco_module <- network$trans_comm(use_col = "module")
meco_module_func <- trans_func$new(meco_module)
meco_module_func$cal_spe_func(prok_database = "FUNGuild")
# Check if I want the abundance weighted or not
meco_module_func$cal_spe_func_perc(abundance_weighted = T)
meco_module_func$plot_spe_func_perc(order_x = paste0("M", 1:10)) + labs(title="Fungi functional groups in co-occurence networks")
# use show_prok_func to see the detailed information of traits
# Example
t2$show_prok_func("methanotrophy")

#### Functions and environment ####
# Correlation of the percentage data in res_spe_func_perc to environmental variables.
Metadata<-read.table("Fungi_Metadata.csv",h=T,sep=",",row.names = 1)
# Metadata dataset for Rhizosphere
Metadata_Rhizosphere<-Metadata[Metadata$ProjectFocus=="Rhizosphere",]
t3 <- trans_env$new(dataset = meco_datasetMPJ, add_data = Metadata_Rhizosphere[,21:30],complete_na = TRUE)
# Use t2 object from function analysis. Reload if needed.
t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman",by_group = "Landcover")
t3$res_cor
write.table(t3$res_cor, "CorrelationFunctionsEnv.csv", sep = "t")
# Visualise
t3$plot_cor(xtext_angle = 90)
# t3$plot_cor(pheatmap = TRUE, main = "Fungi functional groups and environmental parameters")

#### Differential test to see different functions across groups, land cover types ####
# First, it is better to clone the dataset
tmp_mt <- clone(meco_datasetMPJ)
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
str(t2$res_spe_func_perc)
tmp <- as.data.frame(t(t2$res_spe_func_perc), check.names = FALSE)
# assign the table back to taxa_abund list for further analysis
tmp_mt$taxa_abund$func <- tmp
tmp
# select the "func" in taxa_abund list in trans_diff
t2 <- trans_diff$new(dataset = tmp_mt, method = "anova", group = "Landcover", taxa_level = "func")
# See results of relative abundance
t2$res_abund
write.csv(t2$res_abund, file = "abundance_table.csv")
# See anova results
t2$res_diff
t2$plot_diff_abund(add_sig = TRUE) + ggplot2::ylab("Relative abundance (%)")+ggplot2::labs(title = "Fungi functional groups")
