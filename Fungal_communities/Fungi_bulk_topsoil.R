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
# install.packages("recipes")
# install_github("zdk123/SpiecEasi")
# install.packages('BAT')

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
library("recipes")
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
ps_BulkRhizo2YO <- subset_samples(ps_DE, ProjectFocus == "Rhizosphere")
# Bulk only for BulkRhizo1_2YO comparison in 2yo and 10yo plantation (only E2 and I1 transects)
ps_BulkRhizo1_2YO<-subset_samples(ps_DE, original_sample_id %in% 
                                    c("E2NA", "E2NB", "E2NC", "E2ND","I1SA", "I1SB", "I1SC", "I1SD"))
# Bulk only for BulkRhizo2YO comparison in 2yo plantation (only E2 transect)
ps_BulkRhizo2YO<-subset_samples(ps_DE, original_sample_id %in% c("E2NA", "E2NB", "E2NC", "E2ND"))
# Bulk only for BulkRhizo10YO comparison in 10yo plantation (only I1 transect)
ps_BulkRhizo10YO<-subset_samples(ps_DE, original_sample_id %in% c("I1SA", "I1SB", "I1SC", "I1SD"))

#### Processing sample-rarefying for each of the BulkRhizo1_2YO dataset #### 
# This step is done to ensure that the samples are comparable, having equal readings-depth
# (done for Main Project)
# orders samples from highest to lowest
sample_sums(ps_BulkRhizo1_2YO)[order(sample_sums(ps_BulkRhizo1_2YO))]
# Observe the number of reads, retain the minimum number similar to others (eg. if present 
# 18-490-23897-23786-[], retain 23897).
# We see that the sample 148 has really low reads (549). To discard. From Mahogany.
# Minimum to retain: 63,341
# To display the plots to recognize better were to cut:
asv_before_BulkRhizo1_2YO <- as(otu_table(ps_BulkRhizo1_2YO), "matrix")
out<-vegan::rarecurve(asv_before_BulkRhizo1_2YO, step=100,lwd=2, ylab="ASV Richness", xlab="Sequence Sample Size", main="INSERT_GENE_TARGET rRNA", label=F)
# Set seed to reproduce the data, since the rarefaction will sample
set.seed(100)
# Rarefy the samples without replacement. 
# Rarefaction is used to simulate an even number of reads per sample. 
ps_rarefiedBulkRhizo1_2YO <- rarefy_even_depth(ps_BulkRhizo2YO, sample.size =   63341 , replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# Allocation of objects after rarefying
asv_counts_rarefiedBulkRhizo1_2YO <- as.matrix(otu_table(ps_rarefiedBulkRhizo1_2YO))
After_rar_sequence_readsBulkRhizo1_2YO <- sum(asv_counts_rarefiedBulkRhizo1_2YO)
After_rar_sequence_readsBulkRhizo1_2YO
# Taxa
after_rar_num_asvsBulkBulkRhizo1_2YO<- ntaxa(ps_rarefiedBulkRhizo1_2YO)
after_rar_num_asvsBulkBulkRhizo1_2YO
# Check
summarize_phyloseq(ps_rarefiedBulkRhizo1_2YO)
ps_rarefiedBulkRhizo1_2YO
# Check the structure
str(sample_data(ps_rarefiedBulkRhizo1_2YO))

#### Processing sample-rarefying for each of the BulkRhizo2YO dataset #### 
# This step is done to ensure that the samples are comparable, having equal readings-depth
# (done for Main Project)
# orders samples from highest to lowest
sample_sums(ps_BulkRhizo2YO)[order(sample_sums(ps_BulkRhizo2YO))]
# Observe the number of reads, retain the minimum number similar to others (eg. if present 
# 18-490-23897-23786-[], retain 23897).
# We see that the sample 148 has really low reads (549). To discard. From Mahogany.
# Minimum to retain: 70842
# To display the plots to recognize better were to cut:
asv_before_BulkRhizo2YO <- as(otu_table(ps_BulkRhizo2YO), "matrix")
out<-vegan::rarecurve(asv_before_BulkRhizo2YO, step=100,lwd=2, ylab="ASV Richness", xlab="Sequence Sample Size", main="INSERT_GENE_TARGET rRNA", label=F)
# Set seed to reproduce the data, since the rarefaction will sample
set.seed(100)
# Rarefy the samples without replacement. 
# Rarefaction is used to simulate an even number of reads per sample. 
ps_rarefiedBulkRhizo2YO <- rarefy_even_depth(ps_BulkRhizo2YO, sample.size =   70842 , replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# Allocation of objects after rarefying
asv_counts_rarefiedBulkRhizo2YO <- as.matrix(otu_table(ps_rarefiedBulkRhizo2YO))
After_rar_sequence_readsBulkRhizo2YO <- sum(asv_counts_rarefiedBulkRhizo2YO)
After_rar_sequence_readsBulkRhizo2YO
# Taxa
after_rar_num_asvsBulkRhizo2YO<- ntaxa(ps_rarefiedBulkRhizo2YO)
after_rar_num_asvsBulkRhizo2YO
# Check
summarize_phyloseq(ps_rarefiedBulkRhizo2YO)
ps_rarefiedBulkRhizo2YO
# Check the structure
str(sample_data(ps_rarefiedBulkRhizo2YO))

#### Processing sample-rarefying for each of the BulkRhizo10YO dataset #### 
# This step is done to ensure that the samples are comparable, having equal readings-depth
# (done for Main Project)
# orders samples from highest to lowest
sample_sums(ps_BulkRhizo10YO)[order(sample_sums(ps_BulkRhizo10YO))]
# Observe the number of reads, retain the minimum number similar to others (eg. if present 
# 18-490-23897-23786-[], retain 23897).
# We see that the sample 148 has really low reads (549). To discard. From Mahogany.
# Minimum to retain: 63341
# To display the plots to recognize better were to cut:
asv_before_BulkRhizo10YO <- as(otu_table(ps_BulkRhizo10YO), "matrix")
out<-vegan::rarecurve(asv_before_BulkRhizo10YO, step=100,lwd=2, ylab="ASV Richness", xlab="Sequence Sample Size", main="INSERT_GENE_TARGET rRNA", label=F)
# Set seed to reproduce the data, since the rarefaction will sample
set.seed(100)
# Rarefy the samples without replacement. 
# Rarefaction is used to simulate an even number of reads per sample. 
ps_rarefiedBulkRhizo10YO <- rarefy_even_depth(ps_BulkRhizo10YO, sample.size =   63341 , replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# Allocation of objects after rarefying
asv_counts_rarefiedBulkRhizo10YO <- as.matrix(otu_table(ps_rarefiedBulkRhizo10YO))
After_rar_sequence_readsBulkRhizo10YO <- sum(asv_counts_rarefiedBulkRhizo10YO)
After_rar_sequence_readsBulkRhizo10YO
# Taxa
after_rar_num_asvsBulkRhizo10YO<- ntaxa(ps_rarefiedBulkRhizo10YO)
after_rar_num_asvsBulkRhizo10YO
# Check
summarize_phyloseq(ps_rarefiedBulkRhizo10YO)
ps_rarefiedBulkRhizo10YO
# Check the structure
str(sample_data(ps_rarefiedBulkRhizo10YO))


