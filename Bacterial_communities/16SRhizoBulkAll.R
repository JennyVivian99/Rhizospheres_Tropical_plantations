# Comparisons rhizospheres with bulk soil #
#### Load library ####
library(microeco)
library(ape)
library(vegan)
library(phyloseq)
library(ggtree)
library(dplyr)
library(dunn.test)
library(ggplot2)
library(indicspecies)

#### Set seed ####
set.seed(100)

#### Excel preparation ####

# Copy the sample table of the rhizosphere soil
# in each of the corresponding tables of the bulk soil analysis.
# Exclude from sample table and feature table the data not belonging to the 2 years-old plantation

#### Load data and clean ####

otu_dataBulk <- read.csv("feature_table.csv")
nrow(otu_dataBulk)
summary(otu_dataBulk)
otu_dataRhizosphere <- read.csv("feature_tableR.csv")
otuAll<-merge(otu_dataBulk,otu_dataRhizosphere, by="ID", all = TRUE)
nrow(otuAll)

# Change into Zero the NAs present
otuAll[is.na(otuAll)] <- 0
summary(otuAll)

# Make the first column the name of the rows
rownames(otuAll) <- otuAll$ID
summary(otuAll)

# Remove the ID column
otuAll<-otuAll[,-1]
summary(otuAll)

# Load taxa tables and merge them
taxa_dataBulk <- read.csv("tax_table.csv")
nrow(taxa_dataBulk)
taxa_dataRhizosphere <- read.csv("tax_tableR.csv")
taxaAll<-merge(taxa_dataBulk,taxa_dataRhizosphere,by = "ID", all = TRUE)
nrow(taxaAll)

# Make the first column the name of the rows
rownames(taxaAll) <- taxaAll$ID
summary(taxaAll)

# Remove the ID column
taxaAll<-taxaAll[,-1]
summary(taxaAll)

# Get rid of the duplicated columns
# Get the unique column prefixes (e.g., "Kingdom", "Phylum")
taxa_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Loop through each rank and merge the .x and .y columns
taxaAll_final <- taxaAll

for (rank in taxa_ranks) {
  # Create a new column with the combined name
  taxaAll_final[[rank]] <- coalesce(taxaAll_final[[paste0(rank, ".x")]], taxaAll_final[[paste0(rank, ".y")]])
}

# Remove the old .x and .y columns
taxaAll_final <- taxaAll_final %>%
  select(all_of(taxa_ranks))
# Check
print(taxaAll_final)

# Load sample table
sample_dataAll<- read.csv("sample_table.csv", row.names = 1)
summary(sample_dataAll)
# Adjust names for consistency with otu_data
rownames(sample_dataAll) <- gsub("-", ".", rownames(sample_dataAll))
sample_dataAll$sample_id <- gsub("-", ".",sample_dataAll$sample_id)

#### Create the microtable object ####
dataset <- microtable$new(
  otu_table = otuAll,
  tax_table = taxaAll_final,
  sample_table = sample_dataAll
)

#### Alpha diversity calculation ####
# To add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow, and the tree has to be provided
# If not specified, different indexes of diversity are calculated
dataset$cal_alphadiv(PD = F)
# return alpha_diversity in the object
class(dataset$alpha_diversity)
# save alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")
# See levels
dataset$sample_table$ProjectFocus <-as.factor(dataset$sample_table$ProjectFocus)
# Check
dataset$sample_table$ProjectFocus
# Creating a trans_alpha object can return two data.frame with the prefix ‘data_’: 
# data_alpha and data_stat. data_alpha is used for subsequent differential test 
# and visualization.
t1 <- trans_alpha$new(dataset = dataset, group = "ProjectFocus")
t1$data_alpha
# Save table alpha diversity with ASV counts, number of ASVs for each sample
write.csv(t1$data_alpha, "alpha_BulkRhizosphere.csv", row.names = TRUE)
# return t1$data_stat
head(t1$data_stat)
t1$data_stat
# Then, we test the differences among groups using Kruskal-Wallis Rank Sum Test (overall test when groups > 2),
# Wilcoxon Rank Sum Tests (for paired groups), 
# Dunn’s Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2) and Anova can also be used with multiple comparisons.
# Normality check
shapiro.test(t1$data_alpha$Value[t1$data_alpha$Measure=="Observed"])
t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)
t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
head(t1$res_diff)
t1$res_diff
# The result is stored in object$res_diff
# more options
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
# head(t1$res_diff)
t1$cal_diff(method = "wilcox")
head(t1$res_diff)
# t1$cal_diff(method = "t.test")
# head(t1$res_diff)
# Try KW_dunn
t1$cal_diff(method = "KW_dunn")
t1$plot_alpha(measure = "Observed", shape = "ProjectFocus")+labs(title="Bacteria ASVs richness 2yo")
# Aggregate to Phylum level
t1 <- trans_alpha$new(dataset = dataset$merge_taxa("Phylum"), group = "ProjectFocus")
t1$data_alpha
t1$data_stat
# Try KW_dunn
t1$cal_diff(method = "KW_dunn")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "ProjectFocus")+labs(title="Bacteria Phylum richness 2yo")

#### Beta diversity calculation without Phylogenetic tree ####
dataset$cal_betadiv(unifrac = F)
#  return beta_diversity list in the object
class(dataset$beta_diversity)
#  save beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")
# create trans_beta object
# For PCoA and NMDS, measure parameter must be provided.
# measure parameter should be either one of names(mt$beta_diversity) or a customized symmetric matrix
t1 <- trans_beta$new(dataset = dataset, group = "ProjectFocus", measure = "jaccard")
t1$dataset
# Visualisation of the beta diversity with NMDS
t1$cal_ordination(method = "NMDS")
specific_colors <- c("Field"= "#CDC9A5","Saccharum spp." = "darkorchid1", "Acacia mangiumY" = "deepskyblue1", 
                     "Coccus nucifera" = "chocolate4", "Imperata cylindrica"="red")
plot<-t1$plot_ordination(plot_color = "ProjectFocus", plot_type = "point", color_values = specific_colors)
plot + labs(title="Bacteria communities 2yo")+theme_bw()
# Visualisation of the beta diversity
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
t1$res_ordination
# plot the PCoA result with confidence ellipse
# Define colors
specific_colors <- c("Field2"= "#CDC9A5","Saccharum spp." = "darkorchid1", "Acacia mangiumY" = "deepskyblue1", 
                     "Coccus nucifera" = "chocolate4", "Imperata cylindrica"="red","Field10"= "#CDC9A5","Acacia mangiumO" = "deepskyblue1", "Mahogany (Swietenia macrophylla)" = "blue", 
                     "Narra (Pterocarpus indicus)" = "green3")
# Define custom shapes (e.g., 15=square, 16=circle, 17=triangle, 18=diamond, 4=cross)
specific_shapes <- c("Field2"= 15, "Saccharum spp." = 15, "Acacia mangiumY" = 15,
                     "Coccus nucifera" = 15, "Imperata cylindrica"= 15, "Field10"=16, "Acacia mangiumO" = 16, "Mahogany (Swietenia macrophylla)" = 16, 
                     "Narra (Pterocarpus indicus)" = 16)
plot<-t1$plot_ordination(plot_color = "ProjectFocus", plot_shape = "ProjectFocus",plot_type = "point", color_values = specific_colors, point_size = 4)
plot + labs(title="Bacteria communities 2 and 10yo")+theme_bw()+theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 14), 
                                                                   # Control the size of the tick mark labels on the Y-axis
                                                                   axis.text.y = element_text(size = 14))+scale_shape_manual(values = specific_shapes)

# Return exploring the previous general dataset and compute KW
t1 <- trans_beta$new(dataset = dataset, group = "ProjectFocus", measure = "jaccard")
t1$cal_manova(manova_all = F)
t1$res_manova
write.csv(t1$res_manova,"permanova16SRhizoBulk.csv",row.names = T)

# Test distances
t1$cal_group_distance()
t1$cal_group_distance_diff(method="KW_dunn")
t1$res_group_distance_diff

# Anosim analysis
t1$cal_anosim(group = "ProjectFocus")
t1$res_anosim
t1$cal_anosim(group = "ProjectFocus", paired = TRUE)
t1$res_anosim

# PERMDISP analysis
# PERMDISP(Anderson et al. 2011) is implemented to test multivariate homogeneity of groups dispersions (variances) based on the betadisper function of vegan package.
# for the whole comparison and for each paired groups
t1$cal_betadisper()
t1$res_betadisper

#### Composition-based class exploration with plots ####
# Composition-based class. These analyses are to visualise the taxonomic abundance, considering both the different ASVs,
# and the number of reads for each of them.
# The trans_abund class and trans_venn class are organised into the section ‘Composition-based class’, 
# since they are mainly used to show the composition information of communities.
# create trans_abund object
# Phyla level (ntaxa=): select top 5 abundant Phyla.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 5)
# t1 object now include the transformed abundance data t1$abund_data 
# and other elements for the following plotting
# Plots
# Adjusting landcover factors levels in the right order
custom_phylum_colors <- c(
  "Glomeromycota" = "purple",
  "Mucoromycota" = "steelblue",
  "Ascomycota" = "forestgreen",
  "Rozellomycota"="lightblue",
  "Basidiomycota" = "orange",
  "Others" = "grey70"
)
t1$plot_bar(color_values = custom_phylum_colors,others_color = "grey70", facet = "ProjectFocus", xtext_keep = FALSE, legend_text_italic = FALSE)+
  labs(title = "Relative abundance of fungi phyla Bulk-Rhizo 2-10yo")

# See important taxa for differences, Glomeromycota for fungi
# Initialize the differential abundance analysis object ---
# Using "KW_dunn" to perform Kruskal-Wallis test followed by Dunn's post-hoc.
# The analysis is performed at the Phylum level.
t_diff <- trans_diff$new(
  dataset = dataset, 
  method = "KW_dunn", 
  group = "ProjectFocus",
  taxa_level = "Phylum"
)
# See results
t_diff$res_diff
# Use the plot_diff_abund method, passing the ID to the 'taxa_names' argument.
specific_colors_named <- c(
  "Field2"= "#CDC9A5",
  "Saccharum spp." = "darkorchid1", 
  "Acacia mangiumY" = "deepskyblue1",
  "Coccus nucifera" = "chocolate4", 
  "Imperata cylindrica"="red",
  "Field10"= "#CDC9F9",
  "Acacia mangiumO" = "pink", 
  "Mahogany (Swietenia macrophylla)" = "darkblue", 
  "Narra (Pterocarpus indicus)" = "green3"
)
plot1 <- t_diff$plot_diff_abund(group_names = levels(t_diff$group_factor),
                                # The KW_dunn method automatically calculates group letters which are shown with add_sig = TRUE
                                add_sig = TRUE, use_number = 1:6,
                                color_values = specific_colors_named)
# Customize with ggplot2 functions
plot <- plot1 +
  labs(title = paste("Most abundant phylum")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(plot)
