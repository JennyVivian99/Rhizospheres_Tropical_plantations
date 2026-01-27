# Comparisons rhizospheres with bulk soil #
#### Load library ####
library(microeco)
library(ape)
library(vegan)
library(phyloseq)
library (ggtree)
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
t1$plot_alpha(measure = "Observed", shape = "ProjectFocus")+labs(title="Fungi ASVs richness 2yo")
# Aggregate to Phylum level
t1 <- trans_alpha$new(dataset = dataset$merge_taxa("Phylum"), group = "ProjectFocus")
t1$data_alpha
t1$data_stat
# Try KW_dunn
t1$cal_diff(method = "KW_dunn")
head(t1$res_diff)
t1$plot_alpha(measure = "Observed", shape = "ProjectFocus")+labs(title="Fungi Phylum richness 2yo")

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
plot + labs(title="Fungi communities 2yo")+theme_bw()
# Visualisation of the beta diversity
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
t1$res_ordination
# plot the PCoA result with confidence ellipse
specific_colors <- c("Field"= "#CDC9A5","Saccharum spp." = "darkorchid1", "Acacia mangiumY" = "deepskyblue1", 
                     "Coccus nucifera" = "chocolate4", "Imperata cylindrica"="red")
plot<-t1$plot_ordination(plot_color = "ProjectFocus", plot_type = "point", color_values = specific_colors, point_size = 4)
plot + labs(title="Fungi communities 2yo")+theme_bw()+theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 14), 
                                                            # Control the size of the tick mark labels on the Y-axis
                                                            axis.text.y = element_text(size = 14))

# Return exploring the previous general dataset and compute KW
t1 <- trans_beta$new(dataset = dataset, group = "ProjectFocus", measure = "jaccard")
t1$cal_manova(manova_all = F)
t1$res_manova

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
# Using specific colors for each taxa
custom_phylum_colors <- c(
  "Glomeromycota" = "purple",
  "Mucoromycota" = "steelblue",
  "Ascomycota" = "forestgreen",
  "Rozellomycota"="lightblue",
  "Basidiomycota" = "orange",
  "Others" = "grey70"
)
t1$plot_bar(color_values = custom_phylum_colors,others_color = "grey70", facet = "ProjectFocus", xtext_keep = FALSE, legend_text_italic = FALSE)+
  labs(title = "Relative abundance of fungi phyla Bulk-Rhizo 2yo")
# Without costomised colors
t1$plot_bar(others_color = "grey70", facet = "ProjectFocus", xtext_keep = FALSE, legend_text_italic = FALSE)+
  labs(title = "Relative abundance of fungi phyla Bulk-Rhizo 2yo")
# ASV level (ntaxa=): select top 5 abundant Phyla
t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 5)
# t1 object now include the transformed abundance data t1$abund_data 
# and other elements for the following plotting
# Plots
# Adjusting landcover factors levels in the right order
t1$plot_bar(others_color = "grey70", facet = "ProjectFocus", xtext_keep = FALSE, legend_text_italic = FALSE)+
  labs(title = "Relative abundance of fungi ASVs Bulk-Rhizo 2yo")

#### Test to see community composition differences among sites ####
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW_dunn", group = "ProjectFocus", taxa_level = "Phylum")
t1$res_abund
t1$plot_diff_abund()+
  labs(title = "Relative abundance of fungi phyla Bulk-Rhizo 2yo with KW")
t1$res_diff
# Focus on the phylum displaying difference
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW_dunn", group = "ProjectFocus", taxa_level = "Phylum")
t1$res_abund
t1$plot_diff_abund(use_number = 1:3)+
  labs(title = "Relative abundance of fungi phyla Bulk-Rhizo 2yo with KW")
# See the results
t1$res_diff
dataset
#### Indicator Species Analysis (IndVal) ####
# Examine the indicator genera for each substrate
# Aggregate Data (OTU -> Genus)
# Extract ASVs and taxa tables
otu_raw <- as.data.frame(dataset$otu_table)
tax_raw <- as.data.frame(dataset$tax_table)
# Merge them to link counts with Taxonomy
# and create a temporary column 'OTU_ID' to ensure safe merging
otu_raw$OTU_ID <- rownames(otu_raw)
tax_raw$OTU_ID <- rownames(tax_raw)
merged_data <- left_join(otu_raw, tax_raw, by = "OTU_ID")
# Filter and Aggregate
# IMPORTANT: filter out entries where Genus is unknown (e.g., "g__", NA, or empty).
# If not, "Unclassified" will likely become a statistically significant indicator 
# simply because it's a huge trash-bin category.
genus_aggregated <- merged_data %>%
  # Filter out unclassified genera - Adjust "g__" to match your specific taxonomy format
  filter(Genus != "g__" & Genus != "" & !is.na(Genus) & Genus != "Unclassified") %>% 
  # Group by the Genus name
  group_by(Genus) %>%
  # Sum the counts for all numeric columns (the samples)
  summarise(across(where(is.numeric), sum))
# Format for IndVal (Samples x Genera)
# Move Genus names back to row names and transpose
genus_matrix <- as.data.frame(genus_aggregated)
rownames(genus_matrix) <- genus_matrix$Genus
genus_matrix <- genus_matrix[, -1] # Remove the text 'Genus' column
# Transpose: IndVal expects Rows = Samples, Columns = Genera
comm_data_genus <- t(genus_matrix)
# Run IndVal on the Genus Data
group_variable <- "ProjectFocus" 
groups_factor <- factor(dataset$sample_table[[group_variable]])
# Run multipatt on the NEW aggregated matrix
indval_genus_results <- multipatt(
  x = comm_data_genus,      
  cluster = groups_factor, 
  func = "r.g",            
  control = how(nperm = 9999),
  duleg = TRUE # true if only single-group indicators are needed
)

summary(indval_genus_results)

# Extract and Clean Results
genus_sign <- as.data.frame(indval_genus_results$sign)
# Filter Significant results
alpha_threshold <- 0.05
genus_sign_clean <- genus_sign %>%
  filter(p.value <= alpha_threshold) %>%
  mutate(
    Genus = rownames(.),
    # Map group index to group name
    Group_Name = names(indval_genus_results$comb[1,])[index]
  ) %>%
  # Optional: Filter out combinations if you only want strict specialists
  filter(Group_Name %in% levels(groups_factor)) %>%
  select(Genus, Group_Name, stat, p.value) %>%
  arrange(Group_Name, desc(stat))
# See briefly the results
print("Significant Genus-Level Indicators:")
print(genus_sign_clean)
# Save the results
write.csv(genus_sign_clean, "significant_indicators_genus_level.csv")


# Old analysis
# Prepare Data for IndVal Analysis ---
# The abundance matrix (samples as rows, taxa as columns), and 
# a factor vector defining the groups (e.g., Rhizosphere vs. Bulk)
# Check the column name in your sample table that defines the group, here "ProjectFocus"
group_variable <- "ProjectFocus" 
# Extract the normalized OTU table (microeco stores samples in columns, we need rows)
# Use the raw OTU table here, as IndVal works best on counts or presence/absence,
# transposing it to fit the 'indicspecies' requirement (Samples x Taxa).
comm_data <- t(dataset$otu_table) %>% as.data.frame()
# Extract the grouping factor
groups_factor <- factor(dataset$sample_table[[group_variable]])
# Run the Indicator Value (IndVal) Analysis
# func = "r.g" calculates the Indicator Value (IndVal) based on abundance and frequency
# nperm = 9999 for a robust permutation test (adjust based on computing power)
indval_results <- multipatt(
  x = comm_data,      
  cluster = groups_factor, 
  func = "r.g",           
  control = how(nperm = 9999) 
)
# Process and Visualize Results
# Directly access the raw data frame of results from the multipatt object.
indval_df <- as.data.frame(indval_results$sign)
indval_df
# Filter the raw data frame for the significance threshold (P < 0.05)
alpha_threshold <- 0.05
indval_df_significant <- indval_df %>%
  filter(p.value <= alpha_threshold)
indval_df_significant
# Handle case where no significant indicators are found
if (nrow(indval_df_significant) == 0) {
  print(paste("!!! ATTENTION: No statistically significant indicator species found (P <= ", alpha_threshold, ").", sep=""))
  print("Skipping detailed merging and plotting steps.")
  
  # Print the console summary report (what you saw before) for reference:
  print(summary(indval_results, alpha = alpha_threshold))
  
  final_indicators_clean <- data.frame() 
  
} 
# Case where there are significant results
row.names(indval_df_significant)
indval_df_significant$ID <- rownames(indval_df_significant)
indval_df_significant$ID
group_names_map <- names(indval_results$comb[1,])
group_names_map
# Explicitly create ID and Group columns outside of the main pipeline
# Explicitly create the 'ID' column from the Taxa row names
indval_df_significant$ID <- rownames(indval_df_significant)
# Explicitly map the numerical 'index' column to the actual string group name
group_names_map <- names(indval_results$comb[1,])
indval_df_significant$group <- group_names_map[indval_df_significant$index] 
indval_df_significant$group
# Prepare Taxonomy Table: Ensure taxonomy row names are explicitly an 'ID' column
taxa_info_with_id <- dataset$tax_table %>% 
  as.data.frame() %>%
  mutate(ID = rownames(.))
# Merge and Select Clean Columns
final_indicators_clean <- indval_df_significant %>%
  # Use left_join (dplyr) for a clean merge with the taxonomy data
  left_join(taxa_info_with_id, by = "ID") %>%
  dplyr::select(
    Taxon_ID = ID, 
    Indicator_Value = stat,
    Indicator_Group = group, 
    P_Value = p.value,
    Kingdom, Phylum, Class, Order, Family, Genus
  ) %>%
  # Filter for species that are indicators for only single groups (not combinations of groups).
  filter(Indicator_Group %in% levels(groups_factor)) %>%
  arrange(desc(Indicator_Value))
print("Significant Indicator Species Table (P <= 0.05)")
print(final_indicators_clean)
# Save the table
write.csv(final_indicators_clean, "final_indicators_clean.csv", row.names = TRUE)

# Plot the abundance of the top indicators
# Find the Taxon IDs of the top 5 indicators overall
top_taxa_ids <- final_indicators_clean$Taxon_ID[1:min(5, nrow(final_indicators_clean))]
top_taxa_ids  
# Use the trans_diff class from microeco for convenient visualization
t1 <- trans_diff$new(dataset = dataset, method = "KW_dunn", group = group_variable)
# Get the detected group names for correct plotting order
plot_group_order <- unique(final_indicators_clean$Indicator_Group)
# Plot the abundance of the top indicator taxa across the groups
p_indval_abund <- t1$plot_diff_abund(taxa_names = top_taxa_ids, group_names = levels(groups_factor),
     group_order = plot_group_order)
# Print the plot 
print(p_indval_abund)
# Save the plot
# ggsave(p_indval_abund, filename = "Top_Indicator_Abundance_Plot.pdf", width = 8, height = 6)

#### Functional assignment ####
# Paragraph 7.2 of the tutorial (https://chiliubio.github.io/microeco_tutorial/explainable-class.html#trans_func-class)
# create object of trans_func
t2 <- trans_func$new(dataset)
colnames(dataset$tax_table)
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
t2$res_spe_func_perc
t2$res_spe_func_raw_funguild
# To see how many were addressed  in terms of functions
NumberAssigned<-t2$res_spe_func_raw_funguild
NumberAssigned[NumberAssigned== ""] <- NA
NumberAssigned
unassigned_rows <- NumberAssigned[rowSums(is.na(NumberAssigned[, -1])) == (ncol(NumberAssigned) - 1), ]
number_of_unassigned_rows <- nrow(unassigned_rows)
# Print the result
print(paste("The number of rows with empty columns besides the first one is:", number_of_unassigned_rows))
print(paste("on a total number of:", nrow(NumberAssigned)))
# Save the result
write.table(t2$res_spe_func_perc, "FunctionFungiBulkRhizo2YO.csv", row.names = TRUE)
# We used the saved table to retrieve the mean and SD relative abundance of functional groups for each species's rhizosphere assessed
# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list<-list(t2$func_group_list,group=t2$sample_table$Rhizosphere)
t2$func_group_list
# To see
t2$trans_spe_func_perc()
t2$res_spe_func_perc_trans
t2$plot_spe_func_perc()

#### Multivariate plot with functional data ####
# t9 is created from the FunctionFungiBulkRhizo2YO.csv dataset
# t10 is created from the FunctionFungiBulkRhizo2YO.csv functions visible
t9<-read.csv("FunctionFungi_MultiPlotBulkRhizo2yo.csv", row.names=1)
t10<-read.csv("Functions_names.csv", row.names=1)
m2 <- microtable$new(sample_table = as.data.frame(dataset$sample_table), otu_table = as.data.frame(t(t9)), tax_table = t10)
# Correlation for different Rhizospheres and field values
m2$cal_betadiv(method = "jaccard")
t12 <- trans_beta$new(dataset = m2, group = "ProjectFocus", measure = "jaccard")
# PCoA, PCA and NMDS are available
t12$cal_ordination(method = "PCA")
# ordination result list
class(t12$res_ordination)
write.table(t12$res_ordination$scores, "PCAScoresFunctBulkRhizo2YO.csv", row.names = TRUE)
# plot the PCoA results
t12$plot_ordination(plot_color = "ProjectFocus",plot_type = c("point", "ellipse"))+theme_bw()+geom_text(aes(label = sample_id),vjust = -0.6)
t12$plot_ordination(plot_color = "ProjectFocus",plot_type = c("point"))+theme_bw()+ggtitle("Functional Fungi Bulk and Rhizosphere 2yo")
t12$res_ordination

#### Differential test to see different functions across groups, Rhizospheres and field ####
# First, it is better to clone the dataset
tmp_mt <- clone(dataset)
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
str(t2$res_spe_func_perc)
tmp <- as.data.frame(t(t2$res_spe_func_perc), check.names = FALSE)
# assign the table back to taxa_abund list for further analysis
tmp_mt$taxa_abund$func <- tmp
tmp
# select the "func" in taxa_abund list in trans_diff
t2 <- trans_diff$new(dataset = tmp_mt, method = "KW_dunn", group = "ProjectFocus", p_adjust_method="none", taxa_level = "func")
# See results of relative abundance
t2$res_abund
write.csv(t2$res_abund, file = "FungiFunctAbundanceBulkRhizo2yo_table.csv")
# See KW_dunn results
t2$res_diff
write.csv(t2$res_diff, file = "KW_Functions_abundance_tableBulkRhizo2yo.csv")
t2$plot_diff_abund(add_sig = TRUE, use_number = 1:9) + ggplot2::ylab("Relative abundance (%)")+ggplot2::labs(title = "Fungi functional groups Bulk-Rhizo 2yo")