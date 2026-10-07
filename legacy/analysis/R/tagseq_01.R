# import DESeq2 environment and graphical parameters ####
setwd("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/")
load("../06_DESeq2/DESeq2_environment.RData")

source("../../../01_script/R/graphical_parameters.R")

# create output directory
dir.create("figures", showWarnings = F, recursive = T)

#
# import libraries & functions ####

library(tidyverse)
library(ggbreak)
library(patchwork)
library(pheatmap)
library(svglite)
library(ComplexHeatmap)
library(agricolae)
library(ggbeeswarm)
library(ggpubr)

# import functions
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/superplot_w_boxplot.R")
source("../../../01_script/post_processing/R_functions/DEG_heatmap.R")
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/POSTHOC_Dunn_with_facets.R")
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/graphical_parameters.R")

#
# import clusters based on cdca7 upregulation ####
# 3 clusters with both cdca7 alleles

cdca7_log2FC_clusters_files <- list.files("/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/data/ddm1_upTEs_with_AS_cdca7ab_log2FC_x3_clusters"
                                          , full.names = TRUE, pattern="ddm1_upTEs_with_AS_cdca7ab_log2FC_cluster_._deeptools.bed")

# import data
cdca7_log2FC_clusters <- lapply(X = cdca7_log2FC_clusters_files, FUN = read.delim, head=F, sep="\t", quote="", comment.char="") # read them
names(cdca7_log2FC_clusters) <- basename(cdca7_log2FC_clusters_files)

# Merge all dataframes with cluster information
merged_cdca7_log2FC_clusters <- do.call(rbind, lapply(names(cdca7_log2FC_clusters), function(name) {
  df <- cdca7_log2FC_clusters[[name]]
  # Extract cluster number from filename
  cluster_num <- gsub(".*cluster_([0-9]+)_.*", "\\1", name)
  df$cluster <- paste0("cluster_", cluster_num)
  return(df)
}))

#
# import density of 5mC contexts at TEs ####

# import data
density_TE_files <- list.files("/groups/berger/user/pierre.bourguet/genomics/Araport11/nucleotide_content_TE_TAIR10", full.names = TRUE, pattern="*tsv")
density_TE <- lapply(X = density_TE_files, FUN = read.delim, head=T, sep="\t", quote="", comment.char="") # read them
names(density_TE) <- basename(density_TE_files)

merged_density_TE <- as_tibble(density_TE[[1]][,c(4,6,7,9,16,17)])
names(merged_density_TE) <- c("Geneid", "family", "superfamily", "GC", "length", gsub(".tsv", "_density", gsub("TAIR10_TE_", "", names(density_TE)[1])))

# Add only the 17th column from other dataframes
for (i in 2:length(density_TE)) {
  # Extract the unique column (17th column)
  unique_col <- density_TE[[i]][, 17]
  
  # Name the column based on the filename
  col_name <- gsub(".tsv", "_density", gsub("TAIR10_TE_", "", names(density_TE)[i]))
  
  # Add the column to the base dataframe
  merged_density_TE[[col_name]] <- unique_col
}

# View the merged dataframe
head(merged_density_TE)

# Create sum of all three densities  
merged_density_TE$CG_CHG_CHH_density <- merged_density_TE$CG_density + merged_density_TE$CHG_density + merged_density_TE$CHH_density  

# Create all possible combinations of two densities  
merged_density_TE$CG_CHG_density <- merged_density_TE$CG_density + merged_density_TE$CHG_density  
merged_density_TE$CG_CHH_density <- merged_density_TE$CG_density + merged_density_TE$CHH_density  
merged_density_TE$CHG_CHH_density <- merged_density_TE$CHG_density + merged_density_TE$CHH_density

# # merge with mC data
# dataframe_names <- c(
#   "CG_TE", "CHG_TE", "CHH_TE",
#   "CG_TE_avg", "CHG_TE_avg", "CHH_TE_avg",
#   "CG_TE_avg_norm", "CHG_TE_avg_norm", "CHH_TE_avg_norm",
#   "CG_TE_avg_ratio", "CHG_TE_avg_ratio", "CHH_TE_avg_ratio",
#   "mC_TE_avg", "mC_TE_avg_norm", "mC_TE_avg_ratio"
# )
# 
# # Loop through the list of dataframe names
# for (df_name in dataframe_names) {
#   # Dynamically reassign the updated dataframe back to the same name
#   assign(df_name, get(df_name) %>%
#            left_join(merged_density_TE, by = "Geneid"))
# }

#

# extract TEs upregulated in mutants ####

########################### all mutants

# find TEs upregulated in each mutant
up_TEs <- lapply(FUN = up_TEs_intersect_filter, X = all_res, TE_df="TEs", feature_df="features")

# extract nb of upregulated TEs for each genotype and convert to a dataframe
nb_up_TEs <- lapply(FUN = length, X = up_TEs)
nb_up_TEs_df <- data.frame(
  condition = gsub("res_", "", names(nb_up_TEs)),
  Count = unlist(nb_up_TEs)
)

########################### cdca7

# find TEs up in cdca7
up_TEs_cdca7 <- lapply(FUN = up_TEs_intersect_filter, X = list(all_res$res_ab_1, all_res$res_ab_2), TE_df="TEs", feature_df="features")
# extract unique TE ids
up_TEs_cdca7_ids <- intersect(up_TEs_cdca7[[1]], up_TEs_cdca7[[2]])
up_TEs_each_cdca7_ids <- union(up_TEs_cdca7[[1]], up_TEs_cdca7[[2]])

# find proportion of TEs up in cdca7 that are also up in ddm1
ddm1_G2_up_TEs <- up_TEs_intersect_filter(all_res$res_ddm1_2_G2, TE_df="TEs", feature_df="features")
sum(up_TEs_cdca7[[1]] %in% ddm1_G2_up_TEs) / length(up_TEs_cdca7[[1]])
sum(up_TEs_cdca7[[2]] %in% ddm1_G2_up_TEs) / length(up_TEs_cdca7[[2]])

# filter out ddm1 samples, a-long samples & select TEs upregulated in cdca7
vsd_cdca7_up_TEs <- vsd %>%
  dplyr::filter(Geneid %in% up_TEs_cdca7_ids) %>%
  tidyr::pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>%
  dplyr::mutate(condition = gsub("_R[0-9]", "", sample)) %>%
  dplyr::filter(!str_detect(condition, "ddm1_[^2]|F2")) %>%
  dplyr::filter(!str_detect(condition, "ddm1_2_G5")) %>%
  dplyr::filter(!str_detect(condition, "long")) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))

vsd_avg_cdca7_up_TEs <- vsd_avg %>%
  tidyr::pivot_longer(-Geneid, names_to = "condition", values_to = "value") %>%
  dplyr::filter(!str_detect(condition, "ddm1_[^2]|F2")) %>%
  dplyr::filter(!str_detect(condition, "ddm1_2_G5")) %>%
  dplyr::filter(!str_detect(condition, "long")) %>%
  dplyr::filter(Geneid %in% up_TEs_cdca7_ids) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))

########################### ddm1

# find TEs up in any of the ddm1 control
up_TEs_ddm1_control <- lapply(FUN = up_TEs_intersect_filter, X = list(all_res$res_ddm1_2_G2, all_res$res_ddm1_2_G5, all_res$res_F2_ddm1), TE_df="TEs", feature_df="features")
# extract unique TE ids
up_TEs_ddm1_control_ids <- unique(unlist(up_TEs_ddm1_control))

# Prepare the matrices and other variables based on the normalization method
vsd_ddm1_up_TEs <- vsd %>%
  dplyr::filter(Geneid %in% up_TEs_ddm1_control_ids) %>%
  tidyr::pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>%
  dplyr::mutate(condition = gsub("_R[0-9]", "", sample)) %>%
  dplyr::filter(str_detect(condition, "ddm1|Col_0|ab")) %>%
  dplyr::filter(!str_detect(condition, "long|NBD"))

#
# relative effect size ####

# find how many TEs are up compared with ddm1
length(up_TEs_cdca7[[1]]) / length(ddm1_G2_up_TEs) * 100
length(up_TEs_cdca7[[2]]) / length(ddm1_G2_up_TEs) * 100

# 
log2FC_df %>%
  select(Geneid, ddm1_2_G2, ab_1, ab_2) %>%
  filter(Geneid %in% up_TEs_cdca7_ids) %>%
  pivot_longer(cols = -Geneid, names_to = "condition", values_to = "log2FC") %>%
  group_by(condition) %>%
  summarize(average_log2FC = mean(log2FC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_log2FC = average_log2FC / average_log2FC[condition == "ddm1_2_G2"])

log2FC_df %>%
  select(Geneid, ddm1_2_G2, ab_1, ab_2) %>%
  filter(Geneid %in% up_TEs_each_cdca7_ids) %>%
  pivot_longer(cols = -Geneid, names_to = "condition", values_to = "log2FC") %>%
  group_by(condition) %>%
  summarize(average_log2FC = mean(log2FC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_log2FC = average_log2FC / average_log2FC[condition == "ddm1_2_G2"])

#
# barplot: number of up TEs in cdca7 mutants ####

# look at all the data
ggplot(nb_up_TEs_df, aes(y = condition, x = Count)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  theme_minimal()

# filter out ddm1 samples, a-long samples, and add a WT level
nb_up_TEs_df_no_ddm1 <- nb_up_TEs_df %>%
  filter(!str_detect(condition, "ddm1_[^2]|F2")) %>%
  filter(!str_detect(condition, "ddm1_2_G5")) %>%
  filter(!str_detect(condition, "long")) %>%
  bind_rows(tibble(condition = "WT", Count = NA))

# reorder levels
fig1_nb_up_TEs_df$condition <- factor(fig1_nb_up_TEs_df$condition, levels = rev(c("WT", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")))

fig1_up_TEs_hist_WT <- ggplot(fig1_nb_up_TEs_df, aes(y = condition, x = Count, fill = condition)) +
  geom_col() +
  labs(x = "Count", y = "Genotype") +
  scale_y_discrete(position = "left") +
  scale_fill_manual(values = rev(col_muted_2_replicates[c(19,1:4,7,8,5)])) +
  scale_x_continuous(breaks = c(0,100,200,300), limits=c(0,1110), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(330, 990), ticklabels=c(1000,1100), scales=0.2, space = 0.1, expand = F) +
  theme_horizontal_nature +
  theme(
    #panel.border = element_rect(fill=NA, colour = "black", linewidth = pt_0.5_to_mm),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(margin = margin(t = -5)), # reduce space between axis titles and axis labels
    axis.title.y = element_text(margin = margin(r = -5), angle = 90),
    plot.margin = margin(t = -5, r = -5)
    , legend.position = "none"
  ) ; fig1_up_TEs_hist_WT

# export in svg, works fine but fonts are not exported properly (open the svg file in a text editor to see).
svglite::svglite(filename = "figures/up_TEs_cdca7_histogram_with_WT.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
fig1_up_TEs_hist_WT
dev.off()

#
# heatmap of cdca7 upregulated TEs (fig2) ####

# remove Geneid, reorder columns and convert to matrix
m <- vsd_avg_cdca7_up_TEs %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  select(all_of(c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2"))) %>%
  as.matrix()

matrix_range <- seq(from = range(m)[1], to = range(m)[2], by = 1)
color_scale <- colorRamp2(matrix_range, scico(n = length(matrix_range), direction = -1, palette = "lajolla"))

# Plot the heatmap
Heatmap(m,
        name = "log2\n(VST+1)",
        col = color_scale,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width = 15, height = 15,
        column_title = paste0("cdca7α/β upTEs\nn=", length(up_TEs_cdca7_ids))
)

# remove clustering & sort the heatmap by row sums (TEs with highest expression across genotypes at the top)
m <- m[order(rowSums(m), decreasing = T),] 
fig1_up_TEs_heatmap <- Heatmap(t(m),
                               name = "log2\n(VST+1)",
                               col = color_scale,
                               cluster_rows = F,
                               cluster_columns = F,
                               width=15, height=15,
                               column_title=NULL,
                               row_names_side = "left",
                               row_names_gp = grid::gpar(fontsize = 6),
                               use_raster = F,
                               border= T,
                               heatmap_legend_param = list(
                                 title = "VST", at = seq(0,12,6), 
                                 #labels = c("0", "5", "10"),
                                 legend_height = unit(3, "cm"),
                                 legend_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 6),
                                 title_gp = gpar(fontsize = 6)
                               )
) ; fig1_up_TEs_heatmap

svglite::svglite(filename = "figures/heatmap_cdca7_up_TEs_VST.svg", width = 60*mm_to_inches, height = 40*mm_to_inches)
draw(fig1_up_TEs_heatmap)
dev.off()

# Tukey's HSD
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/ANOVA_Tukey_HSD.R")
ANOVA_Tukey_HSD(data = vsd_avg_cdca7_up_TEs, response = "value", factor = "condition")

# calculate average expression per condition, and difference with Col_0
vsd_avg_cdca7_up_TEs %>%
  group_by(condition) %>%
  summarize(mean = mean(value)) %>%
  mutate(diff = mean - mean[condition == "Col_0"])
#
# superplot: TEs up in cdca7 mutants (VST counts) ####

################## filter normalized counts and format for superplot

# Define the specific order vector
order_vector <- c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2")

# Define the color palette
my_colors <- rev(col_muted_2_replicates[c(19,1:4,7,8,5)])

# run the function and customize the plot
superplot_upTEs_cdca7 <- superplot_w_boxplot(vsd_cdca7_up_TEs, rev(order_vector), my_colors) + 
  coord_flip(ylim=c(2, 13)) +
  xlab("Genotype") + ylab("Transcript levels log2(ESF+1)") + labs(title="upregulated TEs")
superplot_upTEs_cdca7

# write svg output
svglite::svglite(filename = "figures/superplot_TEs_up_cdca7.svg", width = 80*mm_to_inches, height = 50*mm_to_inches)
superplot_upTEs_cdca7 & theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
means <- vsd_cdca7_up_TEs %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = mean(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = means)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
# barplot: number of up TEs in ddm1 cdca7 mutants ####

# filter samples
nb_up_TEs_df_w_ddm1 <- nb_up_TEs_df %>%
  #filter(str_detect(condition, "ddm1|ab")) %>% # use this one to include cdca7ab controls
  filter(str_detect(condition, "ddm1")) %>%
  filter(!str_detect(condition, "long|NBD")) %>%
 # filter(!condition %in% "F2_ddm1") %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
#nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G2", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_2_G5", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))
nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ddm1_2_G2", "ddm1_a_1", "F2_ddm1", "F2_ddm1_a_2", "ddm1_2_G5", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_ddm1 <- ggplot(nb_up_TEs_df_w_ddm1, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,3,5,3,5,3,5,5,5,5)]))) +
  theme(legend.position = "none") ; barplot_nb_up_TEs_ddm1

#

# barplot: number of up TEs in ddm1 cdca7 mutants, using DESeq2 to contrast ddm1 cdca7 mutants with their respective ddm1 control ####
library(DESeq2)

# ddm1_a_1 versus ddm1 G2
ddm1_a_1_vs_ddm1_G2 <- results(dds, contrast=c("condition", "ddm1_a_1", "ddm1_2_G2"))

# F2 ddm1_a_2 versus F2 ddm1
F2_ddm1_a_2_vs_F2_ddm1 <- results(dds, contrast=c("condition", "F2_ddm1_a_2", "F2_ddm1"))

# ddm1_b and ddm1_ab versus ddm1 G5
ddm1_b_1_vs_ddm1_G5 <- results(dds, contrast=c("condition", "ddm1_b_1", "ddm1_2_G5"))
ddm1_b_2_vs_ddm1_G5 <- results(dds, contrast=c("condition", "ddm1_b_2", "ddm1_2_G5"))
ddm1_ab_1_vs_ddm1_G5 <- results(dds, contrast=c("condition", "ddm1_ab_1", "ddm1_2_G5"))
ddm1_ab_2_vs_ddm1_G5 <- results(dds, contrast=c("condition", "ddm1_ab_2", "ddm1_2_G5"))

# create a new df with the number of upregulated TEs in each contrast
up_TEs_ddm1_cdca7 <- data.frame(
  condition = c("ddm1_a_1_vs_ddm1_G2", "F2_ddm1_a_2_vs_F2_ddm1", "ddm1_b_1_vs_ddm1_G5", "ddm1_b_2_vs_ddm1_G5", "ddm1_ab_1_vs_ddm1_G5", "ddm1_ab_2_vs_ddm1_G5"),
  Count = unlist(lapply(FUN=length, X=lapply(FUN = up_TEs_intersect_filter, X = list(ddm1_a_1_vs_ddm1_G2, F2_ddm1_a_2_vs_F2_ddm1, ddm1_b_1_vs_ddm1_G5, ddm1_b_2_vs_ddm1_G5, ddm1_ab_1_vs_ddm1_G5, ddm1_ab_2_vs_ddm1_G5), TE_df="TEs", feature_df="features")))
)

# filter samples
nb_up_TEs_df_w_ddm1 <- nb_up_TEs_df %>%
  filter(str_detect(condition, "ddm1_2_G|F2_ddm1$")) %>% # use this one include cdca7ab controls
  #bind_rows(tibble(condition = "Col_0", Count = 0)) %>%
  bind_rows(up_TEs_ddm1_cdca7)

# reorder levels
nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_2", "F2_ddm1", "ddm1_2_G2", "ddm1_2_G5", "F2_ddm1_a_2_vs_F2_ddm1", "ddm1_a_1_vs_ddm1_G2", "ddm1_b_1_vs_ddm1_G5", "ddm1_b_2_vs_ddm1_G5", "ddm1_ab_1_vs_ddm1_G5", "ddm1_ab_2_vs_ddm1_G5")))

barplot_nb_up_TEs_ddm1 <- ggplot(nb_up_TEs_df_w_ddm1, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,3,3,3)], rep(col_muted[5], 6)))) +
  theme_horizontal_nature +
  theme(
    panel.border = element_rect(fill=NA, colour = "black", linewidth = pt_0.5_to_mm),
    legend.position = "none"
  ) ; barplot_nb_up_TEs_ddm1

svglite::svglite(filename = "figures/barplot_nb_TEs_up_ddm1_cdca7.svg", width = 70*mm_to_inches, height = 45*mm_to_inches)
barplot_nb_up_TEs_ddm1
dev.off()

#
# superplot: TEs up in ddm1 cdca7 mutants (VST counts) ####

# Define the specific order vector
#order_vector <- c("Col_0", "ab_1", "ab_2" , "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")
order_vector <- c("Col_0", "ab_1", "ab_2", "F2_ddm1", "F2_ddm1_a_2", "ddm1_2_G2", "ddm1_a_1", "ddm1_2_G5", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")
#order_vector <- c("Col_0", "ab_2", "F2_ddm1", "F2_ddm1_a_2", "ddm1_2_G2", "ddm1_a_1", "ddm1_2_G5", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")

# Define the color palette
#my_colors <- c("grey55", col_muted[c(4,4,3)], rep(col_muted[5], 6)) # to include cdca7ab controls
my_colors <- c("grey55", col_muted[c(4,4,3,5,3,5,3,5,5,5,5)]) # without

# run the function and customize the plot
superplot_upTEs_ddm1 <- superplot_w_boxplot(vsd_ddm1_up_TEs, rev(order_vector), rev(my_colors)) + 
  coord_flip(ylim=c(2.5, 10.2)) +
  xlab("Genotype") + ylab("variance-stabilized counts") + labs(title="upregulated TEs")
superplot_upTEs_ddm1

# write svg output
svglite::svglite(filename = "figures/superplot_TEs_up_ddm1_cdca7_VST.svg", width = 50*mm_to_inches, height = 60*mm_to_inches)
set.seed(55)
superplot_upTEs_ddm1 + theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
means <- vsd_ddm1_up_TEs %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = mean(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = means)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#

######################### normalize expression values in the mutant by their respective WT or ddm1 controls
library(dplyr)

# Define a mapping of samples to their respective controls, including those normalized by Col_0
control_mapping <- tribble(
  ~condition,     ~control,
  "F2_ddm1_a_2",  "F2_ddm1",
  "ddm1_a_1",     "ddm1_2_G2",
  "ddm1_ab_1",    "ddm1_2_G5",
  "ddm1_ab_2",    "ddm1_2_G5",
  "ddm1_b_1",     "ddm1_2_G5",
  "ddm1_b_2",     "ddm1_2_G5",
  "ab_2",         "Col_0",
  "ddm1_2_G2",    "Col_0",
  "ddm1_2_G5",    "Col_0",
  "F2_ddm1",  "Col_0",
  "Col_0",  "Col_0"
)

# Step 1: Calculate the average control value per Geneid for each control condition
condition_avg <- m %>%
  filter(condition %in% unique(control_mapping$control)) %>%
  group_by(Geneid, condition) %>%
  summarize(control_mean = mean(value, na.rm = TRUE), .groups = 'drop') %>%
  rename(control = condition)
View(condition_avg)

# Step 2: Normalize the expression values by the control mean for the mapped samples
normalized_m <- m %>%
  left_join(control_mapping, by = "condition") %>%
  left_join(condition_avg, by = c("Geneid", "control")) %>%
  mutate(normalized_value = ifelse(!is.na(control), value - control_mean, NA)) %>%
  filter(!condition == "Col_0")

# View the resulting tibble with normalized values
View(normalized_m)

# import the superplot function
source("/groups/berger/user/pierre.bourguet/genomics/scripts/R/superplot_w_boxplot.R")

# Define the specific order vector
order_vector <- c("ab_2", "F2_ddm1", "F2_ddm1_a_2", "ddm1_2_G2", "ddm1_a_1", "ddm1_2_G5", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")

# Define the color palette
my_colors <- c("grey55", col_muted[c(3,5,3,5,3,5,5,5,5)]) # without
# rename the dataframe for compatibility with the function
names(normalized_m)[c(3,7)] <- c("unnormalized_value", "value")

# run the function and customize the plot
superplot_upTEs_ddm1 <- superplot_w_boxplot(normalized_m, rev(order_vector), rev(my_colors)) + 
  coord_flip(ylim=c(-0.5, 1.5)) +
  xlab("Genotype") + ylab("Transcript levels log2(ESF+1)") + labs(title="upregulated TEs")
superplot_upTEs_ddm1

# write svg output
svglite::svglite(filename = "figures/fig_rlog_TEs_ddm1_cdca7.svg", width = 80*mm_to_inches, height = 50*mm_to_inches)
set.seed(55)
barplot_nb_up_TEs_ddm1 + superplot_upTEs_ddm1 + patchwork::plot_layout(axes = 'collect') & theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- normalized_m %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
  
# are TEs up in ddm1 G2 also up in cdca7ab ? ####

########## correlate log2FC in ddm1 G2 and cdca7ab_2

# select TEs up in ddm1 G2, add a boolean to color TEs also up in cdca7ab_2
log2FC_df_ddm1_G2_up_TEs <- log2FC_df %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  mutate(up_in_cdca7ab = Geneid %in% up_TEs_cdca7_ids)

range <- c(-4, 12)

# scatterplot of ddm1 versus ab_2
ggplot(log2FC_df_ddm1_G2_up_TEs, aes(x = `ddm1_2_G2`, y = `ab_2`, color=up_in_cdca7ab)) +
  geom_point(size=(1.5 / log10(length(ddm1_G2_up_TEs)) ))  +
  ylim(range) + xlim(range) +
  geom_abline(intercept = 0, slope = 1, size = ggplot_line_width_1pt/2, linetype = "dotted", color="#4D4D4D", alpha=0.5) +
  labs(y = "cdca7ab-2", x="ddm1") +
  geom_smooth(alpha=0.8, size=0, fill="grey80", method=lm) +
  geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9) +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position="inside",
        legend.position.inside=c(.2,.9)) +
  stat_cor(method = "pearson", label.y = c(11,12), label.x = c(3,3), aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`, `* ")))

# scatterplot of ab_1 versus ab_2
ggplot(log2FC_df_ddm1_G2_up_TEs, aes(x = `ab_1`, y = `ab_2`, color=up_in_cdca7ab)) +
  geom_point(size=(1.5 / log10(length(ddm1_G2_up_TEs)) ))  +
  ylim(range) + xlim(range) +
  geom_abline(intercept = 0, slope = 1, size = ggplot_line_width_1pt/2, linetype = "dotted", color="#4D4D4D", alpha=0.5) +
  labs(y = "cdca7ab-2", x="cdca7ab-1") +
  geom_smooth(alpha=0.8, size=0, fill="grey80", method=lm) +
  geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9) +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position="inside",
        legend.position.inside=c(.2,.9)) +
  stat_cor(method = "pearson", label.y = c(11,12), label.x = c(3,3), aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`, `* ")))

########## correlate VST in ddm1 G2 and cdca7ab_2

# select TEs up in ddm1 G2, add a boolean to color TEs also up in cdca7ab_2
vsd_ddm1_G2_up_TEs <- vsd_avg %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  mutate(up_in_cdca7ab = Geneid %in% up_TEs_cdca7_ids)

range <- range(vsd_ddm1_G2_up_TEs %>%
                 select("ddm1_2_G2", "ab_2"))

# scatterplot of ddm1 versus ab_2
ggplot(vsd_ddm1_G2_up_TEs, aes(x = `ddm1_2_G2`, y = `ab_2`, color=up_in_cdca7ab)) +
  geom_point(size=(1.5 / log10(length(ddm1_G2_up_TEs))))  +
  ylim(range) + xlim(range) +
  geom_abline(intercept = 0, slope = 1, linewidth = ggplot_line_width_1pt/2, linetype = "dotted", color="#4D4D4D", alpha=0.5) +
  labs(y = "cdca7ab-2", x="ddm1") +
  geom_smooth(alpha=0.8, size=0, fill="grey80", method=lm) +
  geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9) +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position="inside",
        legend.position.inside=c(.2,.9)) +
  stat_cor(method = "pearson", label.y = c(9,10), aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`, `* ")))

############### ordered dot plot using log2FC

#### ab_1
log2FC_df_ddm1_G2_up_TEs_sort_ab1 <- log2FC_df_ddm1_G2_up_TEs %>%
  select(ddm1_2_G2, ab_1, ab_2, up_in_cdca7ab) %>%
  arrange(desc(ab_1)) %>%
  mutate(rank = row_number()) %>%  # Create a rank column based on the order
  pivot_longer(cols = -c(up_in_cdca7ab, rank), names_to = "genotype", values_to = "value") %>%
  filter(!genotype == "ab_2")

ordered_log2FC_dotplot_ab1 <- ggplot(log2FC_df_ddm1_G2_up_TEs_sort_ab1, aes(x = rank, y = value, color = genotype)) +
  geom_point(size=(0.1 / log10(length(ddm1_G2_up_TEs))))  +
  scale_color_manual(values = col_muted_2_replicates[c(7,5)]) +
  labs(x = "TEs upregulated in ddm1", y = "log2FC") +
  geom_hline(yintercept = 0, linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D") +
  ylim(-4,12) +
  
  geom_smooth(alpha=0.8, size=1, fill="grey80", method=lm) +
  geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9) +
  stat_cor(method = "spearman", label.y = c(9,10), aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`, `* "))) +
  
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_classic() +
  theme_horizontal_nature +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) ; ordered_log2FC_dotplot_ab1

#### ab_2

# order the data based on their log2FC in ddm1
log2FC_df_ddm1_G2_up_TEs_sort_ab2 <- log2FC_df_ddm1_G2_up_TEs %>%
  select(ddm1_2_G2, ab_1, ab_2, up_in_cdca7ab) %>%
  arrange(desc(ab_2)) %>%
  mutate(rank = row_number()) %>%  # Create a rank column based on the order
  pivot_longer(cols = -c(rank, up_in_cdca7ab), names_to = "genotype", values_to = "value") %>%
  filter(!genotype == "ab_1")

ordered_log2FC_dotplot_ab2 <- ggplot(log2FC_df_ddm1_G2_up_TEs_sort_ab2, aes(x = rank, y = value, color = genotype)) +
  geom_point(size=(0.1 / log10(length(ddm1_G2_up_TEs))))  +
  scale_color_manual(values = col_muted_2_replicates[c(8,5)]) +
  labs(x = "TEs upregulated in ddm1", y = "log2FC") +
  geom_hline(yintercept = 0, linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D") +
  ylim(-4,12) +
  
  geom_smooth(alpha=0.8, size=1, fill="grey80", method=lm) +
  geom_line(stat="smooth", method=lm, linewidth=0.2, alpha=0.9) +
  stat_cor(method = "spearman", label.y = c(9,10), aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "*`, `* "))) +
  
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_classic() +
  theme_horizontal_nature +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) ; ordered_log2FC_dotplot_ab2

# combine plots

ordered_log2FC_dotplot_ab <- ordered_log2FC_dotplot_ab1 + ordered_log2FC_dotplot_ab2 + plot_layout(guides = "collect", axis_titles = "collect")
ggsave("figures/ordered_dotplot_log2FC_ddm1_cdca7ab.svg", ordered_log2FC_dotplot_ab, width = 60*mm_to_inches, height = 40*mm_to_inches)

###### ordered dot plot without ddm1, coloring TEs up in cdca7ab_1 or cdca7ab_2

# select TEs up in ddm1 G2, add a boolean to color TEs also up in cdca7ab_2

log2FC_df_ddm1_G2_up_TEs_ab <- log2FC_df %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  mutate(up_in_cdca7ab = Geneid %in% c(up_TEs_cdca7[[2]])) %>%
  arrange(desc(ab_2)) %>%
  mutate(rank = row_number()) %>%  # Create a rank column based on the order
  select(Geneid, rank, ab_2, up_in_cdca7ab) %>%
  pivot_longer(cols = -c(Geneid, rank, up_in_cdca7ab), names_to = "genotype", values_to = "value")

ordered_log2FC_color_up_cdca7 <- ggplot(log2FC_df_ddm1_G2_up_TEs_ab, aes(x = rank, y = value, color = up_in_cdca7ab)) +
  geom_col(width=(0.1 / log10(length(ddm1_G2_up_TEs))), linewidth = ggplot_line_width_1pt/2) +
  geom_point(size=(0.1 / log10(length(ddm1_G2_up_TEs)))) +
  scale_color_manual(values = col_muted_2_replicates[c(11,8)]) +
  labs(x = "TEs upregulated in ddm1", y = "log2FC") +
  geom_hline(yintercept = 0, linewidth = ggplot_line_width_1pt/2, colour = "#4D4D4D") +
  ylim(range(log2FC_df_ddm1_G2_up_TEs_ab$value)) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_classic() +
  theme_horizontal_nature +
  theme(
    #legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) ; ordered_log2FC_color_up_cdca7

ggsave("figures/ordered_dotplot_log2FC_ddm1color_up_cdca7.svg", ordered_log2FC_color_up_cdca7, width = 50*mm_to_inches, height = 40*mm_to_inches)

############## calculate % of TEs of interest

# find % of TEs up in ddm1 G2 with 0 log2FC or less in cdca7ab_1 & cdca7ab_2
log2FC_df_ddm1_G2_up_TEs %>%
  filter(ab_1 <= 0 & ab_2 <= 0) %>%
  nrow(.)

# find % of TEs up in ddm1 G2, not up in cdca7_ab, with positive log2FC in cdca7ab_1 & cdca7ab_2
log2FC_df %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  # remove TEs not up in both cdca7_ab mutants
  filter(! Geneid %in% intersect(up_TEs_cdca7[[1]], up_TEs_cdca7[[2]])) %>%
  # Add a new column "log2FC_in_ab"
  mutate(log2FC_in_ab = case_when(
    ab_1 <= 0 & ab_2 <= 0 ~ "negative",
    ab_1 > 0 & ab_2 > 0 ~ "positive",
    TRUE ~ NA_character_
  )) %>%
  pull(log2FC_in_ab) %>%
  table(useNA = "always") %>%
  # divide by total
  prop.table()


 ############### find TEs up in ddm1 G2 with 0 log2FC or less in cdca7ab_1 & cdca7ab_2

log2FC_df_ddm1_G2_up_TEs %>%
  select(Geneid, ddm1_2_G2, ab_1, ab_2, up_in_cdca7ab) %>%
  filter(ab_1 > 0 & ab_2 > 0) %>%
  nrow()

# export as a tsv
dir.create("data", showWarnings = FALSE)
write.table(log2FC_df_ddm1_G2_up_TEs# %>%
              #pivot_longer(cols = -c(Geneid, up_in_cdca7ab), names_to = "genotype", values_to = "log2FC")
              , "data/log2FC_ddm1_G2_up_TEs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

########## write three different files with a deeptools-compatible format

# export one tsv in deeptools format for each category of TEs up in ddm1 G2
# First, create the three filtered datasets based on log2FC_in_ab values
ddm1_up_TEs_with_ab <- log2FC_df_ddm1_G2_up_TEs %>%
  mutate(log2FC_in_ab = case_when(
    ab_1 <= 0 & ab_2 <= 0 ~ "negative",
    ab_1 > 0 & ab_2 > 0 ~ "positive",
    TRUE ~ "discordant"
  ))

# Split into three datasets
negative_TEs <- ddm1_up_TEs_with_ab %>% filter(log2FC_in_ab == "negative")
positive_TEs <- ddm1_up_TEs_with_ab %>% filter(log2FC_in_ab == "positive")
discordant_TEs <- ddm1_up_TEs_with_ab %>% filter(log2FC_in_ab == "discordant")

# Function to create BED format from TE data
create_bed_format <- function(te_data, araport_annotations) {
  te_data %>%
    inner_join(araport_annotations, by = "Geneid") %>%
    filter(!str_ends(Geneid, "_AS")) %>%  # Remove entries ending with "_AS"
    mutate(
      chr = str_remove(chr, "^Chr"),  # Remove "Chr" prefix
      name = Geneid,
      score = 0
    ) %>%
    select(chr, start, end, name, score, strand) %>%
    arrange(chr, start)
}

# Create BED files for each category
negative_bed <- create_bed_format(negative_TEs, Araport11_annotations)
positive_bed <- create_bed_format(positive_TEs, Araport11_annotations)
discordant_bed <- create_bed_format(discordant_TEs, Araport11_annotations)

# Write the three BED files
write.table(negative_bed, "data/ddm1_G2_up_TEs_negative_log2FC_in_cdca7ab.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(positive_bed, "data/ddm1_G2_up_TEs_positive_log2FC_in_cdca7ab.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(discordant_bed, "data/ddm1_G2_up_TEs_discordant_log2FC_in_cdca7ab.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



############### ordered dot plot using VST

# order the data based on their log2FC in ddm1
df <- vsd_ddm1_G2_up_TEs %>% select(Geneid, ddm1_2_G2, ab_1, ab_2) %>%
  arrange(desc(ab_2)) %>%
  pivot_longer(cols = -Geneid, names_to = "genotype", values_to = "value") %>%
  mutate(Geneid = factor(Geneid, levels = unique(Geneid)))

ggplot(df %>% filter(!genotype == "ab_1"), aes(x = Geneid, y = value, color = genotype)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = col_muted_2_replicates[c(7,5)]) +  # Customize colors as needed
  labs(x = "TEs upregulated in ddm1", y = "VST") +
  geom_hline(yintercept = 0) +  # Add a horizontal line at y = 0
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ylim(min(df$value), max(df$value))

#
# superfamilies & families of TEs up in cdca7 vs TEs up only in ddm1 ####

# THE LATEST CODE USED IN THE PAPER IS IN THE "23.10_mC_analysis_cdca7.R" SCRIPT

# Merge 3 cluster with TEs based on Geneid
merged_cdca7_log2FC_clusters_TEs <- as_tibble(merged_cdca7_log2FC_clusters) %>%
  dplyr::rename(Geneid = V4, chr = V1, start = V2, end = V3, strand = V6, Cluster = cluster) %>%
  dplyr::select(Geneid, chr, Cluster, strand) %>%
  distinct() %>%
  left_join(TEs, by = "Geneid")

# TE superfamilies - cluster enrichment

TE_superfamily_heatmap <- function (clusters) {
  require("pheatmap")
  TE_families <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_Transposable_Elements.txt", head=T, sep="\t", quote="", comment.char="")
  names(TE_families)[1] <- "Geneid"
  clusters_families <- merge(clusters, TE_families, by="Geneid")
  z <- table(TE_families$Transposon_Super_Family) / nrow(TE_families) * 100 # genome % of each TE superfamily
  for (i in unique(clusters$Cluster)) {
    z <- dplyr::bind_rows(z, # this function binds rows based on column names and importantly, can handle a row with missing columns
                          table(clusters_families$Transposon_Super_Family[clusters_families$Cluster == i]) / sum(clusters_families$Cluster == i) * 100)
  }
  row.names(z) <- c("genome", paste0("", unique(clusters$Cluster)))
  pheatmap(z, cellwidth=20, cellheight=20, cluster_cols = F, cluster_rows = F, show_rownames=T, fontsize=10, main="% of TE superfamilies"
           , color=col_TE_families, border_color = "black", breaks= seq(0, max(z, na.rm = T), length.out=100) )
}

TE_superfamily_heatmap(merged_cdca7_log2FC_clusters_TEs)

# save the heatmap
svglite::svglite(filename = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/figures/TE_superfamily_heatmap_at_ddm1_upTEs_with_AS_cdca7ab_log2FC_x3_clusters.svg"
                 , width = 210*mm_to_inches, height = 100*mm_to_inches)
TE_superfamily_heatmap(merged_cdca7_log2FC_clusters_TEs)
dev.off()

# TE families

TE_family_heatmap <- function(clusters, max_cols_per_plot = 15) {
  require("pheatmap")
  TE_families <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_Transposable_Elements.txt", 
                            head=T, sep="\t", quote="", comment.char="")
  names(TE_families)[1] <- "Geneid"
  clusters_families <- merge(clusters, TE_families, by="Geneid")
  
  # Get families that actually appear in clusters
  families_in_clusters <- unique(clusters_families$Transposon_Family)
  
  # Filter genome data to only include families found in clusters
  TE_families_filtered <- TE_families[TE_families$Transposon_Family %in% families_in_clusters, ]
  
  # Now build the table with only relevant families
  z <- table(TE_families_filtered$Transposon_Family) / nrow(TE_families_filtered) * 100
  
  for (i in unique(clusters$Cluster)) {
    z <- dplyr::bind_rows(z,
                          table(clusters_families$Transposon_Family[clusters_families$Cluster == i]) / 
                            sum(clusters_families$Cluster == i) * 100)
  }
  row.names(z) <- c("genome", paste0("", unique(clusters$Cluster)))
  
  # Rest of the function remains the same...
  z_matrix <- as.matrix(z)
  n_cols <- ncol(z_matrix)
  n_plots <- ceiling(n_cols / max_cols_per_plot)
  breaks <- seq(0, max(z_matrix, na.rm = T), length.out = 100)
  
  plot_list <- list()
  for(plot_num in 1:n_plots) {
    start_col <- (plot_num - 1) * max_cols_per_plot + 1
    end_col <- min(plot_num * max_cols_per_plot, n_cols)
    z_subset <- z_matrix[, start_col:end_col, drop = FALSE]
    
    plot_list[[plot_num]] <- pheatmap(z_subset, 
                                      cellwidth = 10, cellheight = 10, 
                                      cluster_cols = FALSE, cluster_rows = FALSE, 
                                      show_rownames = TRUE, fontsize = 10, 
                                      main = paste0("% of TE families (Part ", plot_num, "/", n_plots, ")"),
                                      color = col_TE_families, border_color = "black", 
                                      breaks = breaks, na_col = "grey90")
  }
  return(plot_list)
}

# This will create multiple plots automatically
plots <- TE_family_heatmap(merged_cdca7_log2FC_clusters_TEs, max_cols_per_plot = 61)

plots[[1]]
plots[[2]]
plots[[3]]

#
# pericentromeric location of TEs up in cdca7 vs TEs up only in ddm1 ####

# proportion in peri and arms according to H3K9me2-rich regions (Bernatavichute 2008 Plos ONE)
cluster_peri_arm_proportion <- function(df, clusters, type) {
  library("GenomicRanges")

  # Read coordinate files
  peri <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_het_peri_coordinates.tsv", header=F)
  arms <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_euchro_coordinates.tsv", header=F)
  
  # Convert to GRanges objects
  peri_gr <- GRanges(seqnames = peri$V1, ranges = IRanges(start = peri$V2, end = peri$V3))
  arms_gr <- GRanges(seqnames = arms$V1, ranges = IRanges(start = arms$V2, end = arms$V3))
  
  # Merge df and clusters
  y <- merge(clusters, df, by="Geneid")
  
  # Initialize proportions dataframe
  proportions <- data.frame(cluster = c("genome", "genome", rep(unique(sort(clusters$Cluster)), each = 2)),
                            location = rep(c("arms", "peris"), length(unique(clusters$Cluster)) + 1),
                            value = NA)
  
  # Read genome annotation
  genome <- read.delim(file="/groups/berger/user/pierre.bourguet/genomics/Araport11/Araport11_GFF3_PCG_TE_TEG.SAF", 
                       head=T, sep="\t", quote="", comment.char="")
  genome <- genome[genome$Type == type,]
  
  # Convert genome to GRanges
  genome_gr <- GRanges(seqnames = genome$Chr,
                       ranges = IRanges(start = genome$Start, end = genome$End))
  
  # Calculate genome-wide overlaps
  genome_arms_overlaps <- findOverlaps(genome_gr, arms_gr)
  genome_peris_overlaps <- findOverlaps(genome_gr, peri_gr)
  
  proportions$value[proportions$cluster=="genome" & proportions$location=="arms"] <- length(genome_arms_overlaps)
  proportions$value[proportions$cluster=="genome" & proportions$location=="peris"] <- length(genome_peris_overlaps)
  
  # Calculate cluster-specific overlaps
  for (i in unique(clusters$Cluster)) {
    y_cluster <- y[y$Cluster == i,]
    
    # Convert cluster data to GRanges
    y_cluster_gr <- GRanges(seqnames = paste0("chr", y_cluster$chr),
                            ranges = IRanges(start = y_cluster$start, end = y_cluster$end))
    
    # Find overlaps
    cluster_arms_overlaps <- findOverlaps(y_cluster_gr, arms_gr)
    cluster_peris_overlaps <- findOverlaps(y_cluster_gr, peri_gr)
    
    proportions$value[proportions$cluster==i & proportions$location=="arms"] <- length(cluster_arms_overlaps)
    proportions$value[proportions$cluster==i & proportions$location=="peris"] <- length(cluster_peris_overlaps)
  }
  
  # Create plot
  ggplot(proportions, aes(fill=location, y=value, x=cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    theme_classic() +
    ggtitle(paste0("distribution of ", type)) +
    scale_fill_manual(values=c("grey86", "grey42")) + 
    labs(y="proportion")
}

plot_peri_arm_freq <- cluster_peri_arm_proportion(df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
                            clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
                            type = "transposable_element")
plot_peri_arm_freq + theme_horizontal_nature

# save the plot
svglite::svglite(filename = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/figures/proportion_arm_peri_of_cdca7ab_log2FC_x3_clusters.svg"
                 , width = 30*mm_to_inches, height = 40*mm_to_inches)
plot_peri_arm_freq + theme_horizontal_nature# + ylim(0,1)
dev.off()

# chi square test to calculate enrichment in pericentromeres
pairwise_chisq_peri_enrichment <- function(df, clusters, type, correction_method = "fdr") {
  library("GenomicRanges")
  
  # Read coordinate files
  peri <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_het_peri_coordinates.tsv", header=F)
  arms <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_euchro_coordinates.tsv", header=F)
  
  # Convert to GRanges objects
  peri_gr <- GRanges(seqnames = peri$V1, ranges = IRanges(start = peri$V2, end = peri$V3))
  arms_gr <- GRanges(seqnames = arms$V1, ranges = IRanges(start = arms$V2, end = arms$V3))
  
  # Merge df and clusters
  y <- merge(clusters, df, by="Geneid")
  
  # Get unique clusters
  cluster_names <- unique(sort(clusters$Cluster))
  
  # Calculate counts for each cluster
  cluster_counts <- data.frame(cluster = cluster_names, arms = 0, peris = 0)
  
  for (i in cluster_names) {
    y_cluster <- y[y$Cluster == i,]
    
    # Convert cluster data to GRanges
    y_cluster_gr <- GRanges(seqnames = paste0("chr", y_cluster$chr),
                            ranges = IRanges(start = y_cluster$start, end = y_cluster$end))
    
    # Find overlaps
    cluster_arms_overlaps <- findOverlaps(y_cluster_gr, arms_gr)
    cluster_peris_overlaps <- findOverlaps(y_cluster_gr, peri_gr)
    
    cluster_counts$arms[cluster_counts$cluster == i] <- length(cluster_arms_overlaps)
    cluster_counts$peris[cluster_counts$cluster == i] <- length(cluster_peris_overlaps)
  }
  
  # Perform pairwise chi-square tests
  results <- list()
  
  # Store all p-values for correction  
  all_p_values <- c()  
  comparison_names <- c()  
  
  # Perform pairwise chi-square tests  
  results <- list()
  
  for (i in 1:(length(cluster_names)-1)) {
    for (j in (i+1):length(cluster_names)) {
      cluster1 <- cluster_names[i]
      cluster2 <- cluster_names[j]
      
      # Create contingency table
      contingency_table <- matrix(c(
        cluster_counts$arms[cluster_counts$cluster == cluster1],
        cluster_counts$peris[cluster_counts$cluster == cluster1],
        cluster_counts$arms[cluster_counts$cluster == cluster2],
        cluster_counts$peris[cluster_counts$cluster == cluster2]
      ), nrow = 2, byrow = TRUE,
      dimnames = list(c(cluster1, cluster2), c("arms", "peris")))
      
      # Perform chi-square test
      chisq_result <- chisq.test(contingency_table)
      
      # Calculate proportions
      prop1_peri <- cluster_counts$peris[cluster_counts$cluster == cluster1] / 
        (cluster_counts$arms[cluster_counts$cluster == cluster1] + 
           cluster_counts$peris[cluster_counts$cluster == cluster1])
      prop2_peri <- cluster_counts$peris[cluster_counts$cluster == cluster2] / 
        (cluster_counts$arms[cluster_counts$cluster == cluster2] + 
           cluster_counts$peris[cluster_counts$cluster == cluster2])
      
      # Store results
      comparison_name <- paste(cluster1, "vs", cluster2, sep = "_")
      all_p_values <- c(all_p_values, chisq_result$p.value)  
      comparison_names <- c(comparison_names, comparison_name)  
      
      results[[comparison_name]] <- list(  
        contingency_table = contingency_table,  
        chi2_statistic = chisq_result$statistic,  
        p_value_raw = chisq_result$p.value,  
        prop_peri_cluster1 = prop1_peri,  
        prop_peri_cluster2 = prop2_peri,  
        difference = prop1_peri - prop2_peri
      )
    }
  }
  
  # Apply multiple testing correction  
  adjusted_p_values <- p.adjust(all_p_values, method = correction_method)  
  
  # Add adjusted p-values to results  
  for (i in 1:length(comparison_names)) {  
    results[[comparison_names[i]]]$p_value_adjusted <- adjusted_p_values[i]  
  }  
  
  return(results)
  
}

chisq_results <- pairwise_chisq_peri_enrichment(
  df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
  clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
  type = "transposable_element"
)

# Print results
for (comparison in names(chisq_results)) {
  cat("\n", comparison, ":\n")  
  cat("  Chi-square statistic:", round(chisq_results[[comparison]]$chi2_statistic, 4), "\n")  
  cat("  P-value (raw):", format(chisq_results[[comparison]]$p_value_raw, scientific = TRUE, digits = 3), "\n")  
  cat("  P-value (FDR-adjusted):", format(chisq_results[[comparison]]$p_value_adjusted, scientific = TRUE, digits = 3), "\n")
  cat("  Pericentromeric proportions:", 
      round(chisq_results[[comparison]]$prop_peri_cluster1, 3), "vs", 
      round(chisq_results[[comparison]]$prop_peri_cluster2, 3), "\n")
  cat("  Contingency table:\n")
  print(chisq_results[[comparison]]$contingency_table)
}

# function to plot the number of TEs in each cluster across chromosomes, with pericentromeric regions and centromeres highlighted
plot_te_density_by_cluster <- function(df, clusters, type, window_size = 500000, colors = NULL, highlight_peri = FALSE, highlight_centromeres = FALSE) {
  library("GenomicRanges")
  library("ggplot2")
  library("dplyr")
  
  # Merge df and clusters
  y <- merge(clusters, df, by="Geneid")
  
  # Define chromosome lengths for Arabidopsis (TAIR10)
  chr_lengths <- c(chr1 = 30427671, chr2 = 19698289, chr3 = 23459830, 
                   chr4 = 18585056, chr5 = 26975502)
  
  # Read pericentromeric coordinates
  peri <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_het_peri_coordinates.tsv", header=F)
  
  # Read centromere coordinates
  centromeres <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_cenH3_coordinates_GSE88907_nochr.tsv", 
                            header = FALSE, col.names = c("chr", "start", "end"))
  # Add chr prefix to match other data
  centromeres$chr <- paste0("chr", centromeres$chr)
  
  # Create windows for all chromosomes
  windows_list <- list()
  cumulative_pos <- 0
  chr_boundaries <- c(0)
  
  for (chr in names(chr_lengths)) {
    n_windows <- ceiling(chr_lengths[chr] / window_size)
    starts <- seq(1, chr_lengths[chr], by = window_size)
    ends <- pmin(starts + window_size - 1, chr_lengths[chr])
    
    windows_df <- data.frame(
      chr = chr,
      start = starts,
      end = ends,
      window_id = paste0(chr, "_", 1:length(starts)),
      genome_pos = cumulative_pos + starts + (window_size/2)  # midpoint for plotting
    )
    
    windows_list[[chr]] <- windows_df
    cumulative_pos <- cumulative_pos + chr_lengths[chr]
    chr_boundaries <- c(chr_boundaries, cumulative_pos)
  }
  
  all_windows <- do.call(rbind, windows_list)
  
  # Convert pericentromeric coordinates to genome positions
  peri_regions <- data.frame()
  cumulative <- 0
  for (chr in names(chr_lengths)) {
    chr_peri <- peri[peri$V1 == chr, ]
    if (nrow(chr_peri) > 0) {
      chr_peri$genome_start <- cumulative + chr_peri$V2
      chr_peri$genome_end <- cumulative + chr_peri$V3
      peri_regions <- rbind(peri_regions, chr_peri)
    }
    cumulative <- cumulative + chr_lengths[chr]
  }
  
  # Convert centromere coordinates to genome positions
  centromere_regions <- data.frame()
  cumulative <- 0
  for (chr in names(chr_lengths)) {
    chr_centromeres <- centromeres[centromeres$chr == chr, ]
    if (nrow(chr_centromeres) > 0) {
      chr_centromeres$genome_start <- cumulative + chr_centromeres$start
      chr_centromeres$genome_end <- cumulative + chr_centromeres$end
      centromere_regions <- rbind(centromere_regions, chr_centromeres)
    }
    cumulative <- cumulative + chr_lengths[chr]
  }
  
  # Convert windows to GRanges
  windows_gr <- GRanges(seqnames = all_windows$chr,
                        ranges = IRanges(start = all_windows$start, end = all_windows$end))
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Process each cluster
  for (cluster in unique(clusters$Cluster)) {
    y_cluster <- y[y$Cluster == cluster,]
    
    # Convert cluster data to GRanges
    y_cluster_gr <- GRanges(seqnames = paste0("chr", y_cluster$chr),
                            ranges = IRanges(start = y_cluster$start, end = y_cluster$end))
    
    # Find overlaps between TEs and windows
    overlaps <- findOverlaps(y_cluster_gr, windows_gr)
    
    # Count TEs per window
    te_counts <- table(subjectHits(overlaps))
    
    # Create results for this cluster
    cluster_results <- all_windows
    cluster_results$cluster <- cluster
    cluster_results$te_count <- 0
    cluster_results$te_count[as.numeric(names(te_counts))] <- as.numeric(te_counts)
    
    results <- rbind(results, cluster_results)
  }
  
  # Fix chromosome boundary and label calculations
  chr_boundaries_plot <- chr_boundaries[-length(chr_boundaries)]  # Remove last boundary
  
  # Calculate chromosome center positions correctly
  chr_centers <- numeric(length(chr_lengths))
  cumulative <- 0
  for (i in 1:length(chr_lengths)) {
    chr_centers[i] <- cumulative + chr_lengths[i]/2
    cumulative <- cumulative + chr_lengths[i]
  }
  
  chr_labels <- data.frame(
    chr = names(chr_lengths),
    pos = chr_centers
  )
  
  # Create adaptive y-axis label based on window size
  if (window_size >= 1000000) {
    window_label <- paste0(window_size / 1000000, "Mb")
  } else if (window_size >= 1000) {
    window_label <- paste0(window_size / 1000, "kb")
  } else {
    window_label <- paste0(window_size, "bp")
  }
  
  # Start building the plot
  p <- ggplot(results, aes(x = genome_pos, y = te_count))
  
  # Add pericentromeric highlighting FIRST (behind bars) if requested
  if (highlight_peri && nrow(peri_regions) > 0) {
    p <- p + geom_rect(data = peri_regions, 
                       aes(xmin = genome_start, xmax = genome_end, 
                           ymin = -Inf, ymax = Inf),
                       fill = "grey70", alpha = 0.3, inherit.aes = FALSE)
  }
  
  # Add centromere highlighting (on top of pericentromeric regions) if requested
  if (highlight_centromeres && nrow(centromere_regions) > 0) {
    p <- p + geom_rect(data = centromere_regions, 
                       aes(xmin = genome_start, xmax = genome_end, 
                           ymin = -Inf, ymax = Inf),
                       fill = "black", alpha = 0.5, inherit.aes = FALSE)
  }
  
  # Add bars on top of rectangles
  p <- p + geom_col(aes(fill = cluster), width = window_size * 0.8, alpha = 0.8) +
    facet_wrap(~cluster, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = chr_boundaries_plot[-1], linetype = "dashed", alpha = 0.5, size=ggplot_line_width_1pt/4) +
    scale_x_continuous(
      breaks = chr_labels$pos,
      labels = chr_labels$chr,
      expand = c(0.01, 0)
    ) +
    theme_classic() +
    labs(
      x = "Chromosome",
      y = paste0("Number of ", type, " per ", window_label, " window"),
      title = paste0("Distribution of ", type, " across genome by cluster")
    ) +
    theme_horizontal_nature +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"  # Remove legend since clusters are in facets
    )
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  return(p)
}

# With pericentromeric regions highlighted
plot_te_density_by_cluster(df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
                           clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
                           type = "transposable_element",
                           window_size = 250000,
                           highlight_peri = TRUE)

# With custom colors and highlighting
plot_te_density_by_cluster(df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
                           clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
                           type = "transposable_element",
                           window_size = 500000,
                           colors = c("#E68B50", "#CB514A", "#71372E"),
                           highlight_peri = TRUE)

# save the plot
svglite::svglite(filename = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/figures/chr_location_of_ddm1_upTEs_with_AS_cdca7ab_log2FC_x3_clusters.svg"
                 , width = 100*mm_to_inches, height = 75*mm_to_inches)
plot_te_density_by_cluster(df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
                           clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
                           type = "transposable_element",
                           window_size = 500000,
                           colors = c("#E68B50", "#CB514A", "#71372E"),
                           highlight_peri = TRUE,
                           highlight_centromeres = TRUE)
dev.off()

# THE LATEST CODE FOR TE FAMILY & SUPERFAMILY ANALYSIS (USED IN THE PAPER) IS IN THE "23.10_mC_analysis_cdca7.R" SCRIPT

#
# distance to centromere of TEs up in cdca7 vs TEs up only in ddm1 ####

cluster_centromere_distance <- function(df, clusters, type, ylim=NULL) {
  library("GenomicRanges")
  library("ggplot2")
  library("ggbeeswarm")
  library("dplyr")
  
  # Read centromere coordinates
  centromeres <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_cenH3_coordinates_GSE88907_nochr.tsv", 
                            header = FALSE, col.names = c("chr", "start", "end"))
  
  # Calculate centromere midpoints
  centromeres$midpoint <- (centromeres$start + centromeres$end) / 2
  
  # Read all TEs for genome background
  all_TEs <- read.delim("/groups/berger/user/pierre.bourguet/genomics/Araport11/TAIR10_Transposable_Elements.txt",
                        header = TRUE, sep = "\t", quote = "", comment.char = "", 
                        col.names = c("Geneid", "sense", "start", "end", "family", "superfamily"))
  
  # Extract chromosome from Geneid (assuming format like AT1TE12345)
  all_TEs$chr <- as.numeric(substr(all_TEs$Geneid, 3, 3))
  
  # Calculate TE midpoints for all TEs
  all_TEs$midpoint <- (all_TEs$start + all_TEs$end) / 2
  
  # Calculate distances to centromeres for all TEs
  all_TEs$distance_to_centromere <- NA
  for(i in 1:nrow(all_TEs)) {
    chr <- all_TEs$chr[i]
    te_midpoint <- all_TEs$midpoint[i]
    centromere_midpoint <- centromeres$midpoint[centromeres$chr == chr]
    
    if(length(centromere_midpoint) > 0) {
      all_TEs$distance_to_centromere[i] <- abs(te_midpoint - centromere_midpoint)
    }
  }
  
  # Create genome background data
  genome_data <- data.frame(
    Geneid = all_TEs$Geneid,
    cluster = "genome",
    chr = all_TEs$chr,
    start = all_TEs$start,
    end = all_TEs$end,
    midpoint = all_TEs$midpoint,
    distance_to_centromere = all_TEs$distance_to_centromere
  )
  
  # Merge df and clusters for cluster-specific TEs
  y <- merge(clusters, df, by = "Geneid")
  
  # Extract chromosome from Geneid for cluster TEs
  y$chr_num <- as.numeric(substr(y$Geneid, 3, 3))
  
  # Calculate TE midpoints for cluster TEs
  y$midpoint <- (y$start + y$end) / 2
  
  # Calculate distances to centromeres for cluster TEs
  y$distance_to_centromere <- NA
  for(i in 1:nrow(y)) {
    chr <- y$chr_num[i]
    te_midpoint <- y$midpoint[i]
    centromere_midpoint <- centromeres$midpoint[centromeres$chr == chr]
    
    if(length(centromere_midpoint) > 0) {
      y$distance_to_centromere[i] <- abs(te_midpoint - centromere_midpoint)
    }
  }
  
  # Create cluster data
  cluster_data <- data.frame(
    Geneid = y$Geneid,
    cluster = y$Cluster,
    chr = y$chr_num,
    start = y$start,
    end = y$end,
    midpoint = y$midpoint,
    distance_to_centromere = y$distance_to_centromere
  )
  
  # Combine all data
  plot_data <- rbind(genome_data, cluster_data)
  
  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$distance_to_centromere), ]
  
  # Create cluster-only data for geom_quasirandom
  cluster_only_data <- plot_data[plot_data$cluster != "genome", ]
  
  # Set up colors
  n_clusters <- length(unique(cluster_only_data$cluster))
  cluster_colors <- c("#E68B50", "#CB514A", "#71372E")[1:n_clusters]
  names(cluster_colors) <- unique(sort(cluster_only_data$cluster))
  all_colors <- c("genome" = "grey50", cluster_colors)
  
  # Create plot
  centromere_plot <- ggplot(plot_data, aes(x = cluster, y = log10(distance_to_centromere), color = cluster)) +
    {if(nrow(cluster_only_data) > 0) 
      geom_quasirandom(data = cluster_only_data, alpha = 1, size = 0.01, color = "grey90")} +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, linewidth = ggplot_line_width_1pt/2) +
    scale_color_manual(values = all_colors) +
    {if(!is.null(ylim)) coord_cartesian(ylim = ylim)} +
    theme_classic() +
    labs(x = "Cluster", 
         y = "Distance to centromere (log10(bp))") +
    ggtitle(paste0("Distance to centromere - ", type)) +
    theme(legend.position = "none",
          plot.margin = margin(5, 1, 5, 1),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(plot = centromere_plot, data = plot_data))
}

# Execute the function
result <- cluster_centromere_distance(
  df = merged_cdca7_log2FC_clusters_TEs %>% select(-Cluster), 
  clusters = merged_cdca7_log2FC_clusters_TEs[,c(1,3)], 
  type = "transposable_element",
  ylim= c(4.5, 7.2)  # Adjust ylim as needed
)

# Display plot
print(result$plot)

# Access data
distance_data <- result$data
head(distance_data)

# Dunn post hoc tests
perform_dunn_test_fsa <- function(distance_data) {
  library(FSA)
  library(dplyr)
  
  # Filter and transform data
  cluster_data <- distance_data %>%
    filter(cluster != "genome") %>%
    filter(!is.na(distance_to_centromere)) %>%
    filter(distance_to_centromere > 0) %>%
    mutate(log_distance = log10(distance_to_centromere))
  
  # Kruskal-Wallis test first
  kw_test <- kruskal.test(log_distance ~ cluster, data = cluster_data)
  cat("Kruskal-Wallis test:\n")
  print(kw_test)
  cat("\n")
  
  # Dunn test with BH correction
  dunn_test <- dunnTest(log_distance ~ cluster, 
                        data = cluster_data,
                        method = "bh")
  
  cat("Dunn test results (BH corrected):\n")
  print(dunn_test)
  
  return(list(
    kruskal_wallis = kw_test,
    dunn_test = dunn_test,
    test_data = cluster_data
  ))
}

dunn_results_fsa <- perform_dunn_test_fsa(distance_data)

# save the plot with svglite
svglite::svglite(filename = "/groups/berger/user/pierre.bourguet/genomics/scripts/3_prime_tag-seq/pipeline_vikas/04_output/tagseq_01_cdca7_mutants_AtRTD3_ATTE_150bp_5M_min_50bp/07_analysis/figures/distance_to_centromere_cdca7ab_log2FC_x3_clusters.svg"
                 , width = 25*mm_to_inches, height = 50*mm_to_inches)
result$plot + theme_horizontal_nature + theme(legend.position = "none")
dev.off()

# cdca7 long mutants: number of up TEs, heatmap and superplots ####

# filter samples
nb_up_TEs_df_a_long_no_ddm1 <- nb_up_TEs_df %>% filter(str_detect(condition, "long|a_1|a_2|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

nb_up_TEs_df_a_long_only_ab <- nb_up_TEs_df %>% filter(str_detect(condition, "a_long_b|ab")) %>%
  filter(!str_detect(condition, "ddm1")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0)) %>%
  mutate(condition = factor(condition, levels = rev(c("Col_0", "ab_1", "ab_2", "a_long_b")))) %>%
  arrange(condition)
  
nb_up_TEs_df_a_long_w_ddm1 <- nb_up_TEs_df %>% filter(str_detect(condition, "ddm1") & str_detect(condition, "G5|long|ab")) %>%
  bind_rows(tibble(condition = "Col_0", Count = 0))

# reorder levels
# nb_up_TEs_df_w_ddm1$condition <- factor(nb_up_TEs_df_w_ddm1$condition, levels = rev(c("Col_0", "ab_1", "ab_2", "ddm1_2_G5", "ddm1_a_1", "F2_ddm1_a_2", "ddm1_b_1", "ddm1_b_2", "ddm1_ab_1", "ddm1_ab_2")))

barplot_nb_up_TEs_a_long <- ggplot(rbind(nb_up_TEs_df_a_long_no_ddm1, nb_up_TEs_df_a_long_w_ddm1), aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  #scale_fill_manual(values = rev(c(col_muted[c(10,4,4,3,5)], rep(col_muted[5], 5)))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long

barplot_nb_up_TEs_a_long_b <- ggplot(nb_up_TEs_df_a_long_only_ab, aes(y = condition, x = Count, fill=condition)) +
  geom_col() +
  labs(title = "Number of upregulated TEs in each genotype", x = "Count", y = "Genotype") +
  scale_fill_manual(values = rev(c(col_muted[c(10,4,4,6)]))) +
  theme_minimal() + theme(legend.position = "none") ; barplot_nb_up_TEs_a_long_b

## heatmap at TEs up in cdca7a/b
cdca7_ab_up_TEs <- intersect(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))
DEG_heatmap(cdca7_ab_up_TEs, "cdca7_ab_up_TEs_cdca7a_long_mutants"
            , z = rld_df %>% dplyr::select(matches("long_b|ab|Geneid") & -matches("ddm1"))
            , n = "rlog")
wide_rld_avg

# superplot
#cdca7_ab_up_TEs_union <- union(up(all_res$res_ab_1, "TEs"), up(all_res$res_ab_2, "TEs"))

a_long_b_rld <- rld_df %>%
  dplyr::filter(Geneid %in% cdca7_ab_up_TEs) %>%
  dplyr::select(matches("long_b|ab|Geneid|Col_0") & -matches("ddm1")) %>%
  tidyr::pivot_longer(cols = -Geneid, names_to = "sample", values_to = "value") %>%
  dplyr::mutate(condition=stringr::str_replace(sample, "_R[123]", ""))

superplot_a_long_b <- superplot_w_boxplot(data = a_long_b_rld, condition_order = rev(c("Col_0", "ab_1", "ab_2", "a_long_b"))
                                          , colors = rev(c("grey55", col_muted[c(4,4,6)]))) + 
  coord_flip(ylim=c(1.2,9)) +
  xlab("Genotype") + ylab("Transcript levels (rlog)") + labs(title="upregulated TEs")

# write svg output
svglite::svglite(filename = "figures/superplot_a_long_b.svg", width = 3, height = 1.5)
barplot_nb_up_TEs_a_long_b + superplot_a_long_b + plot_layout(axes = "collect", widths = c(2,3)) #& theme_horizontal_nature
dev.off()

# Tukey HSD
# Calculate median values for each sample
medians <- a_long_b_rld %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarize(value = median(value))

# Perform ANOVA
anova_result <- aov(value ~ condition, data = medians)
summary(anova_result)

# Perform Tukey's HSD test
tukey_result <- HSD.test(anova_result, "condition")
print(tukey_result)

#
# how do cdca7-upregulated TEs compare with those upregulated in h2a.w suvh456 & ddm1? ####

# import the data
w_s456_RPM <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_c23_s456_a56_m1_10d_in_vitro_qseq/RPM.tsv"), sep="\t", dec=".", quote=F)
w_s456_RPM_AS <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_c23_s456_a56_m1_10d_in_vitro_qseq/RPM_AS.tsv"), sep="\t", dec=".", quote=F)

# add a _AS suffix to the Geneid column and append both dataframes
w_s456_RPM_AS$Geneid <- paste0(w_s456_RPM_AS$Geneid, "_AS")
w_s456_RPM_merged <- rbind(w_s456_RPM, w_s456_RPM_AS)

#### first heatmap: select TEs up in ddm1 in this tagseq01 cdca7 dataset, then plot a heatmap to compare cdca7 and w_s456 datasets
# after z-score normalization

# select TEs up in ddm1
ddm1_G2_up_TEs <- up(all_res$res_ddm1_2_G2, "TEs")

# select samples of interest
m_tagseq01 <- RPM_merged_avg %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  select(Col_0, ddm1_2_G2, ddm1_2_G5, ab_1, ab_2, a_1, a_2) %>%
  as.matrix(.) %>%
  # apply row-wise z-score normalization
  t() %>% scale() %>% t()

m_w_s456 <- w_s456_RPM_merged %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  select(Col_0, ddm1, met1, suvh456, w_suvh456) %>%
  as.matrix(.) %>%
  # apply row-wise z-score normalization
  t() %>% scale() %>% t()

# make a heatmap using complexheatmap
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
tagseq01_h <- Heatmap(m_tagseq01,
        name = "z-score",
        col = col_zscore,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15,
        km = 4, # number of k-means clusters
        row_km_repeats = 100,
        border = TRUE
)
w_s456_h <- Heatmap(m_w_s456,
        name = "z-score",
        col = col_zscore,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15,
        km = 4, # number of k-means clusters
        row_km_repeats = 100,
        border = TRUE
)
tagseq01_h + w_s456_h
w_s456_h + tagseq01_h


#============ now using clusters of ddm1 up TEs from Nat Plants revision 1

# First, merge your expression data with cluster information
m_w_s456 <- w_s456_RPM_merged %>%
  filter(Geneid %in% merged_cdca7_log2FC_clusters$V4) %>%  # Use V4 which contains the TE IDs
  left_join(merged_cdca7_log2FC_clusters %>% select(V4, cluster), by = c("Geneid" = "V4")) %>%
  arrange(cluster) %>%  # Sort by cluster to group them together
  select(Col_0, ddm1, met1, atxr56, cmt23, suvh456, w_suvh456) %>%
  as.matrix(.) %>%
  # apply row-wise z-score normalization
  t() %>% scale() %>% t()

# Create cluster annotation
cluster_info <- w_s456_RPM_merged %>%
  filter(Geneid %in% merged_cdca7_log2FC_clusters$V4) %>%
  left_join(merged_cdca7_log2FC_clusters %>% select(V4, cluster), by = c("Geneid" = "V4")) %>%
  arrange(cluster) %>%
  pull(cluster)

# Create row annotation for clusters
row_ha <- rowAnnotation(
  Cluster = cluster_info,
  col = list(Cluster = c("cluster_1" = "#E31A1C", 
                         "cluster_2" = "#1F78B4", 
                         "cluster_3" = "#33A02C"))
)

# Create the heatmap with your predefined clusters
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B"))

w_s456_h <- Heatmap(m_w_s456,
                    name = "z-score",
                    col = col_zscore,
                    cluster_rows = FALSE,  # Don't cluster since we're using predefined clusters
                    cluster_columns = FALSE,
                    #width = unit(15, "cm"), 
                    #height = unit(15, "cm"),
                    border = TRUE,
                    left_annotation = row_ha,  # Add cluster annotation
                    row_split = cluster_info,  # Split rows by cluster
                    row_title_rot = 0
)
w_s456_h


#============= with boxplots ===

boxplot_data <- w_s456_RPM_merged %>%
  filter(Geneid %in% merged_cdca7_log2FC_clusters$V4) %>%
  left_join(merged_cdca7_log2FC_clusters %>% select(V4, cluster), by = c("Geneid" = "V4")) %>%
  select(Geneid, cluster, Col_0, ddm1, met1, atxr56, cmt23, suvh456, w_suvh456) %>%
  pivot_longer(cols = c(Col_0, ddm1, met1, atxr56, cmt23, suvh456, w_suvh456), 
               names_to = "Sample", 
               values_to = "Expression") %>%
  # Rename w_suvh456 to "h2a.w suvh456"
  mutate(Sample = ifelse(Sample == "w_suvh456", "h2a.w suvh456", Sample)) %>%
  # Set factor order for samples
  mutate(Sample = factor(Sample, levels = c("Col_0", "cmt23", "atxr56", "suvh456", "h2a.w suvh456", "met1", "ddm1"))) %>%
  # Optional: log transform if needed for better visualization
  mutate(log_Expression = log2(Expression + 1))  # +1 to handle zeros

# Boxplot comparing samples within each cluster (no legend)
p2 <- ggplot(boxplot_data, aes(x = Sample, y = log_Expression, color = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.shape=NA) +
  facet_wrap(~cluster, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("#E68B50", "#CB514A", "#71372E")) +
  labs(title = "TE Expression by Sample Across Clusters",
       x = "",
       y = "Log2 (RPM + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 12))

print(p2)

# ===================== now use published data instead (polyA w2 suvh456, lacks many samples: atxr56, ddm1, met1)

# import the data
w_s456_TPM_polyA <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3_10d_seedlings/w2_h1_suvh456_cmt3_polyA/TPM.tsv"), sep="\t", dec=".", quote=F)
w_s456_TPM_AS_polyA <- as_tibble(read.delim("/groups/berger/user/pierre.bourguet/genomics/RNAseq/w2_h1_suvh456_cmt3_10d_seedlings/w2_h1_suvh456_cmt3_polyA/TPM_AS.tsv"), sep="\t", dec=".", quote=F)

# add a _AS suffix to the Geneid column and append both dataframes
w_s456_TPM_AS_polyA$Geneid <- paste0(w_s456_TPM_AS_polyA$Geneid, "_AS")
w_s456_TPM_polyA_merged <- rbind(w_s456_TPM_polyA, w_s456_TPM_AS_polyA)

boxplot_data_polyA <- w_s456_TPM_polyA_merged %>%
  filter(Geneid %in% merged_cdca7_log2FC_clusters$V4) %>%
  left_join(merged_cdca7_log2FC_clusters %>% select(V4, cluster), by = c("Geneid" = "V4")) %>%
  select(Geneid, cluster, WT, w2, h1, w2_h1, cmt3, w2_cmt3, cmt23, w2_cmt23, suvh456, w2_suvh456) %>%
  pivot_longer(cols = c(WT, w2, h1, w2_h1, cmt3, w2_cmt3, cmt23, w2_cmt23, suvh456, w2_suvh456), 
               names_to = "Sample", 
               values_to = "Expression") %>%
  # Rename w_suvh456 to "h2a.w suvh456"
  mutate(Sample = ifelse(Sample == "w2_suvh456", "w2_suvh456", Sample)) %>%
  # Set factor order for samples
  mutate(Sample = factor(Sample, levels = c("WT", "w2", "h1", "w2_h1", "cmt3", "w2_cmt3", "cmt23", "w2_cmt23", "suvh456", "w2_suvh456"))) %>%
  # Optional: log transform if needed for better visualization
  mutate(log_Expression = log2(Expression + 1))  # +1 to handle zeros

# Boxplot comparing samples within each cluster (no legend)
p2 <- ggplot(boxplot_data_polyA %>% filter(Sample %in% c("WT", "h1", "cmt3", "cmt23", "suvh456")), aes(x = Sample, y = log_Expression, color = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.shape=NA) +
  facet_wrap(~cluster, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("#E68B50", "#CB514A", "#71372E")) +
  labs(title = "TE Expression by Sample Across Clusters",
       x = "",
       y = "Log2 (TPM + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        legend.position = "none")  +
  coord_cartesian(ylim = c(0, 1))

print(p2)

#
# import cmt2 cmt3 drm12 transcriptome data from Stroud et al 2014 NSMB to look at cluster expression in there ####

setwd("/groups/berger/user/pierre.bourguet/genomics/bigwig/2014_NSMB_Stroud_cmt2/tar/mRNA_bigwigs/avg_over_TEs")
files_TE <- c(list.files(pattern = "_rep..tsv")) # list of all target files
Stroud_2014_list <- lapply(X = files_TE, FUN = read.delim, head=F, sep="\t", quote="", comment.char="") # read them
names(Stroud_2014_list) <- files_TE
Stroud_2014_list_data <- lapply(X = Stroud_2014_list, FUN = function(x) { return(x[, c(1,6) ]) } )
# merge together regions covered in all samples.
Stroud_2014_RPM <- as_tibble(Stroud_2014_list_data %>% purrr::reduce(inner_join, by='V1')) # merge all dataframes together
# renaming columns
names(Stroud_2014_RPM) <- c("Geneid", gsub(x = names(Stroud_2014_list_data), pattern = ".tsv", replacement=""))




#
### NOT UPDATED from there ####

grDevices::cairo_pdf("figures/CDCA7-ab_log2_RPM_at_TEs.pdf", width = 4, height = 2)
set.seed(130)

# heatmap of ddm1 upregulated TEs ####

## using VST at ddm1 upTEs

# select samples of interest
m <- vsd_avg %>%
  filter(Geneid %in% ddm1_G2_up_TEs) %>%
  select(any_of(c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2"))) %>%
  as.matrix(.)
# see data spread to scale the heatmap
summary(m)
col_VST <- colorRamp2(seq(from=3, to=12, by=1), scico(n=10, direction=-1, palette="lajolla"))
Heatmap(m,
        name = "VST",
        col = col_VST,
        cluster_rows = T,
        cluster_columns = FALSE,
        width=15, height=15,
        column_title=paste0("ddm1-2 upTEs\nn=", length(ddm1_G2_up_TEs))
)

# sort the heatmap by row sums (TEs with highest expression across genotypes at the top)
m <- m[order(rowSums(m), decreasing = T),] 
ddm1_up_TEs_heatmap <- Heatmap(t(m),
                               name = "VST",
                               col = col_VST,
                               cluster_rows = F,
                               cluster_columns = F,
                               width=15, height=15,
                               column_title=NULL,
                               row_names_side = "left",
                               row_names_gp = grid::gpar(fontsize = 6),
                               use_raster = F,
                               border= T,
                               heatmap_legend_param = list(
                                 title = "VST", at = c(0, 5, 10), 
                                 labels = c("0", "5", "10"),
                                 legend_height = unit(1, "cm"),
                                 legend_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 6),
                                 title_gp = gpar(fontsize = 6)
                               )
) ; ddm1_up_TEs_heatmap

## Z-score heatmap, scaling on VST
z <- VST_merged_avg %>%
  select(-matches("ddm1|F2|mom1|long"))
col_order <- colnames(m)[c(1:3,6,7,4,5)]

m <- as.matrix(z[z$Geneid %in% cdca7_ab_up_TEs,-1])
m <- t(scale(t(m)))
summary(m)
col_zscore <- colorRamp2(c(-2,0,2), c("#2166AC", "#F7F7F7", "#B2182B")) # an alternative from https://personal.sron.nl/~pault/
# Create the heatmap
Heatmap(m,
        name = "z-score",
        col = col_zscore,
        column_order = col_order, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        width=15, height=15
)

# Save the heatmap to a PDF file
ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
file_name <- paste0(output_dir, n, "/heatmap_", y, ".pdf")
pdf(file_name, width = 15, height = 7)
draw(ht, heatmap_legend_side = "right")
dev.off()

# heatmap at ddm1 up TEs
ddm1_G5_up_TEs <- up(all_res$res_ddm1_2_G5, "TEs")
# with all samples
DEG_heatmap(ddm1_G5_up_TEs, "ddm1_G5", z = VST_merged, n = normalization)

#

# make linear correlations between samples using VST at cdca7 up TEs ####

### between cdca7-ab mutants
cdca7_ab_up_PCGs <- intersect(up(all_res$res_ab_1, "features"), up(all_res$res_ab_2, "features"))
ddm1_up_PCGs <- up(all_res$res_ddm1_2_G2, "features")

# linear correlation between samples
cor(vsd_avg_cdca7_up_TEs %>% dplyr::select(-Geneid), method = "pearson")
# scatterplot of cdca7-ab mutants
ggplot(cdca7_ab_vsd, aes(x = ab_1, y = a_long_b)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  geom_abline() +
  labs(title = "cdca7-ab mutants", x = "VST(ab_1)", y = "VST(ab_2)") +
  coord_fixed(ylim=c(2,11), xlim=c(2,11))

data <- cdca7_ab_vsd %>% dplyr::select(c("Col_0", "ab_1", "a_long_b", "ab_2", "ddm1_2_G2", "ddm1_2_G5"))
data <- cdca7_ab_vsd %>% dplyr::select(c("Col_0", "a_1", "a_2", "a_long_1", "a_long_2", "a_long_3"))
pairs(data)
library(GGally)
ggpairs(data, title = "Scatter Plot Matrix for mtcars Dataset", axisLabels = "show")
#coord_fixed(ylim=c(1,11), xlim=c(1,11)) +
#geom_abline(intercept = 0, slope =1)

#
# plot heatmaps ####
suppressPackageStartupMessages(library("RColorBrewer")) ; suppressPackageStartupMessages(library("pheatmap"))
DEG_heatmap <- function(x, y, z, n) {
  if (is.vector(x)==T) {
    if (length(x) > 3) {
      if (n=="rlog") {
        sampleDistMatrix <- as.matrix(subset(z, subset=row.names(z) %in% x))
        title <- paste0(y, "\nn=", length(x),"\n", n)
      }
      else {
        sampleDistMatrix <- as.matrix(log2(subset(z, subset=row.names(z) %in% x)+1))
        title <- paste0(y, "\nn=", length(x),"\nlog2(", n, "+1)")
      }
      rownames(sampleDistMatrix) <- NULL
      colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
      ifelse(!dir.exists(paste0(output_dir, n)), dir.create(paste0(output_dir, n)), FALSE)
      pheatmap(sampleDistMatrix, col=colors, filename=paste0(output_dir, n, "/heatmap_", y, ".pdf"), main=title, cluster_cols = F)
    }
  }
}

mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_summary_norm_mean[,9:(9+length(conditions))], n = normalization) )
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_summary_norm[,9:(8+nrow(samples))], n = normalization) )


# write file with RPM / RPM values for all samples, averaged over replicates ####
write.table(x=cts_summary_norm_mean, file=paste0(args[1], normalization, strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)

# function that adds annotations to a df using its row.names as Geneid
annotate_df <- function(df) { 
  df$Geneid <- row.names(df)
  df <- merge(x = cts_summary[, which(!names(cts_summary) %in% names(cts_summary)[sample_columns])], y = df, by="Geneid")
  row.names(df) <- df$Geneid
  return(df)
}
# median of ratios (from DESeq2, median of ratios to geometric mean). Write tables and heatmaps ####
cts_MoR <- annotate_df( as.data.frame(counts(dds, normalized=T)) )
cts_MoR_mean <- average_replicates(cts_MoR)
write.table(x=cts_MoR, file=paste0(args[1], "ESF", strand, ".tsv"), quote = F, sep="\t", row.names=F, col.names=T)
write.table(x=cts_MoR_mean, file=paste0(args[1], "ESF", strand, "_mean.tsv"), quote = F, sep="\t", row.names=F, col.names=T)
# write heatmaps
mapply(FUN = DEG_heatmap, x=DEGs, y=names(DEGs), MoreArgs = list(z=cts_MoR[,9:(8+nrow(samples))], n = "MoR") )
mapply(FUN = DEG_heatmap, x=DEGs, y=paste0(names(DEGs), "_mean"), MoreArgs = list(z=cts_MoR_mean[,9:(9+length(conditions))], n = "MoR") )
# exploring differences between median of ratios normalization by DESeq and RPM / RPM ####
if (FALSE==TRUE) { # just protecting this code so it's not executed when i run the script
  test <- cts_summary[,sample_columns]
  par(pty = "s")
  plot(x=colSums(test) / colSums(test)[6], y=sizeFactors(dds), xlim=c(0.5,1.7), ylim=c(0.5,1.7), xlab=normalization)
  abline(a = 0, b=1)
  text(x=colSums(test) / colSums(test)[6], sizeFactors(dds), labels=names(sizeFactors(dds)), cex= 0.5, pos=1)
}


# quick single gene check ####

GOI <- ESF %>%
  filter(Geneid == "AT4G37110") %>%
  pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>%
  mutate(condition = gsub("_R.$", "", sample)) %>%
  filter(condition %in% c("a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "Col_0")) %>%
  # average values across replicates
  group_by(Geneid, condition) %>%
  summarise(mean = mean(value)) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2"))) %>%
  arrange(condition)

ggplot(GOI, aes(x=condition, y=mean, color=condition)) +
  geom_quasirandom() +
  theme_classic() +
  coord_flip()

#
#### plot single gene & perform Tukey's HSD #########

plot_gene_expression <- function(data, 
                                 gene_id, 
                                 genotypes, 
                                 colors = NULL, 
                                 genotype_labels = NULL,
                                 alpha = 0.05,
                                 reverse_order = TRUE) {
  
  # Load required libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(agricolae)
  
  # Check if gene exists in the data
  if (!gene_id %in% data$Geneid) {
    available_genes <- head(unique(data$Geneid), 10)  # Show first 10 genes as examples
    stop(paste0("Gene '", gene_id, "' not found in the data.\n",
                "Available genes include: ", paste(available_genes, collapse = ", "), 
                ifelse(length(unique(data$Geneid)) > 10, "...", ""),
                "\nTotal genes in dataset: ", length(unique(data$Geneid))))
  }
  
  # Helper function for Tukey test
  perform_tukey_test <- function(data, value_col, group_col, alpha = 0.05) {
    formula_str <- paste(value_col, "~", group_col)
    aov_result <- aov(as.formula(formula_str), data = data)
    tukey_result <- TukeyHSD(aov_result, conf.level = 1 - alpha)
    
    tukey_df <- as.data.frame(tukey_result[[group_col]])
    tukey_df$comparison <- rownames(tukey_df)
    rownames(tukey_df) <- NULL
    tukey_df <- tukey_df[, c("comparison", "diff", "lwr", "upr", "p adj")]
    tukey_df$significant <- tukey_df$`p adj` < alpha
    
    return(list(
      anova = summary(aov_result),
      tukey = tukey_df,
      tukey_object = tukey_result
    ))
  }
  
  # Helper function for statistical groups
  create_stat_groups <- function(data, value_col, group_col, alpha = 0.05) {
    formula_str <- paste(value_col, "~", group_col)
    aov_result <- aov(as.formula(formula_str), data = data)
    hsd_result <- HSD.test(aov_result, group_col, alpha = alpha)
    
    groups_df <- data.frame(
      condition = rownames(hsd_result$groups),
      mean = hsd_result$groups[, 1],
      groups = hsd_result$groups$groups,
      stringsAsFactors = FALSE
    )
    
    return(list(
      groups = groups_df,
      hsd_result = hsd_result
    ))
  }
  
  # Process the data
  GOI <- data %>%
    filter(Geneid == gene_id) %>%
    pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>%
    mutate(condition = gsub("_R.$", "", sample)) %>%
    filter(condition %in% genotypes) %>%
    mutate(condition = factor(condition, levels = if(reverse_order) rev(genotypes) else genotypes))
  
  # Check if any data remains after filtering
  if (nrow(GOI) == 0) {
    available_conditions <- data %>%
      filter(Geneid == gene_id) %>%
      pivot_longer(-Geneid, names_to = "sample", values_to = "value") %>%
      mutate(condition = gsub("_R.$", "", sample)) %>%
      pull(condition) %>%
      unique()
    
    stop(paste0("No data found for gene '", gene_id, "' with the specified genotypes.\n",
                "Requested genotypes: ", paste(genotypes, collapse = ", "), "\n",
                "Available conditions for this gene: ", paste(available_conditions, collapse = ", ")))
  }
  
  # Calculate means
  GOI_means <- GOI %>%
    group_by(Geneid, condition) %>%
    summarise(mean = mean(value), .groups = "drop")
  
  # Perform statistical analysis
  tukey_results <- perform_tukey_test(GOI, "value", "condition", alpha)
  stat_groups <- create_stat_groups(GOI, "value", "condition", alpha)
  
  # Merge with statistical groups
  GOI_means_with_groups <- merge(GOI_means, stat_groups$groups, 
                                 by.x = "condition", by.y = "condition") %>%  
    select(condition, Geneid, mean = mean.x, groups)
  
  # Apply custom labels if provided
  if (!is.null(genotype_labels)) {
    if (length(genotype_labels) != length(genotypes)) {
      stop("genotype_labels must be the same length as genotypes")
    }
    
    # Create mapping
    label_mapping <- setNames(genotype_labels, genotypes)
    
    # Apply to data
    GOI$condition_label <- factor(label_mapping[as.character(GOI$condition)], 
                                  levels = if(reverse_order) rev(genotype_labels) else genotype_labels)
    GOI_means_with_groups$condition_label <- factor(label_mapping[as.character(GOI_means_with_groups$condition)], 
                                                    levels = if(reverse_order) rev(genotype_labels) else genotype_labels)
    
    x_var <- "condition_label"
  } else {
    x_var <- "condition"
  }
  
  # Create the plot
  p <- ggplot() +
    geom_col(data = GOI_means_with_groups, 
             aes_string(x = x_var, y = "mean", fill = x_var), 
             alpha = 0.7) +
    geom_quasirandom(data = GOI, 
                     aes_string(x = x_var, y = "value", color = x_var), 
                     width = 0.2, size = 2, alpha = 1) +
    geom_text(data = GOI_means_with_groups, 
              aes_string(x = x_var, y = "mean + max(mean) * 0.05", label = "groups"),
              size = 4) +
    theme_classic() +
    coord_flip() +
    labs(x = "", y = paste0(gene_id, " normalized expression"))
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    if (length(colors) != length(genotypes)) {
      stop("colors must be the same length as genotypes")
    }
    
    color_values <- if(reverse_order) rev(colors) else colors
    p <- p + 
      scale_fill_manual(values = color_values) +
      scale_color_manual(values = color_values)
  }
  
  # Return list with plot and statistical results
  return(list(
    plot = p,
    tukey_results = tukey_results,
    stat_groups = stat_groups$groups,
    data = GOI,
    means = GOI_means_with_groups
  ))
}

#  CDCA7a with custom colors
result <- plot_gene_expression(ESF, "AT4G37110", 
                               c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2", "ddm1_2_G5"),
                               colors = c(col_muted[10], col_muted_2_replicates[c(1:4,7:8,5,6)]))

# Display the plot
print(result$plot)
# Access statistical results
print(result$stat_groups)

# save the plot
svglite::svglite("figures/CDCA7a_expression.svg", width = 50*mm_to_inches, height = 50*mm_to_inches)
result$plot + theme_horizontal_nature  +  theme(legend.position = "none") 
dev.off()

# CDCA7b: error because not detected (as max 1 / 2 reads in only a few samples, was filtered out)
result <- plot_gene_expression(ESF, "AT2G23530", 
                               c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2", "ddm1_2_G5"),
                               colors = c(col_muted[10], col_muted_2_replicates[c(1:4,7:8,5,6)]))

# DDM1: also not detected
result <- plot_gene_expression(ESF, "AT5G66750", 
                               c("Col_0", "a_1", "a_2", "b_1", "b_2", "ab_1", "ab_2", "ddm1_2_G2", "ddm1_2_G5"),
                               colors = c(col_muted[10], col_muted_2_replicates[c(1:4,7:8,5,6)]))
