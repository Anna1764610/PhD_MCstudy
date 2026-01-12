### Description: Statistics and Graphing 16S dataset        ###
### Author: Anna Li, date: 08012026                         ###
### Email: anna.li@adelaide.edu.au                          ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====Load phyloseq objects====
physeq_raw <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw.rds"))
physeq_prop <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_prop.rds"))
physeq_raw_tree <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw_tree.rds"))
map_data <-readRDS(file.path(project_dir, "phyloseq_objects", "map_data.rds"))

#====Calculate alpha diversity====
alpha_div <- estimate_richness(physeq_raw, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_div)
library(dplyr) #for data manipulation/component of tidyverse -  for cleanup of alpha_div
rownames(alpha_div) <- sub("^X", "", rownames(alpha_div))
rownames(alpha_div) <- gsub("\\.", "-", rownames(alpha_div))  # replace dots with dashes
#=====Graphing alpha diversity====
#Create object without donor samples for alpha diversity plots
Recipient_data.df <- PC_data.df[!rownames(PC_data.df)%in% c("25-06236","25-06294","25-06237","25-06295"),]

#Graphing shannon diversity across timepoints
Shannon_timeline_plot <- ggplot(Recipient_data.df, aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Mouse_ID), alpha = 0.2, linewidth = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity per Timepoint",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme(legend.position = "right")

Shannon_timeline_plot

#====Calculate Beta diversity====
#----Bray-curtis----
## All samples
braycurtis_dist = phyloseq::distance(physeq_prop,method="bray") #input physeq_prop for bray curtis analysis as it calculates using composition
braycurtis_pc = pcoa(braycurtis_dist)

## Subset by study ID (distances are computed within the study run only)
physeq_MBC1 <- subset_samples(physeq_prop, Study_ID == "MBC1")
physeq_MBC2 <- subset_samples(physeq_prop, Study_ID == "MBC2")

map_data_MBC1 <-sample_data(physeq_MBC1)
map_data_MBC2 <-sample_data(physeq_MBC2)

## MBC1
bray_MBC1 <- phyloseq::distance(physeq_MBC1, method="bray")
pc_MBC1 <- pcoa(bray_MBC1)

## MBC2
bray_MBC2 <- phyloseq::distance(physeq_MBC2, method="bray")
pc_MBC2 <- pcoa(bray_MBC2)

#-----unifrac----
#~~~~~~
##OPTION: root the tree or export from QIIME...
library(phytools)
# extract tree
tree <- phy_tree(physeq_raw_tree)
# midpoint root
tree_rooted <- midpoint.root(tree)
# put the rooted tree back into phyloseq
phy_tree(physeq_raw_tree) <- tree_rooted
#~~~~~~

#unifrac
unifrac_dist = phyloseq::distance(physeq_raw_tree,method="unifrac")
unifrac_pc = pcoa(unifrac_dist)

#weighted Unifrac
wunifrac_dist = phyloseq::distance(physeq_raw_tree,method="wunifrac")
wunifrac_pc = pcoa(wunifrac_dist)

#====Creating dataframes with metadata for graphing====
## All data
PC_data.df = data.frame(map_data,alpha_div,
                        "BC1" = braycurtis_pc$vectors[,1],
                        "BC2" = braycurtis_pc$vectors[,2],
                        "UF1" = unifrac_pc$vectors[,1],
                        "UF2" = unifrac_pc$vectors[,2],
                        "WU1" = wunifrac_pc$vectors[,1],
                        "WU2" = wunifrac_pc$vectors[,2])
group_cols <- c(
  "allo"                = "#1b9e77",   # green
  "syn"                 = "#d95f02",   # orange
  "allo_donor"          = "#7570b3",   # purple
  "syn_donor"           = "#e7298a"    # pink
)

# Add an ordered factor column 
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint,
                                       #levels = timepoint_levels,
                                       ordered = TRUE)

#Define allo_donor and syn_donor as 'timepoint_ordered' for graphing by time
PC_data.df$Timepoint_ordered <- as.character(PC_data.df$Timepoint_ordered)
PC_data.df$Timepoint_ordered[PC_data.df$SampleID %in% c("25-06236", "25-06294")] <- "allo_donor"
PC_data.df$Timepoint_ordered[PC_data.df$SampleID %in% c("25-06237", "25-06295")] <- "syn_donor"
PC_data.df$Timepoint_ordered <- factor(PC_data.df$Timepoint_ordered, 
                                       levels = c("syn_donor","d0","d4","d7","d14","d21","d28","d35","allo_donor"))

#----Repeat for subsets---- 
##Create dataframe for MBC1
PC_data_MBC1 <- data.frame(map_data_MBC1,
                           "BC1" = pc_MBC1$vectors[,1],
                           "BC2" = pc_MBC1$vectors[,2])
group_cols <- c(
  "allo"                = "#1b9e77",   # green
  "syn"                 = "#d95f02",   # orange
  "allo_donor"          = "#7570b3",   # purple
  "syn_donor"           = "#e7298a"    # pink
)

PC_data_MBC1$Timepoint_ordered <- factor(PC_data_MBC1$Timepoint,
                                         #levels = timepoint_levels,
                                         ordered = TRUE)
PC_data_MBC1$Timepoint_ordered <- as.character(PC_data_MBC1$Timepoint_ordered)
PC_data_MBC1$Timepoint_ordered[PC_data_MBC1$SampleID %in% c("25-06236")] <- "allo_donor"
PC_data_MBC1$Timepoint_ordered[PC_data_MBC1$SampleID %in% c("25-06237")] <- "syn_donor"
PC_data_MBC1$Timepoint_ordered <- factor(PC_data_MBC1$Timepoint_ordered, 
                                         levels = c("syn_donor","d0","d4","d7","d14","d21","d28","d35","allo_donor"))

##Create dataframe for MBC2
PC_data_MBC2 <- data.frame(map_data_MBC2,
                           "BC1" = pc_MBC2$vectors[,1],
                           "BC2" = pc_MBC2$vectors[,2])
group_cols <- c(
  "allo"                = "#1b9e77",   # green
  "syn"                 = "#d95f02",   # orange
  "allo_donor"          = "#7570b3",   # purple
  "syn_donor"           = "#e7298a"    # pink
)

PC_data_MBC2$Timepoint_ordered <- factor(PC_data_MBC2$Timepoint,
                                         #levels = timepoint_levels,
                                         ordered = TRUE)
PC_data_MBC2$Timepoint_ordered <- as.character(PC_data_MBC2$Timepoint_ordered)
PC_data_MBC2$Timepoint_ordered[PC_data_MBC2$SampleID %in% c("25-06294")] <- "allo_donor"
PC_data_MBC2$Timepoint_ordered[PC_data_MBC2$SampleID %in% c("25-06295")] <- "syn_donor"
PC_data_MBC2$Timepoint_ordered <- factor(PC_data_MBC2$Timepoint_ordered, 
                                         levels = c("syn_donor","d0","d4","d7","d14","d21","d28","d35","allo_donor"))

##=====Graphing beta diversity====
#Bray curtis
#plot of all samples coloured by timepoint
BC_plot_original <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                            colour = Study_ID, #Change to Timepoint_ordered or Study_ID
                                            size = Group)) + #Inside aes() you don’t set exact colors/sizes; instead, you tell ggplot which variable controls them
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(
    "allo"           = 2,
    "syn"            = 2,
    "allo_donor"     = 6,
    "syn_donor"      = 6
  ), guide = "none") + 
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA (coloured by Study_ID)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    colour = "Study_ID"
  )

BC_plot_original

#Identify donor samples
donor_samples <- c("25-06236","25-06237","25-06294","25-06295")

# Keep only real timepoints (for plotting individual graphs)
timepoints <- setdiff(unique(PC_data.df$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))

timepoints <- setdiff(levels(PC_data.df$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))

BC_PLOTS <- list()

for(tp in timepoints){
  df_tp <- PC_data.df %>% filter(Timepoint_ordered == tp)
  
  BC_PLOTS[[tp]] <- ggplot(df_tp, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
    geom_point(alpha = 0.8, stroke = 0.5) +
    # ellipses for each group at that timepoint
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay donor samples: keep group color, but larger + black outline
    geom_point(data = PC_data.df %>% filter(SampleID %in% donor_samples),
               aes(x = BC1, y = BC2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "allo"           = 2,
      "syn"            = 2,
      "allo_donor"     = 6,
      "syn_donor"      = 6
    ), guide = "none") + 
    theme_minimal(base_size = 14) +
    labs(title = paste("Bray–Curtis PCoA:", tp),
         x = "PCoA Axis 1",
         y = "PCoA Axis 2")
}
BC_PLOTS

library(cowplot)
BC_PLOTS_all <- plot_grid(plotlist = BC_PLOTS[timepoints],
                           ncol = 2, labels = "AUTO")
BC_PLOTS_all

#By study ID
#subset by study ID
donor_samples_MBC1 <- c("25-06236","25-06237")
donor_samples_MBC2 <- c("25-06294","25-06295")

# Keep only real timepoints (for plotting individual graphs)
#MBC1
timepoints <- setdiff(unique(PC_data_MBC1$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))

timepoints <- setdiff(levels(PC_data_MBC1$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))
#MBC2
timepoints <- setdiff(unique(PC_data_MBC2$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))

timepoints <- setdiff(levels(PC_data_MBC2$Timepoint_ordered),
                      c("allo_donor", "syn_donor"))


#graphs
BC_PLOTS_MBC1 <- list()
BC_PLOTS_MBC2 <- list()

#MBC1
for(tp in timepoints){
  df_tp <- PC_data_MBC1 %>% filter(Timepoint_ordered == tp)
  
  BC_PLOTS_MBC1[[tp]] <- ggplot(df_tp, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
    geom_point(alpha = 0.8, stroke = 0.5) +
    # ellipses for each group at that timepoint
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay donor samples: keep group color, but larger + black outline
    geom_point(data = PC_data_MBC1 %>% filter(SampleID %in% donor_samples_MBC1),
               aes(x = BC1, y = BC2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "allo"           = 2,
      "syn"            = 2,
      "allo_donor"     = 6,
      "syn_donor"      = 6
    ), guide = "none") + 
    theme_minimal(base_size = 14) +
    labs(title = paste("Bray–Curtis PCoA:", tp),
         x = "PCoA Axis 1",
         y = "PCoA Axis 2")
}
BC_PLOTS_MBC1

library(cowplot)
BC_PLOTS_MBC1_all <- plot_grid(plotlist = BC_PLOTS_MBC1[timepoints],
                           ncol = 2, labels = "AUTO")
BC_PLOTS_MBC1_all


#MBC2
for(tp in timepoints){
  df_tp <- PC_data_MBC2 %>% filter(Timepoint_ordered == tp)
  
  BC_PLOTS_MBC2[[tp]] <- ggplot(df_tp, aes(x = BC1, y = BC2, colour = Group, size = Group)) +
    geom_point(alpha = 0.8, stroke = 0.5) +
    # ellipses for each group at that timepoint
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay donor samples: keep group color, but larger + black outline
    geom_point(data = PC_data_MBC2 %>% filter(SampleID %in% donor_samples_MBC2),
               aes(x = BC1, y = BC2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "allo"           = 2,
      "syn"            = 2,
      "allo_donor"     = 6,
      "syn_donor"      = 6
    ), guide = "none") + 
    theme_minimal(base_size = 14) +
    labs(title = paste("Bray–Curtis PCoA:", tp),
         x = "PCoA Axis 1",
         y = "PCoA Axis 2")
}
BC_PLOTS_MBC2

library(cowplot)
BC_PLOTS_MBC2_all <- plot_grid(plotlist = BC_PLOTS_MBC2[timepoints],
                               ncol = 2, labels = "AUTO")
BC_PLOTS_MBC2_all

#====save graphs for export====
#Create consistent filepath and parameters for saving graphs

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

save_plot(Shannon_timeline_plot, "Shannon_timeline_plot.pdf")
save_plot(BC_plot_original, "BC_plot_all_study-ID.pdf")
save_plot(BC_PLOTS_all,"BrayCurtis_PCoA_allTimepoints.pdf", width = 12, height = 10)
save_plot(BC_PLOTS_MBC1_all, "BrayCurtis_PCoA_MBC1.pdf", width = 12, height = 10)
save_plot(BC_PLOTS_MBC2_all, "BrayCurtis_PCoA_MBC2.pdf", width = 12, height = 10)

#====save the current workspace====
saveRDS(PC_data.df, file=file.path(project_dir, "Phyloseq_objects","PC_data.df.rds"))
saveRDS(PC_data_MBC1, file=file.path(project_dir, "Phyloseq_objects","PC_data.MBC1.rds"))
saveRDS(PC_data_MBC2, file=file.path(project_dir, "Phyloseq_objects","PC_data.MBC2.rds"))
