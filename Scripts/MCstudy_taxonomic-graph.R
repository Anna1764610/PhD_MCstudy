### Description: Statistics and Graphing 16S dataset        ###
### Author: Anna Li, date: 12012026                         ###
### Email: anna.li@adelaide.edu.au                          ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(tidyr)

#additional
library(scales)
library(BiomeHue)
library(legendry)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====Load filtered counts, taxonomy, metadata, phyloseq objects====
physeq_raw <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw.rds"))
physeq_prop <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_prop.rds"))
physeq_raw_tree <- readRDS(file.path(project_dir, "phyloseq_objects", "physeq_raw_tree.rds"))
map_data <-readRDS(file.path(project_dir, "phyloseq_objects", "map_data.rds"))
#plus PC data if necessary
PC_data.df <-readRDS(file.path(project_dir, "phyloseq_objects", "PC_data.df.rds"))


#define timepoints
tp_order <- c("d0", "d4", "d7", "d14", "d21", "d28", "d35")
#fixing physeq_prop phyloseq object to have ordered timepoints
sample_data(physeq_prop)$Timepoint <- factor(
  sample_data(physeq_prop)$Timepoint,
  levels = tp_order,
  ordered = TRUE
)

#fixing physeq_raw phyloseq object to have ordered timepoints
sample_data(physeq_raw)$Timepoint <- factor(
  sample_data(physeq_raw)$Timepoint,
  levels = tp_order,
  ordered = TRUE
)

##====Creating taxonomic bar plots for donor samples====

#Create subset of physeq object physeq_prop with only the donor samples- taxonomy, asv and metadata included
donor_samples <- c("25-06236","25-06237","25-06294","25-06295")
physeq_donor <- prune_samples(donor_samples, physeq_prop)
#check 
sample_names(physeq_donor)

#subset by study# option - have not run yet!!!
physeq_MBC1 <- subset_samples(physeq_prop, Study_ID == "MBC1")

#plot
bar_plot_donor_genus <- plot_bar(physeq_donor, fill="Genus") +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=8)) +
  labs(title = "Taxonomic composition at genus for donors")
bar_plot_donor_genus #print

#agglomerate to family - more readable
physeq_family_donor <- tax_glom(physeq_donor, taxrank = "Family")

bar_plot_donor_family <- plot_bar(physeq_family_donor, fill="Family") +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=8)) +
  labs(title = "Taxonomic composition at family for donors")
bar_plot_donor_family #print

##====Creating taxonomic bar graphs for timepoints====
physeq_day0 <- subset_samples(physeq_prop, Timepoint == "d0")
physeq_day4 <- subset_samples(physeq_prop, Timepoint == "d4")
physeq_day7 <- subset_samples(physeq_prop, Timepoint == "d7")
physeq_day14 <- subset_samples(physeq_prop, Timepoint == "d14")
physeq_day21 <- subset_samples(physeq_prop, Timepoint == "d21")
physeq_day28 <- subset_samples(physeq_prop, Timepoint == "d28")
physeq_day35 <- subset_samples(physeq_prop, Timepoint == "d35")
#may enter in additional definitions^ for only recipients group = "allo" or "syn"
#check
sample_names(physeq_day4)

#plot at genus level - highest specificity
bar_plot_day28_genus <- plot_bar(physeq_day28, fill="Genus") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=60, hjust=1, size=8),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(title = "Taxonomic composition at genus level for Day 28")

bar_plot_day28_genus #print

#agglomerate to family - more readable
physeq_family_day28 <- tax_glom(physeq_day28, taxrank = "Family")
#plot at family level
bar_plot_day28_family <- plot_bar(physeq_family_day28, fill = "Family") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=60, hjust=1, size=8),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(title = "Taxonomic composition at family for day 28")

bar_plot_day28_family #print

#====Defining unique genera====
#MBC1----
MBC1_physeq_donors <- subset_samples(physeq_raw, SampleID %in% c("25-06236", "25-06237")) #use physeq_raw 
MBC1_physeq_donors_genus <- tax_glom(MBC1_physeq_donors, taxrank = "Genus") #agglomerate
#Checking raw count table
otu_genus_donor_MBC1 <- as.data.frame(otu_table(MBC1_physeq_donors_genus)) #raw count table
#Create taxonomy dataframe for reference
tax_df <- as.data.frame(tax_table(MBC1_physeq_donors_genus))
otu_genus_donor_MBC1$Genus <-tax_df$Genus #add genus identity to raw count table
#Check file otu_genus_donor_MBC1 to see differences

#MBC2----
MBC2_physeq_donors <- subset_samples(physeq_raw, SampleID %in% c("25-06294", "25-06295")) 
MBC2_physeq_donors_genus <- tax_glom(MBC2_physeq_donors, taxrank = "Genus") #agglomerate
#Checking raw count table
otu_genus_donor_MBC2 <- as.data.frame(otu_table(MBC2_physeq_donors_genus)) #raw count table
#Create taxonomy dataframe for reference
tax_df <- as.data.frame(tax_table(MBC2_physeq_donors_genus))
otu_genus_donor_MBC2$Genus <-tax_df$Genus #add genus identity to raw count table
#Check file otu_genus_donor_MBC2 to see differences

#Decision point----
#Identified for MBC1 - akkermansia is the genus present in allo donor but not in syn donor - but not consistent
#Identified for both MBC1 and MBC2 - Atopobiaceae as the genus present in allo donor but not in syn donor
#Identified for both MBC1 and MBC2 - Desulfovibrionaceae and Ruminococcus as present in syn donor but not in allo donor
#Is this true for recipients? not sure yet 

#====Plotting abundance of specific genus over time====
#Agglomerate
physeq_genus <- tax_glom(physeq_raw, taxrank = "Genus")

#Allo-donor specific genus: Atopobiaceae----
physeq_atop <- subset_taxa(physeq_genus, Genus == "Atopobiaceae") #filter to genus of choice
physeq_atop #checking

#convert to long dataframe, can check to see that at d0 abundance is low in recipients and syn donor
df_atop <- psmelt(physeq_atop)

#for this df file specifically --
df_atop$Timepoint <- factor(df_atop$Timepoint, levels = tp_order, ordered = TRUE)
levels(df_atop$Timepoint) #this checks if timepoint is ordered

#define means...ignore for now --
#df_atop_mean <- df_atop %>%
  group_by(Group, Timepoint) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

#plot x and y graph
atop_plot <- ggplot(df_atop, aes(x = Timepoint, y = Abundance, colour = Group)) +
  geom_point(aes(size = Group, shape = Study_ID), alpha = 0.8) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.5) +
  scale_colour_manual(values = group_cols) +   # your custom colours
  scale_size_manual(
    values = c(
      "allo" = 3,
      "syn" = 3,
      "allo_donor" = 6,   # bigger
      "syn_donor" = 6     # bigger
    )) +
  labs(
    title = "Atopobiaceae abundance over time",
    x = "Timepoint",
    y = "Raw counts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  )

atop_plot

#Syn-donor specific genus: Desulfovibrionaceae----
physeq_desulf <- subset_taxa(physeq_genus, Genus == "Desulfovibrionaceae") #filter to genus of choice
physeq_desulf #checking
df_desulf <- psmelt(physeq_desulf) #convert to long dataframe

#plot x and y graph
desulf_plot <- ggplot(df_desulf, aes(x = Timepoint, y = Abundance, colour = Group)) +
  geom_point(aes(size = Group, shape = Study_ID), alpha = 0.8) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.5) +
  scale_colour_manual(values = group_cols) +   # your custom colours
  scale_size_manual(
    values = c(
      "allo" = 3,
      "syn" = 3,
      "allo_donor" = 6,   # bigger
      "syn_donor" = 6     # bigger
    )) +
  labs(
    title = "Desulfovibrionaceae abundance over time",
    x = "Timepoint",
    y = "Raw counts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  )

desulf_plot

#Syn-donor specific genus: 	Ruminococcus----
physeq_rumi <- subset_taxa(physeq_genus, Genus == "Ruminococcus") #filter to genus of choice
physeq_rumi #checking
df_rumi <- psmelt(physeq_rumi) #convert to long dataframe

#plot x and y graph
rumi_plot <- ggplot(df_rumi, aes(x = Timepoint, y = Abundance, colour = Group)) +
  geom_point(aes(size = Group, shape = Study_ID), alpha = 0.8) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.5) +
  scale_colour_manual(values = group_cols) +   # your custom colours
  scale_size_manual(
    values = c(
      "allo" = 3,
      "syn" = 3,
      "allo_donor" = 6,   # bigger
      "syn_donor" = 6     # bigger
    )) +
  labs(
    title = "Ruminococcus abundance over time",
    x = "Timepoint",
    y = "Raw counts"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  )

rumi_plot

#====Saving graphs====
#Create consistent file path
save_plot <- function(plot_obj, filename, width = 16, height = 8) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

#Saving graphs
save_plot(bar_plot_donor_genus, "bar_plot_donor_genus.pdf")
save_plot(bar_plot_donor_family, "bar_plot_donor_family.pdf")

save_plot(bar_plot_day0_genus, "bar_plot_day0_genus.pdf")
save_plot(bar_plot_day0_family, "bar_plot_day0_family.pdf")

save_plot(bar_plot_day4_genus, "bar_plot_day4_genus.pdf")
save_plot(bar_plot_day4_family, "bar_plot_day4_family.pdf")

save_plot(bar_plot_day7_genus, "bar_plot_day7_genus.pdf")
save_plot(bar_plot_day7_family, "bar_plot_day7_family.pdf")

save_plot(bar_plot_day28_genus, "bar_plot_day28_genus.pdf")
save_plot(bar_plot_day28_family, "bar_plot_day28_family.pdf")

save_plot(atop_plot, "atop_plot.pdf")
save_plot(desulf_plot, "desulf_plot.pdf")
save_plot(rumi_plot, "rumi_plot.pdf")
