### Description: Statistics and Graphing 16S dataset        ###
### Author: Anna Li, date: 19042026                         ###
### Email: anna.li@adelaide.edu.au                          ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(tidyr)
library(DESeq2)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====taxonomy graph of donors====
#Create subset of physeq object physeq_prop with only the donor samples- taxonomy, asv and metadata included
donor_samples <- c("25-06236","25-06237","25-06294","25-06295")
physeq_donor <- prune_samples(donor_samples, physeq_prop)
#check 
sample_names(physeq_donor)

## Define at family level
barplot_donor <- tax_glom(physeq_donor, taxrank = "Family")
barplot_donor_family <- psmelt(barplot_donor) #Melt phyloseq table into a data frame for plotting

barplot_donor_family$Family[is.na(barplot_donor_family$Family)] <- "Unclassified"

#list of genus
unique_family <- unique(taxa_filtered$Family)
unique_family

#Assign colours
family_colour = c(
  "Lactobacillaceae" = "red4",
  "Prevotellaceae" = "blue4",
  "Bacteroidaceae" = "forestgreen",
  "Muribaculaceae" = "goldenrod2",
  "Marinifilaceae" = "darkorange2",
  "Desulfovibrionaceae" = "purple4",
  "Rikenellaceae" = "cyan4",
  "Tannerellaceae" = "magenta4",
  "Lachnospiraceae" = "salmon2",
  "Sutterellaceae" = "powderblue",
  "Ruminococcaceae" = "#3C280D",
  "Gastranaerophilales" = "chartreuse4",
  "Erysipelatoclostridiaceae" = "deeppink3",
  "Erysipelotrichaceae" = "orchid4",
  "[Eubacterium]_coprostanoligenes_group" = "sienna3",
  "Acholeplasmataceae" = "turquoise4",
  "Clostridia_UCG-014" = "slateblue3",
  "Oscillospiraceae" = "peachpuff1",
  "Enterobacteriaceae" = "lightskyblue4",
  "Staphylococcaceae" = "firebrick3",
  "Unclassified" = "grey70",
  "Anaerovoracaceae" = "springgreen4",
  "Clostridiaceae" = "tan4",
  "Butyricicoccaceae" = "seagreen3",
  "Eggerthellaceae" = "wheat4",
  "Clostridia_vadinBB60_group" = "darkgoldenrod3",
  "Akkermansiaceae" = "dodgerblue3",
  "Peptococcaceae" = "olivedrab4",
  "Atopobiaceae" = "coral4",
  "UCG-010" = "mediumpurple4",
  "Christensenellaceae" = "aquamarine4",
  "RF39" = "navy",
  "Saccharimonadaceae" = "palevioletred4",
  "Monoglobaceae" = "darkseagreen4",
  "Peptostreptococcaceae" = "indianred3",
  "Defluviitaleaceae" = "khaki4",
  "Enterococcaceae" = "violetred4"
)

##Plot family bar graph
barplot_donor <- ggplot(barplot_donor_family, aes(x = Unique_mouse_ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Donor_ID, scales = "free_x") +
  scale_fill_manual(values = family_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance (Family level)"), x = "sample", y = "Relative Abundance")

barplot_donor

#====Creating taxonomic bar graphs for all timepoints====
#define timepoints
tp_order <- c("d0", "d4", "d7", "d14", "d21", "d28", "d35")
#fixing physeq_prop phyloseq object to have ordered timepoints
sample_data(physeq_prop)$Timepoint <- factor(
  sample_data(physeq_prop)$Timepoint,
  levels = tp_order,
  ordered = TRUE
)

#fixing physeq_raw phyloseq object to have ordered timepoints - just in case
sample_data(physeq_raw)$Timepoint <- factor(
  sample_data(physeq_raw)$Timepoint,
  levels = tp_order,
  ordered = TRUE
)

#subset physeq to timepoints
physeq_day0 <- subset_samples(physeq_prop, Timepoint == "d0")
physeq_day4 <- subset_samples(physeq_prop, Timepoint == "d4")
physeq_day7 <- subset_samples(physeq_prop, Timepoint == "d7")
physeq_day14 <- subset_samples(physeq_prop, Timepoint == "d14")
physeq_day21 <- subset_samples(physeq_prop, Timepoint == "d21")
physeq_day28 <- subset_samples(physeq_prop, Timepoint == "d28")
physeq_day35 <- subset_samples(physeq_prop, Timepoint == "d35")

## Define as family level - adjust for physeq object ***here*** as needed
barplot <- tax_glom(physeq_day21, taxrank = "Family")
barplot_recipient_family <- psmelt(barplot) #Melt phyloseq table into a data frame for plotting
#define unclassified
barplot_recipient_family$Family[is.na(barplot_donor_family$Family)] <- "Unclassified"

##Plot family bar graph
barplot <- ggplot(barplot_recipient_family, aes(x = Unique_mouse_ID, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Donor_ID, scales = "free_x") +
  scale_fill_manual(values = family_colour) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = paste("Relative Abundance by donor ID - Day21"), x = "sample", y = "Relative Abundance")

barplot

#====Differential abundance

#Create subset of physeq object physeq_raw with only the donor samples- taxonomy, asv and metadata included
#physeq_donor <- prune_samples(donor_samples, physeq_raw)
#deseq? or use LINDA somehow? 


#====Identify differential abundant genus in each donor (log fold change)====
#create a dataframe that contains ***genus*** from physeq_prop
genus_ps <- tax_glom(physeq_prop, taxrank = "Genus")
genus_df<-psmelt(genus_ps)
#create ordered time column
genus_df$Timepoint_ordered <- factor(
  genus_df$Timepoint,
  levels = c("d0","d4","d7","d14","d21","d28","d35"),
  ordered = TRUE
)

#create a dataframe that contains **family** from physeq_prop
family_ps <- tax_glom(physeq_prop, taxrank = "Family")
family_df <-psmelt(family_ps)
#create ordered times
family_df$Timepoint_ordered <- factor(
  family_df$Timepoint,
  levels = c("d0","d4","d7","d14","d21","d28","d35"),
  ordered = TRUE
)

#Identify fold change of each genus by bone marrow donor sample
#positive means enriched in allo mouse sample
#negative means enriched in syn mouse sample

#MBC1
MBC1donor_genus_df <-genus_df %>%
  filter(SampleID %in% c("25-06236", "25-06237")) %>%
  select(Genus, SampleID, Abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Abundance)

MBC1donor_genus_df<-MBC1donor_genus_df %>%
  mutate(log2FC = log2((`25-06236` + 1e-6) / (`25-06237` + 1e-6)))

#MBC2
MBC2donor_genus_df <-genus_df %>%
  filter(SampleID %in% c("25-06294", "25-06295")) %>%
  select(Genus, SampleID, Abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Abundance)

MBC2donor_genus_df<-MBC2donor_genus_df %>%
  mutate(log2FC = log2((`25-06294` + 1e-6) / (`25-06295` + 1e-6)))


#====longitudinal abundance of certain taxa over time====

#define taxa of interest ***Plug in taxa and run code from here***
taxa_of_interest <- "Lactobacillaceae"

#CHOOSE:
#1. genus 
genus_target_df <- genus_df %>%
  filter(Genus == taxa_of_interest)

#OR

#2. family
family_target_df <- family_df %>%
  filter(Family == taxa_of_interest)

# 1 and 2. define donor value 

donor_value <-family_target_df %>% #change to genus_target_df or family_target_df
  filter(Sample %in% c("25-06294")) %>% #change donor
  summarise(mean_abund = mean(Abundance)) %>%
  pull(mean_abund)

#1. plot genus trajectory

genus_trajetory <-ggplot(genus_target_df %>%
                           filter(Group %in% c("allo", "syn"), 
                                  Study_ID%in% ("MBC1"),
                                  Timepoint%in% c("d0","d4","d7","d14","d21","d28","d35"),),
                         aes(x = Timepoint_ordered,
                             y = Abundance,
                             group = Unique_mouse_ID,
                             color = Group)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = Group),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  geom_hline(yintercept = donor_value, linetype = "dashed")+
  theme_classic() +
  labs(title = paste("MBC1 trajectory of", taxa_of_interest),
       y = "Relative abundance")

genus_trajetory

#2. plot family trajectory

family_trajetory <-ggplot(family_target_df %>%
                            filter(Group %in% c("allo", "syn"), 
                                   Study_ID%in% ("MBC2"),
                                   Timepoint%in% c("d0","d4","d7","d14","d21","d28","d35"),),
                          aes(x = Timepoint_ordered,
                              y = Abundance,
                              group = Unique_mouse_ID,
                              color = Group)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = Group),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  geom_hline(yintercept = donor_value, linetype = "dashed")+
  theme_classic() +
  labs(title = paste("MBC2 Family level trajectory of", taxa_of_interest),
       y = "Relative abundance")

family_trajetory


#====Saving graphs====
#Create consistent file path
save_plot <- function(plot_obj, filename, width = 12, height = 8) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

#Saving graphs
save_plot(barplot_donor, "barplot_donor.pdf")
save_plot(barplot, "barplot_recipients.pdf")
save_plot(genus_trajetory, "MBC2 trajectory of Eggerthellaceae.pdf", width= 6, height=4)
save_plot(family_trajetory, "MBC2 Family level trajectory of Lactobacillaceae.pdf", width= 6, height=4)
