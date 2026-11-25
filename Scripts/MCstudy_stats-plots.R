### Description: Statistics and Graphing 16S dataset        ###
### Author: Anna Li, date: 08012026 to 19042026             ###
### Email: anna.li@adelaide.edu.au                          ###

## Load Libraries
library(ape)              # Required for read.tree
library(phyloseq)         # Required for convert relative abundance to phyloseq object
library(vegan)            # Required for adonis2 - beta diversity
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(gridExtra)        # Required for arranging plots (grid.arrange)
library(dplyr)
library(pairwiseAdonis)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#====Calculate alpha diversity====
alpha_div <- estimate_richness(physeq_raw, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_div)
library(dplyr) #for data manipulation/component of tidyverse -  for cleanup of alpha_div
rownames(alpha_div) <- sub("^X", "", rownames(alpha_div))
rownames(alpha_div) <- gsub("\\.", "-", rownames(alpha_div))  # replace dots with dashes
#====Calculate Beta diversity====
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

#----Aitchison----
ps_clr <- microbiome::transform(physeq_prop, "clr") #input is proportional
clr_counts <- data.frame(
  otu_table(ps_clr),
  check.names = FALSE
)
aitchison_dist = dist(t(clr_counts),method = "euclidean") #The Aitchison distance is just Euclidean on CLR transformed

# MDS (PCoA on Aitchison distance)
mds <- cmdscale(aitchison_dist, k = 2, eig = TRUE)

mds_points <- mds$points
eigvals <- mds$eig

# % variation explained (Aitchison should be Euclidean -> typically all +ve,
# but we guard just in case)
pos_eig <- eigvals[eigvals > 0]
pct_var <- pos_eig / sum(pos_eig) * 100

#define labels
xlab_mds1 <- paste0("MDS1 (", round(pct_var[1], 1), "%)")
ylab_mds2 <- paste0("MDS2 (", round(pct_var[2], 1), "%)")

#----Unifrac----
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

#====Creating dataframe with metadata for graphing====
## All data
PC_data.df = data.frame(map_data,alpha_div,
                        "BC1" = braycurtis_pc$vectors[,1],
                        "BC2" = braycurtis_pc$vectors[,2],
                        "UF1" = unifrac_pc$vectors[,1],
                        "UF2" = unifrac_pc$vectors[,2],
                        "WU1" = wunifrac_pc$vectors[,1],
                        "WU2" = wunifrac_pc$vectors[,2],
                        "MDS1" = mds_points[,1],
                        "MDS2" = mds_points[,2])
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

#====Graphing alpha diversity====
#Create object without donor samples for alpha diversity plots
Recipient_data.df <- PC_data.df[!rownames(PC_data.df)%in% c("25-06236","25-06294","25-06237","25-06295"),]

Recipient_data.df$Unique_mouse_ID <- factor(
  Recipient_data.df$Unique_mouse_ID,
  levels = sort(unique(as.numeric(as.character(
    Recipient_data.df$Unique_mouse_ID
  ))))
)

#Graphing shannon diversity across timepoints
Shannon_timeline_plot <- ggplot(Recipient_data.df,
                                aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(shape = Study_ID, colour = Group), size = 2.5, alpha = 0.9) +
  theme_minimal(base_size = 14) +
#  scale_colour_manual(values = mouse_colours) +
  labs(title = "Shannon Diversity per Timepoint",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  theme(legend.position = "right")

Shannon_timeline_plot

#Graphing shannon diversity across timepoints by study ID
Shannon_timeline_plot_a <- ggplot(Recipient_data.df %>%
                                  filter(Study_ID %in% ("MBC1")), 
                                aes(x = Timepoint_ordered, y = Shannon, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(colour = Unique_mouse_ID, shape = Group), size = 2.5, alpha = 0.9) +
  #  scale_colour_manual(values = mouse_colours) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.line = element_line(
          linewidth = 0.5,
          linetype = "solid"
        ))+
  labs(title = "Study MBC1 diversity",
       x = "Timepoint",
       y = "Shannon diversity Index") +
  guides(fill = "none")

Shannon_timeline_plot_a

#Graphing observed diversity
Observed_timeline_plot <- ggplot(Recipient_data.df, 
                                 aes(x = Timepoint_ordered, y = Observed, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(shape = Study_ID, colour = Group), size = 2.5, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(title = "Observed Diversity per Timepoint",
       x = "Timepoint",
       y = "Observed Diversity") +
  theme(legend.position = "right")

Observed_timeline_plot

Observed_timeline_plot_a <- ggplot(Recipient_data.df %>%
                                    filter(Study_ID %in% ("MBC2")), 
                                  aes(x = Timepoint_ordered, y = Observed, fill = Timepoint_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = Unique_mouse_ID), alpha = 0.2, linewidth = 0.2) +
  geom_point(aes(colour = Unique_mouse_ID, shape = Group), size = 2.5, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  #  scale_colour_manual(values = mouse_colours) +
  labs(title = "Shannon Diversity MBC2",
       x = "Timepoint",
       y = "Shannon Diversity Index") +
  guides(fill = "none") +
  theme(legend.position = "right")

Observed_timeline_plot_a

#====Graphing beta diversity====
#----Bray curtis and Aitchison global statistics----

#create metadata that has no donor samples
metadata_recipients <- metadata_clean %>%
  filter(Group %in% c("allo", "syn"))
#cleanup - make the first column 'SampleID' as the rownames
#library(tibble)
#metadata_recipients <- metadata_recipients %>%
#  column_to_rownames("SampleID")

#align with braycurtis_dist
braycurtis_matrix <- as.matrix(braycurtis_dist)
samples <- intersect(rownames(braycurtis_matrix), rownames(metadata_recipients))
metadata_recipients <- metadata_recipients[samples, ]
braycurtis_subset <- as.dist(braycurtis_matrix[samples, samples])

#Repeat for Aitchison - metadata_recipients with aitchison distance 
aitchison_matrix<- as.matrix(aitchison_dist)
samples <- intersect(rownames(aitchison_matrix), rownames(metadata_recipients))
metadata_recipients <- metadata_recipients[samples, ]
aitchison_subset <- as.dist(aitchison_matrix[samples, samples])

#make sure timepoint is a factor, not numeric
metadata_recipients$Timepoint <- as.factor(metadata_recipients$Timepoint)
#Double check group/time balance - both groups are present across time - yes
table(metadata_recipients$Group, metadata_recipients$Timepoint)

#make sure unique mouse ID are in order when graphed
metadata_recipients$Unique_mouse_ID <- factor(
  metadata_recipients$Unique_mouse_ID,
  levels = sort(unique(as.numeric(as.character(
    metadata_recipients$Unique_mouse_ID
  ))))
)

###Stats:(Question: is variation in gut microbiome explained by time? By Study_ID)
adonis2(braycurtis_subset ~ Study_ID, #replace here with aitchison_subset or braycurtis_subset
        data = metadata_recipients,
        permutations = 999,
        strata = metadata_recipients$MouseID) #permutations are constrained within each mouse/repeated measures

###Stats: (Question: does allo vs syn lead to different longitudinal changes?)
adonis2(aitchison_subset ~ Group * Timepoint, 
        data = metadata_recipients,
        permutations = 999,
        strata = metadata_recipients$MouseID, 
        by = "margin") #term level output

#checking inter dispersion # replace here with aitchison_subset
betadisper(braycurtis_subset, metadata_recipients$Group) %>% anova()
betadisper(braycurtis_subset, Recipient_data.df$Timepoint_ordered) %>% anova()

#Plotting dispersion of Timepoint_ordered
bd_timepoint_ordered <-betadisper(braycurtis_subset, Recipient_data.df$Timepoint_ordered)

boxplot(bd_timepoint_ordered)
plot(bd_timepoint_ordered)

#plotting dispersion of Group
bd_group <- betadisper(braycurtis_subset, Recipient_data.df$Group)

boxplot(bd_group)
plot(bd_group)

#can I say that donor type/group influences how the microbiome changes over time? - conclusion: Not at a global level

#----Bray curtis PCOA plot----
#Here we pair the global statistics with graphs
###plot: all samples coloured by timepoint
BC_plot_original <- ggplot(PC_data.df, aes(x = BC1, y = BC2, 
                                            colour = Timepoint_ordered, #Can change to Timepoint_ordered or Study_ID or X for what you want to show
                                            size = Group)) + 
  geom_point() +
  scale_size_manual(values = c(
    "allo"           = 2,
    "syn"            = 2,
    "allo_donor"     = 6,
    "syn_donor"      = 6
  ), guide = "none") + 
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA (coloured by Timepoint)",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    colour = "Timepoint_ordered"
  )

BC_plot_original

###Plot: individual timepoints
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
    geom_point(alpha = 0.8, stroke = 0.5, aes(shape = Study_ID)) +
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
                           ncol = 2,
                          labels = c("i", "ii", "iii", "iv", "v", "vi", "Vii"))
BC_PLOTS_all

###Plot:By study ID
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

#----Aitchisons MDS plot----
###plot: all samples coloured by timepoint
MDS_plot_original <- ggplot(PC_data.df, aes(x = MDS1, y = MDS2, 
                                           colour = Timepoint_ordered, #Can change to Timepoint_ordered or Study_ID or X for what you want to show
                                           size = Group)) + 
  geom_point() +
  scale_size_manual(values = c(
    "allo"           = 1,
    "syn"            = 1,
    "allo_donor"     = 1,
    "syn_donor"      = 1
  ), guide = "none") + 
  geom_point(
    data = PC_data.df %>% 
      filter(Group %in% c("allo","syn","allo_donor","syn_donor")),
    aes(shape = Group),
    size = 3,          # size has to be +1 from previous to overlay the shape
    stroke = 1
  )+
  scale_shape_manual(values = c(
    "allo" = 16,            # circle
    "syn" = 17,             # triangle
    "allo_donor" = 15,
    "syn_donor" = 14
  )) +
  scale_color_viridis_d(option = "plasma") +  # discrete palette for ordered factors such as (Timepoint_ordered)
  theme_minimal(base_size = 14) +
  labs(
    title = "Aitchisons MDS (study ID)",
    x = xlab_mds1,
    y = ylab_mds2,
    colour = "Timepoint_ordered"
  )

MDS_plot_original

#relative taxon ratios are vastly different between studies...means it is not a good measure?
library(cowplot)

A_PLOTS <- list()

for(tp in timepoints){
  df_tp <- PC_data.df %>% filter(Timepoint_ordered == tp)
  
  A_PLOTS[[tp]] <- ggplot(df_tp, aes(x = MDS1, y = MDS2, colour = Group, size = Group)) +
    geom_point(alpha = 0.8, stroke = 0.5, aes(shape = Study_ID)) +
    stat_ellipse(level = 0.95, size = 0.5) +
    # overlay donor samples: keep group color, but larger + black outline
    geom_point(data = PC_data.df %>% filter(SampleID %in% donor_samples),
               aes(x = MDS1, y = MDS2, colour = Group), 
               size = 5, shape = 21, stroke = 1.5, fill = NA) +
    scale_color_manual(values = group_cols) +
    scale_size_manual(values = c(
      "allo"           = 2,
      "syn"            = 2,
      "allo_donor"     = 6,
      "syn_donor"      = 6
    ), guide = "none") + 
    theme_bw(base_size = 14) +
    labs(title = paste("Aitchison MDS:", tp),
         x = xlab_mds1,
         y = ylab_mds2)
}
A_PLOTS

# Extract legend from first plot
shared_legend <- get_legend(
  A_PLOTS[[1]] +
    theme(legend.position = "bottom")
)

A_PLOTS_nolegend <- lapply(A_PLOTS, function(p) {
  p + theme(legend.position = "none")
})

A_PLOTS_all <- plot_grid(
  plot_grid(
    plotlist = A_PLOTS_nolegend[timepoints],
    ncol = 2,
    labels = c("i", "ii", "iii", "iv", "v", "vi", "vii")
  ),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

A_PLOTS_all

#----Distance to donor microbiome----
donor_samples <- c("25-06236","25-06237","25-06294","25-06295")
#Assign in metadata_clean - could do on excel directly
metadata_clean$Donor_ID[metadata_clean$SampleID == "25-06236"] <- "A"
metadata_clean$Donor_ID[metadata_clean$SampleID == "25-06237"] <- "B"
metadata_clean$Donor_ID[metadata_clean$SampleID == "25-06294"] <- "C"
metadata_clean$Donor_ID[metadata_clean$SampleID == "25-06295"] <- "D"

#cleanup - make the first column 'SampleID' as the rownames
library(tibble)
metadata_clean <- metadata_clean %>%
  column_to_rownames("SampleID")

#check match up 
table(metadata_clean$Group, metadata_clean$Donor_ID)
#check no missing values (should return FALSE)
any(is.na(metadata_recipients$dist_to_donor))

### Extract dist_to_donor calculations from _matrix and add to object "metadata_recipients"
#use _matrix that was generated during global statistics
#i.e. braycurtis_matrix <- as.matrix(braycurtis_dist)

#function uses both metadata_recipients and metadata_clean (contains donors) to identify recipients and donor
metadata_recipients$dist_to_donor <- sapply(rownames(metadata_recipients), function(sid) {
  
  donor_id <- metadata_recipients[sid, "Donor_ID"]
  
  donor_samples <- rownames(metadata_clean)[
    metadata_clean$Donor_ID == donor_id & 
      metadata_clean$Group %in% c("allo_donor", "syn_donor")
  ]
  
aitchison_matrix[sid, donor_samples]
})

#sanity check - did I extract properly? - go check the metadata_clean names of samples vs the output for 25-06238 in dist_to_donor matches braycurtis_matrix
braycurtis_matrix["25-06238","25-06236"] #donor A and matching recipient - correct value = 0.7631768
braycurtis_matrix["25-06237","25-06242"] #donor B and matching recipient - correct value = 0.5467520

#define time
metadata_recipients$Timepoint_ordered <- factor(metadata_recipients$Timepoint, 
                                       levels = c("d0","d4","d7","d14","d21","d28","d35"))

#Plot: distance to donor
distance_to_donor_plot <-ggplot(metadata_recipients, aes(Timepoint_ordered, dist_to_donor, 
                                color = Donor_ID, group = Unique_mouse_ID)) + #swap colour between Group and Donor_ID
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Donor_ID), fun = mean, geom = "line", linewidth = 1.5) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.line = element_line(
          linewidth = 0.5,
          linetype = "solid"
        ))+
  labs (title = "Distance to bone marrow donor over time",
        x = NULL,
        y = "Aitchison's distance")

distance_to_donor_plot

#Plot: distance to donor - individual

#Define distinguished colours for each mouse
mouse_ids <- unique(metadata_recipients$Unique_mouse_ID)

mouse_colours <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A","brown",
  "black", "#FDBF6F", "deeppink1", "#CAB2D6", "gold1",
  "maroon", "orchid1", "#FB9A99" , "blue1", "darkturquoise",
  "green1", "yellow4", "darkorange4"
)

distance_to_donor_plot_a <-ggplot(metadata_recipients %>%
                                    filter(Group %in% ("allo")),
                                  aes(Timepoint_ordered, dist_to_donor, 
                                                         color = Unique_mouse_ID, group = Unique_mouse_ID)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  scale_colour_manual(values = mouse_colours) +
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.line = element_line(
          linewidth = 0.5,
          linetype = "solid"
        ))+
  labs (title = "allo-BMT recipient mice distance to respective donor over time",
        x = NULL,
        y = "Bray-Curtis distance")

distance_to_donor_plot_a


###model it

library(lme4)

distance_to_donor_model <- lmer(dist_to_donor ~ Group * Timepoint + Study_ID + (1 | Unique_mouse_ID), 
              data = metadata_recipients)

summary(distance_to_donor_model)

#Simple linear model because Unique_mouse_ID variance = 0 
#(random effect contributes nothing, 
#consistent trajectory after accounting for Group, Timepoint and Study ID,
#not enough variability between mice)
lm(dist_to_donor ~ Group * Timepoint + Study_ID, data = metadata_recipients)

#----Distance to baseline microbiome---- 

## generate distance to baseline calculations
metadata_recipients$dist_to_baseline <- sapply(rownames(metadata_recipients), function(sid) {
  
  mouse_id <- metadata_recipients[sid, "Unique_mouse_ID"]
  
  # find that mouse's baseline sample (d0)
  baseline_sample <- rownames(metadata_recipients)[
    metadata_recipients$Unique_mouse_ID == mouse_id & 
      metadata_recipients$Timepoint == "d0"
  ]
  
  # safety check
  if (length(baseline_sample) != 1) return(NA)
  
  aitchison_matrix[sid, baseline_sample]
})

###Plot: distance to baseline

distance_to_baseline_plot <-ggplot(metadata_recipients, aes(Timepoint_ordered, dist_to_baseline, 
                                                         color = Donor_ID, group = Unique_mouse_ID)) + #swap colour between Group and Donor_ID
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Donor_ID), fun = mean, geom = "line", linewidth = 1.5) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.line = element_line(
          linewidth = 0.5,
          linetype = "solid"
        ))+
  labs (title = "Distnace to baseline over time",
        y= "Aitchison's distance")

distance_to_baseline_plot

###Stats
distance_to_baseline_model <- lmer(dist_to_baseline ~ Group * Timepoint + Study_ID + (1 | Unique_mouse_ID), 
                       data = metadata_recipients)

summary(distance_to_baseline_model)

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

#save plots in order of production
save_plot(Shannon_timeline_plot, "Shannon_timeline_plot.pdf")
save_plot(Shannon_timeline_plot_a, "Shannon_timeline_plot_study-MBC1.pdf", width = 6, height = 4)
save_plot(Observed_timeline_plot, "Observed_timeline_plot.pdf")

save_plot(BC_plot_original, "BC_plot_all_timepoint.pdf")
save_plot(BC_PLOTS_all,"BrayCurtis_PCoA_allTimepoints.pdf", width = 10, height = 12)
save_plot(BC_PLOTS_MBC1_all, "BrayCurtis_PCoA_MBC1.pdf", width = 12, height = 10)
save_plot(BC_PLOTS_MBC2_all, "BrayCurtis_PCoA_MBC2.pdf", width = 12, height = 10)

save_plot(MDS_plot_original, "MDS_plot_original.pdf")
save_plot(A_PLOTS_all, "A_PLOTS_all.pdf", width = 10, height = 12)

save_plot(distance_to_donor_plot, "Distance_to_donor_plot.pdf", width = 6, height = 4)
save_plot(distance_to_donor_plot_a, "allo_distance_to_donor_plot.pdf", width= 6, height = 3)
save_plot(distance_to_baseline_plot, "Distance_to_baseline_plot.pdf", width = 6, height = 4)
