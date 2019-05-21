## Script for generating figures for the article

library(ggplot2)
library(dplyr)
library(ggtree)
library(tibble)
library(vegan)
library(cluster)
library(ape)

## Total dendrogram
cgMLST_data <- read.table("data/chewbbaca/complete_results/clean_cgMLST_data.txt",
                          sep = "\t",
                          header = TRUE,
                          colClasses = "factor")

total_tree <- as.phylo(hclust(daisy(cgMLST_data, metric = "gower"), method = "average"))


tree_metadata <- read.table("data/chewbbaca/complete_results/tree_metadata.txt",
                            sep = "\t",
                            header = TRUE)

tree_heatmap_data <- read.table("data/chewbbaca/complete_results/tree_heatmap_data.txt",
                                sep = "\t",
                                header = TRUE)


palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43",
             "Wild bird" = "#fdae61")

palette2 <- c("1-A" = "#fc8d59",
              "1-I" = "#80b1d3",
              "0" = "grey95")

total_tree_annotated <- ggtree(total_tree,
                               layout = "circular",
                               color = "#808080",
                               size = 0.5) %<+% tree_metadata +
  geom_tippoint(aes(color = species),
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.01,
               size = 2.1) +
  scale_color_manual(values = palette)

total_tree_opened <- rotate_tree(open_tree(total_tree_annotated, 8), 94)

total_tree_heatmap <- gheatmap(
  total_tree_opened,
  tree_heatmap_data,
  offset = 0.04,
  width = 0.3,
  colnames_offset_y = 2.5,
  colnames_position = "top",
  font.size = 3
) +
  scale_fill_manual(
    values = palette2,
    labels = c(
      "Gene/Mutation absent",
      "Acquired gene present",
      "Intrinsic gene mutated"
    )
  ) +
  theme(
    legend.position = c(0.9, 0.9),
    legend.text = element_text(size = 10)
  ) +
  
  guides(colour = guide_legend(override.aes = list(size=5)))


ggsave(
  "C:/Users/vi1511.VETINST/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/total_dendrogram.tiff",
  total_tree_heatmap,
  device = "tiff",
  dpi = 600,
  units = "cm",
  height = 30,
  width = 30
)


## Dendrogram per species, detailed

metadata_species <- read.table("data/chewbbaca/complete_results/metadata_species_specific_trees.txt",
                               sep = "\t",
                               header = TRUE) %>%
  select(gyrA, parC, parE, PMQR) %>%
  dplyr::rename("GyrA" = gyrA,
                "ParC" = parC,
                "ParE" = parE)

palette3 <- c("0" = "grey95",
              # gyrA
              "S83L" = "#c6dbef",
              "S83A" = "#9ecae1",
              "D87N" = "#6baed6",
              "D87H" = "#4292c6",
              "D87Y" = "#2171b5",
              "D87G" = "#08519c",
              "S83L, D87N" = "#08306b",
              # parC
              "S80I" = "#a1d99b",
              "S80R" = "#74c476",
              "S57T" = "#41ab5d",
              "A56T, S80I" = "#238b45",
              "S80I, E84V" = "#005a32",
              # parE
              "D475E" = "#fdd0a2",
              "D463N" = "#fdae6b",
              "S458A" = "#fd8d3c",
              "L416F" = "#f16913",
              "H516Y" = "#d94801",
              "L488M, A512T" = "#8c2d04",
              # PMQR
              "qnrA1" = "#dadaeb",
              "qnrB19" = "#bcbddc",
              "qnrS1" = "#9e9ac8",
              "qnrS2" = "#807dba",
              "qnrS4" = "#6a51a3",
              "qepA4" = "#4a1486")

# Broiler tree

broiler_tree <- calc_tree("data/chewbbaca/complete_results/cgMLST_broiler.tsv")

broiler_tree_annotated <- ggtree(broiler_tree,
                                 layout = "circular",
                                 color = "#808080",
                                 size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#808080",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.02,
               size = 2)

broiler_complete <- gheatmap(broiler_tree_annotated,
                             metadata_species,
                             offset = 0.1,
                             width = 0.5,
                             colnames = FALSE) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Pig tree
pig_tree <- calc_tree("data/chewbbaca/complete_results/cgMLST_pig.tsv")

pig_tree_annotated <- ggtree(pig_tree,
                             layout = "circular",
                             color = "#808080",
                             size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#808080",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.02,
               size = 2)

pig_complete <- gheatmap(pig_tree_annotated,
                         metadata_species,
                         offset = 0.09,
                         width = 0.5,
                         colnames = FALSE) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Red fox tree
fox_tree <- calc_tree("data/chewbbaca/complete_results/cgMLST_fox.tsv")

fox_tree_annotated <- ggtree(fox_tree,
                             layout = "circular",
                             color = "#808080",
                             size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#808080",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.02,
               size = 2)

fox_complete <- gheatmap(fox_tree_annotated,
                         metadata_species,
                         offset = 0.09,
                         width = 0.5,
                         colnames = FALSE) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Wild bird tree
bird_tree <- calc_tree("data/chewbbaca/complete_results/cgMLST_bird.tsv")

bird_tree_annotated <- ggtree(bird_tree,
                              layout = "circular",
                              color = "#808080",
                              size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#808080",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.02,
               size = 2)

bird_complete <- gheatmap(bird_tree_annotated,
                          metadata_species,
                          offset = 0.09,
                          width = 0.5,
                          colnames = FALSE) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Total plot
tot_plot <- plot_grid(broiler_complete,
                      pig_complete,
                      fox_complete,
                      bird_complete,
                      nrow = 2,
                      ncol = 2,
                      align = "h",
                      labels = c("Broiler","Pig","Red fox","Wild bird"))
ggsave(
  "C:/Users/vi1511.VETINST/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/test.tiff",
  tot_plot,
  device = "tiff",
  dpi = 600,
  units = "cm",
  height = 25,
  width = 25
)


# Other figures

## res per species
gene_names <- names(mut_report)

gene_names <- gene_names[-c(1,2)]

palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43",
             "Wild bird" = "#fdae61")

ggsave(
 "C:/Users/vi1511.VETINST/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/res_per_species.tiff",
  mut_report %>%
   filter(id %not_in% ex_samples) %>%
   left_join(isolate_data, by = "id") %>%
   gather(key, value, gene_names) %>%
   group_by(key, value, species) %>%
   dplyr::count() %>%
   ungroup() %>%
   mutate(value = if_else(value == 1, "Present", "Absent")) %>%
   spread(value, n, fill = 0) %>%
   rowwise() %>%
   mutate(
     Total = Present + Absent,
     Percent = round(Present / Total * 100, 1),
     lwr = round(get_binCI(Present, Total)[1], 1),
     upr = round(get_binCI(Present, Total)[2], 1)
   ) %>%
   mutate(species = factor(species,
                           levels = c("Broiler","Pig","Wild bird","Red fox")),
          key = paste0(toupper(substr(key, 1, 1)),
                       substr(key, 2, nchar(key))),
          group = if_else(key %in% c("GyrA","GyrB","ParC","ParE"), 1, 2)) %>%
   ggplot(aes(key, Percent, fill = species)) +
   geom_col(color = "black",
            position = position_dodge(1)) +
   geom_errorbar(
     aes(ymin = lwr, ymax = upr),
     width = 0.6,
     position = position_dodge(1)
   ) +
   scale_fill_manual(values = palette) +
   scale_x_discrete(limits = c("GyrA","GyrB","ParC","ParE",
                               "MarA","MarR","SoxR","RpoB",
                               "RobA")) +
   labs(y = "Percent (%) of isolates") +
   guides(fill = FALSE) +
   theme_classic() +
   theme(
     axis.text.x = element_text(
       angle = 90,
       size = 12,
       hjust = 1,
       vjust = 0.4
     ),
     axis.title.x = element_blank(),
     axis.text.y = element_text(size = 12),
     axis.title.y = element_text(size = 14),
     strip.text = element_text(size = 12)
   ) +
   facet_wrap(~ species),
  dpi = 600,
  units = "cm",
  device = "tiff",
  height = 20,
  width = 25
)

# acquired per species

qnr_genes <- c("qnrA1","qnrB19","qnrS1","qnrS2","qnrS4","qepA4")

ac_plot <- acquired_report %>%
  left_join(isolate_data, by = "id") %>%
  gather(key, value, qnr_genes) %>%
  group_by(key, value, species) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(value = if_else(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  rowwise() %>%
  mutate(
    Total = Present + Absent,
    Percent = round(Present / Total * 100, 1),
    lwr = round(get_binCI(Present, Total)[1], 1),
    upr = round(get_binCI(Present, Total)[2], 1)
  ) %>%
  ggplot(aes(key, Percent, fill = species)) +
  geom_col(color = "black",
           position = position_dodge(0.9)) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    width = 0.6,
    position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = palette) +
  labs(y = "Percent (%) of isolates",
       fill = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(face = "italic"),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.justification = c(0, 1),
    legend.position = c(0.82, 0.97)
  )

ggsave(
  "C:/Users/vi1511.VETINST/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/acquired_per_species.tiff",
  ac_plot,
  device = "tiff",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 25
)

# NMDS analysis

res_mechanisms <- read.table("data/MDS_data.txt",
                             sep = "\t",
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
  column_to_rownames("id")

bin_dist <- dist(res_mechanisms,
                 method = "binary")

nmds_df <- metaMDS(bin_dist,
                   trymax = 200,
                   trace = FALSE) # silence the output

plot_data <- as.data.frame(nmds_df$points) %>%
  rownames_to_column("ID") %>%
  mutate(ST = sub("_.+","",ID),
         ST = factor(ST, levels = c("10","58","117","155","162","355")))

ggplot(plot_data, 
       aes(x = MDS1,
           y = MDS2,
           fill = ST)) +
  geom_point(position = position_jitter(0.015, 0.015),
             size = 3,
             pch = 21) +
  scale_fill_brewer(type = "qual",
                    palette = 6) +
  labs(x = "Dimension 1",
       y = "Dimension 2",
       fill = "Sequence\nType") +
  theme_bw() +
  theme(axis.line = element_line(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.border = element_blank())

ggsave(
  "C:/Users/vi1511.VETINST/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/nmds_plot.tiff",
  nmds_plot,
  device = "tiff",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 25)
