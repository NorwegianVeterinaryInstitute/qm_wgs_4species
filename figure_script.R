## Script for generating figures for the article

library(ggplot2)
library(dplyr)
library(ggtree)
library(tibble)
library(vegan)
library(cluster)
library(ape)

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

# Phylogenetic trees

library(dplyr)
library(distanceR)
library(ape)
library(treeio)
library(ggtree)
library(phytools)
library(readxl)

isol_info <- read.table("data/id.txt",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE)

ST58_raw <- read.iqtree("data/Phylogeny/STs/ST58/iqtree/iqtree.contree")
ST117_raw <- read.iqtree("data/Phylogeny/STs/ST117/iqtree/iqtree.contree")
ST162_raw <- read.iqtree("data/Phylogeny/STs/ST162/iqtree/iqtree.contree")

raw_trees <- list(ST58_raw,
                  ST117_raw,
                  ST162_raw)

names(raw_trees) <- c("ST58","ST117","ST162")

clean_trees <- lapply(
  raw_trees,
  function(x) fix_tip_labels(
    x,
    isol_info,
    "_pilon_spades.fasta",
    tree_type = "treedata"
    )
  )

location_data <- read_excel("data/st_location_data.xlsx") %>%
  select(id, Kommune, Fylke) %>%
  mutate(fyl_id = group_indices(., Fylke),
         kom_id = group_indices(., Kommune)) %>%
  mutate(Location = paste(fyl_id, kom_id, sep = "-")) %>%
  select(id, Location)

tree_metadata <- read.table("data/chewbbaca/complete_results/tree_metadata.txt",
                            sep = "\t",
                            header = TRUE) %>%
  left_join(location_data, by = "id") %>%
  mutate(Year = substr(id, 1, 4)) %>%
  dplyr::rename("Species" = species)

palette_shape <- c(
  "Broiler" = 15,
  "Pig" = 18,
  "Red fox" = 16,
  "Wild bird" = 17
)

library(RColorBrewer)

palette_color <- brewer.pal(11, "Set3")

year_val <- c(
  "2006",
  "2007",
  "2008",
  "2009",
  "2010",
  "2011",
  "2012",
  "2013",
  "2014",
  "2015",
  "2016")

names(palette_color) <- year_val



library(ggplot2)

p1 <- annotate_tree(
  clean_trees$ST58,
  tree_metadata,
  layout = "rectangular",
  tree_type = "treedata",
  midroot = TRUE,
  line_width = 0.5,
  label_variable = "Location",
  color_variable = "Year",
  shape_variable = "Species",
  clade_label_node = 32,
  clade_label = "9 SNP diff",
  cladelabel_offset = 0.00020,
  bootstrap_lab = TRUE,
  bootlab_size = 2,
  tippoint_size = 4,
  label_offset = 0.00005,
  shape_palette = palette_shape,
  color_palette = palette_color
) +
  xlim(0, 0.002) +
  ggtitle("ST58 Clade")


p2 <- annotate_tree(
  clean_trees$ST117,
  tree_metadata,
  layout = "rectangular",
  tree_type = "treedata",
  midroot = TRUE,
  line_width = 0.5,
  label_variable = "Location",
  color_variable = "Year",
  shape_variable = "Species",
  clade_label_node = 23,
  clade_label = "37 SNP diff",
  cladelabel_offset = 0.00022,
  bootstrap_lab = TRUE,
  bootlab_size = 2,
  tippoint_size = 4,
  label_offset = 0.00005,
  shape_palette = palette_shape,
  color_palette = palette_color
) +
  xlim(0, 0.0025) +
  ggtitle("ST117 Clade")

p3 <- annotate_tree(
  clean_trees$ST162,
  tree_metadata,
  layout = "rectangular",
  tree_type = "treedata",
  midroot = TRUE,
  line_width = 0.5,
  label_variable = "Location",
  color_variable = "Year",
  shape_variable = "Species",
  clade_label_node = 29,
  clade_label = "12 SNP diff",
  cladelabel_offset = 0.0038,
  bootstrap_lab = TRUE,
  bootlab_size = 2,
  tippoint_size = 4,
  label_offset = 0.0005,
  shape_palette = palette_shape,
  color_palette = palette_color
) +
  xlim(0, 0.042) +
  ggtitle("ST162 Clade")



ggsave("ST117.png",
       p2,
       device = "png",
       dpi = 600,
       units = "cm",
       height = 25,
       width = 20)

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
