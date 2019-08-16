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
library(RColorBrewer)
library(ggplot2)

## Total tree
isol_info <- read.table("data/id.txt",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE)

mlst_info <- read.table("data/mlst_results.txt",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE)

isolate_data <- read.table("data/isolate_data.txt",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE)

metadata <- isolate_data[,c("id","species")] %>%
  left_join(mlst_info[,c("id","ST")], by = "id") %>%
  rename("Species" = species)

heatmap_data <- read.table("data/chewbbaca/complete_results/tree_heatmap_data.txt",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE)

iqtree <- read.iqtree("data/Phylogeny/Roary/iqtree_roary.contree")
iqtree_clean <- fix_tip_labels(iqtree, isol_info, "_pilon_spades.fasta", tree_type = "treedata")


palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43",
             "Wild bird" = "#fdae61")



palette2 <- c("> 95" = "#6a51a3",
                  "80 - 94" = "#9e9ac8",
                  "50 - 79" = "#cbc9e2",
                  "< 50" = "#f2f0f7",
                  "1-A" = "#fc8d59",
                  "1-I" = "#80b1d3",
                  "0" = "grey95")

iqtree_clean@phylo <- midpoint(iqtree_clean@phylo, labels = "support")

iqtree_clean@data$Bootstrap <- factor(
  case_when(
    iqtree_clean@data$UFboot >= 95 ~ "> 95",
    iqtree_clean@data$UFboot < 95 &
      iqtree_clean@data$UFboot >= 80 ~ "80 - 94",
    iqtree_clean@data$UFboot < 80 &
      iqtree_clean@data$UFboot >= 50 ~ "50 - 79",
    iqtree_clean@data$UFboot < 50 ~ "< 50"
  ),
  ordered = TRUE,
  levels = c("> 95",
             "80 - 94",
             "50 - 79",
             "< 50")
)

p1 <- ggtree(iqtree_clean,
       layout = "circular",
       size = 0.3) %<+% metadata +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21,
                 color = "white",
                 size = 1) +
  geom_tiplab2(aes(label = ST),
               align = TRUE,
               size = 1.5,
               offset = 0.01,
               linesize = 0.2) +
  geom_tippoint(aes(color = Species),
                size = 1) +
  guides(colour = guide_legend(override.aes = list(size=5)),
         fill = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette2) +
  theme(legend.position = "right")

p2 <- add_heatmap(p1,
            "data/chewbbaca/complete_results/tree_heatmap_data.txt",
            font_size = 2,
            colnames_offset = 3.3,
            heatmap_offset = 0.025,
            heatmap_width = 0.3) +
  scale_fill_manual(values = palette2) +
  geom_treescale(x = 0.05, y = 0, offset = 2, fontsize = 1)


ggsave("figures/total_SNP_tree.png",
       p2,
       device = "png",
       units = "cm",
       dpi = 600,
       height = 25,
       width = 25)




isol_info <- read.table("data/id.txt",
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE)

ST117_raw <- read.iqtree("data/Phylogeny/STs/ST117/iqtree/iqtree.contree")
ST162_raw <- read.iqtree("data/Phylogeny/STs/ST162/iqtree/iqtree.contree")

raw_trees <- list(ST117_raw,
                  ST162_raw)

names(raw_trees) <- c("Clade F","Clade A")

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
  "Broiler" = 21,
  "Pig" = 22,
  "Red fox" = 23,
  "Wild bird" = 24
)

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

boot_palette <- c("> 95" = "#6a51a3",
              "80 - 94" = "#9e9ac8",
              "50 - 79" = "#cbc9e2",
              "< 50" = "#f2f0f7")

clean_trees_new <- lapply(clean_trees, function(x){
  x@data$Bootstrap <- factor(case_when(x@data$UFboot >= 95 ~ "> 95",
                                x@data$UFboot <95 & x@data$UFboot >= 80 ~ "80 - 94",
                                x@data$UFboot < 80 & x@data$UFboot >= 50 ~ "50 - 79",
                                x@data$UFboot < 50 ~ "< 50"),
                             ordered = TRUE,
                             levels = c("> 95", "80 - 94", "50 - 79", "< 50"))
  
  return(x)
})

p1 <- annotate_tree(
  clean_trees_new$`Clade F`,
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
  cladelabel_offset = 0.00025,
  bootstrap_lab = FALSE,
  bootstrap_var = "Bootstrap",
  nodepoint_size = 1.5,
  tippoint_size = 4,
  label_offset = 0.00005,
  color_palette = palette_color,
  node_palette = boot_palette,
  shape_palette = palette_shape
) +
  xlim(0, 0.0025) +
  ggtitle("Clade F, ST117")

p2 <- annotate_tree(
  clean_trees_new$`Clade A`,
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
  cladelabel_offset = 0.0055,
  bootstrap_lab = FALSE,
  bootstrap_var = "Bootstrap",
  nodepoint_size = 1.5,
  tippoint_size = 4,
  label_offset = 0.001,
  color_palette = palette_color,
  node_palette = boot_palette,
  shape_palette = palette_shape
) +
  xlim(0, 0.047) +
  ggtitle("Clade A, ST162")

ggsave("figures/CladeA_ST162.png",
       p2,
       device = "png",
       dpi = 600,
       units = "cm",
       height = 15,
       width = 15)

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
