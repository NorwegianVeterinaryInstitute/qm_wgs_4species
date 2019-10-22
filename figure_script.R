## Script for generating figures for the article

library(ggplot2)
library(dplyr)
library(ggtree)
library(tibble)
library(vegan)
library(cluster)
library(ape)

## res per species
palette <- c("Acquired" = "#FB8C6F",
             "Intrinsic" = "#73607D")

acquired_report <- read.table("data/ariba/resfinder_quin_report.txt",
                              sep = "\t",
                              header = TRUE,
                              stringsAsFactors = FALSE) %>%
  select(id, contains("qnr"), qepA4)

acquired_report <- acquired_report[,colSums(acquired_report != 0) > 0]

mut_report <- read.table("data/ariba/megares_quin_report.txt",
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE)

intrinsic_genes <- names(mut_report)

test <- mut_report %>%
  left_join(acquired_report, by = "id") %>%
  select(-nt_change) %>%
  left_join(isolate_data[, c("id", "species")], by = "id") %>%
  gather(key, value,-c(id, species)) %>%
  group_by(key, value, species) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(value = if_else(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  rowwise() %>%
  mutate(Total = Present + Absent,
         Percent = round(Present / Total * 100, 1)) %>%
  mutate(
    key2 = ifelse(key %in% intrinsic_genes, paste0(toupper(substr(key, 1, 1)), substr(key, 2, nchar(key))), key),
    species = factor(species,
                     levels = c("Broiler", "Pig", "Wild bird", "Red fox")),
    group = if_else(
      key %in% intrinsic_genes,
      "Intrinsic",
      "Acquired"
    )
  )

res_plot <- ggplot(test, aes(key2, Percent, fill = group)) +
  geom_col(color = "black",
           size = 1,
           position = position_dodge(1)) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(
    limits = c(
      "GyrA",
      "GyrB",
      "ParC",
      "ParE",
      "MarA",
      "MarR",
      "SoxR",
      "RpoB",
      "RobA",
      "qnrA1",
      "qnrB19",
      "qnrS1",
      "qnrS2",
      "qnrS4",
      "qepA4"
    ),
    labels = c(
      "GyrA",
      "GyrB",
      "ParC",
      "ParE",
      "MarA",
      "MarR",
      "SoxR",
      "RpoB",
      "RobA",
      expression(italic("qnrA1")),
      expression(italic("qnrB19")),
      expression(italic("qnrS1")),
      expression(italic("qnrS2")),
      expression(italic("qnrS4")),
      expression(italic("qepA4"))
    )
  ) +
  labs(y = "Percent (%) of isolates") +
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
    legend.title = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap( ~ species)

ggsave(
  "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/res_per_species.png",
  res_plot,
  dpi = 600,
  units = "cm",
  device = "png",
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
library(phangorn)
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



palette2 <- c(">= 95" = "white",
                  "< 95" = "black",
                  "1-A" = "#FB8C6F",
                  "1-I" = "#73607D",
                  "0" = "grey95")

iqtree_clean@phylo <- midpoint(iqtree_clean@phylo, labels = "support")

iqtree_clean@data$Bootstrap <- factor(
  case_when(
    iqtree_clean@data$UFboot >= 95 ~ ">= 95",
    iqtree_clean@data$UFboot < 95 ~ "< 95"
  ),
  ordered = TRUE,
  levels = c(">= 95",
             "< 95")
)

p1 <- ggtree(iqtree_clean,
       layout = "circular",
       size = 0.3) %<+% metadata +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21,
                 color = "black",
                 size = 0.5) +
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

heatmap_data <- read.table("data/chewbbaca/complete_results/tree_heatmap_data.txt",
                           sep = "\t",
                           header = TRUE)


rot_tree <- rotate_tree(open_tree(p1, 10), 95)



p2 <- gheatmap(rot_tree,
         heatmap_data,
         font.size = 2,
         offset = 0.025,
         colnames_offset_y = 3.3,
         colnames_position = "top",
         width = 0.3,
         color = NULL) +
  scale_fill_manual(values = palette2) +
  geom_treescale(x = 0.1, y = 286, offset = -7, fontsize = 3, linesize = 0.3)


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
ST162_B_raw <- read.iqtree("data/Phylogeny/STs/ST162_B/iqtree/iqtree.contree")

raw_trees <- list(ST117_raw,
                  ST162_B_raw)

names(raw_trees) <- c("Clade F","Clade B")

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

boot_palette <- c(">= 95" = "white",
                  "< 95" = "black")

palette_color <- c(palette_color, boot_palette)

clean_trees_new <- lapply(clean_trees, function(x){
  x@data$Bootstrap <- factor(case_when(x@data$UFboot >= 95 ~ ">= 95",
                                x@data$UFboot <95 ~ "< 95"),
                             ordered = TRUE,
                             levels = c(">= 95", "< 95"))
  
  return(x)
})


p1 <- ggtree(clean_trees_new$`Clade F`,
       layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.0001) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(26, label = "3 SNPs",
                  offset = 0.0005,
                  color = "grey50") +
  geom_treescale() +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette_color) +
  theme(legend.position = "right") +
  xlim(0, 0.004) +
  ggtitle("Clade F, ST117")

p2 <- ggtree(clean_trees_new$`Clade B`,
             layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.0003) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(16, label = "12 SNPs",
                  offset = 0.0018,
                  color = "grey50") +
  geom_treescale() +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette_color) +
  theme(legend.position = "right") +
  xlim(0, 0.015) +
  ggtitle("Clade B, ST162")


ggsave("figures/CladeF_ST117.png",
       p1,
       device = "png",
       dpi = 600,
       units = "cm",
       height = 15,
       width = 15)


## Suppl. trees
ST10_raw <- read.iqtree("data/Phylogeny/STs/ST10/iqtree/iqtree.contree")
ST162_A_raw <- read.iqtree("data/Phylogeny/STs/ST162_A/iqtree/iqtree.contree")
ST355_raw <- read.iqtree("data/Phylogeny/STs/ST355/iqtree/iqtree.contree")
ST744_raw <- read.iqtree("data/Phylogeny/STs/ST744/iqtree/iqtree.contree")

raw_trees <- list(ST10_raw,
                  ST162_A_raw,
                  ST355_raw,
                  ST744_raw)

names(raw_trees) <- c("Clade D","Clade A","Clade E","Clade C")

# Clean tree tip labels
clean_trees <- lapply(
  raw_trees,
  function(x) fix_tip_labels(
    x,
    id_data,
    "_pilon_spades.fasta",
    tree_type = "treedata"
  )
)

clean_trees <- lapply(clean_trees,
                      function(x) {
                        x@data$Bootstrap <- factor(
                          case_when(
                            x@data$UFboot >= 95 ~ ">= 95",
                            x@data$UFboot < 95 ~ "< 95"
                          ),
                          ordered = TRUE,
                          levels = c(">= 95",
                                     "< 95")
                        )
                        return(x)
                      })

# Specify color and shape palettes
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

boot_palette <- c(">= 95" = "white",
                  "< 95" = "black")

palette <- c(palette_color, boot_palette)


# Generate trees
options(scipen=999)

t_ST10 <- ggtree(clean_trees$`Clade D`,
       layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.0003) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(24, label = "40 SNPs",
                  offset = 0.0008,
                  color = "grey50") +
  geom_treescale() +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette) +
  theme(legend.position = "right") +
  xlim(0, 0.007)


t_162A <- ggtree(clean_trees$`Clade A`,
       layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.0003) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(26, label = "18 SNPs",
                  offset = 0.0016,
                  color = "grey50") +
  geom_treescale(x = 0.008) +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette) +
  theme(legend.position = "right") +
  xlim(0, 0.013)


t_ST355 <- ggtree(clean_trees$`Clade E`,
       layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.0003) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(26, label = "6 SNPs",
                  offset = 0.0025,
                  color = "grey50") +
  geom_treescale(x = 0.008) +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette) +
  theme(legend.position = "right") +
  xlim(0, 0.032)


t_ST744 <- ggtree(clean_trees$`Clade C`,
       layout = "rectangular") %<+% tree_metadata +
  geom_tiplab(aes(label = Location),
              offset = 0.00005) +
  geom_tippoint(aes(fill = Year,
                    shape = Species),
                size = 2.5) +
  geom_nodepoint(aes(fill = Bootstrap),
                 pch = 21) +
  geom_cladelabel(14, label = "50 SNPs",
                  offset = 0.0003,
                  color = "grey50") +
  geom_treescale() +
  scale_shape_manual(values = c("Broiler" = 21,
                                "Pig" = 22,
                                "Red fox" = 23,
                                "Wild bird" = 24)) +
  scale_fill_manual(values = palette) +
  theme(legend.position = "right") +
  xlim(0, 0.0032)


ggsave("figures/CladeC_ST744.png",
       t_ST744,
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
