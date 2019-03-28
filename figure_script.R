## Script for generating figures for the article

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
                                header = TRUE) %>%
  select(-gadX)

names(tree_heatmap_data) <- c("GyrA","ParC","ParE","RpoB","MarA","MarR","SoxR","qepA4","qnr")

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
               size = 2.5) +
  scale_color_manual(values = palette)

total_tree_opened <- rotate_tree(open_tree(total_tree_annotated, 8), 93)

ggsave(
  "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/total_dendrogram.tiff",
  gheatmap(
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
    ),
  device = "tiff",
  dpi = 600,
  units = "cm",
  height = 30,
  width = 30
)


## Dendrogram per species, detailed

metadata_species <- read.table("data/chewbbaca/complete_results/metadata_species_specific_trees.txt",
                               sep = "\t",
                               header = TRUE)

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
              "qnrA1" = "#762a83",
              "qnrB19" = "#ef8a62",
              "qnrS1" = "#a6dba0",
              "qnrS2" = "#5aae61",
              "qnrS4" = "#1b7837",
              "qepA4" = "#80cdc1",
              # num
              "1" = "#440154FF",
              "2" = "#3B528BFF",
              "3" = "#21908CFF",
              "4" = "#5DC863FF",
              "5" = "#FDE725FF")

# Broiler tree
broiler_cgMLST <- read.table("data/chewbbaca/complete_results/broiler_cgmlst.txt",
                             sep = "\t",
                             header = TRUE,
                             colClasses = "factor")

broiler_tree <- as.phylo(hclust(daisy(broiler_cgMLST, metric = "gower"), method = "average"))

broiler_tree_annotated <- ggtree(broiler_tree,
                                 layout = "circular",
                                 color = "#808080",
                                 size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#4575b4",
                size = 1) +
  geom_tiplab2(aes(label = ST),
               offset = 0.03,
               size = 1.6)

broiler_tree_opened <- rotate_tree(open_tree(broiler_tree_annotated, 12), 93.5)

broiler_complete <- gheatmap(broiler_tree_opened,
                             metadata_species,
                             offset = 0.12,
                             width = 0.9,
                             colnames_offset_y = 0.9,
                             colnames_position = "top",
                             font.size = 2.5) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Pig tree
pig_cgMLST <- read.table("data/chewbbaca/complete_results/pig_cgmlst.txt",
                         sep = "\t",
                         header = TRUE,
                         colClasses = "factor")

pig_tree <- as.phylo(hclust(daisy(pig_cgMLST, metric = "gower"), method = "average"))

pig_tree_annotated <- ggtree(pig_tree,
                             layout = "circular",
                             color = "#808080",
                             size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#74add1",
                size = 1.3) +
  geom_tiplab2(aes(label = ST),
               offset = 0.03,
               size = 2)

pig_tree_opened <- rotate_tree(open_tree(pig_tree_annotated, 12), 94)

pig_complete <- gheatmap(pig_tree_opened,
                         metadata_species,
                         offset = 0.14,
                         width = 0.9,
                         colnames_offset_y = 0.69,
                         colnames_position = "top",
                         font.size = 2.5) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Red fox tree
fox_cgMLST <- read.table("data/chewbbaca/complete_results/fox_cgmlst.txt",
                         sep = "\t",
                         header = TRUE,
                         colClasses = "factor")

fox_tree <- as.phylo(hclust(daisy(fox_cgMLST, metric = "gower"), method = "average"))

fox_tree_annotated <- ggtree(fox_tree,
                             layout = "circular",
                             color = "#808080",
                             size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#f46d43",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.03,
               size = 2)

fox_tree_opened <- rotate_tree(open_tree(fox_tree_annotated, 12), 93.5)

fox_complete <- gheatmap(fox_tree_opened,
                         metadata_species,
                         offset = 0.15,
                         width = 0.9,
                         colnames_offset_y = 0.34,
                         colnames_position = "top",
                         font.size = 2.5) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

# Wild bird tree
bird_cgMLST <- read.table("data/chewbbaca/complete_results/bird_cgmlst.txt",
                          sep = "\t",
                          header = TRUE,
                          colClasses = "factor")

bird_tree <- as.phylo(hclust(daisy(bird_cgMLST, metric = "gower"), method = "average"))

bird_tree_annotated <- ggtree(bird_tree,
                              layout = "circular",
                              color = "#808080",
                              size = 0.5) %<+% tree_metadata +
  geom_tippoint(color = "#fdae61",
                size = 1.5) +
  geom_tiplab2(aes(label = ST),
               offset = 0.03,
               size = 2)

bird_tree_opened <- rotate_tree(open_tree(bird_tree_annotated, 12), 93.5)

bird_complete <- gheatmap(bird_tree_opened,
                          metadata_species,
                          offset = 0.15,
                          width = 0.9,
                          colnames_offset_y = 0.54,
                          colnames_position = "top",
                          font.size = 2.5) +
  scale_fill_manual(values = palette3) +
  guides(fill = FALSE)

tot_plot <- plot_grid(broiler_complete,
                      pig_complete,
                      fox_complete,
                      bird_complete,
                      nrow = 2,
                      ncol = 2,
                      align = "h",
                      labels = c("Broiler","Pig","Red fox","Wild bird"))

ggsave(
  "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/species_dendrogram.tiff",
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
 "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/res_per_species.tiff",
 mut_report %>%
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

qnr_genes <- c("qnrA1","qnrB19","qnrS1","qnrS2","qnrS4")

ggsave(
  "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/acquired_per_species.tiff",
  acquired_report %>%
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
    ),
  device = "tiff",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 25
)



# Mechanisms per animal species and MIC

AddCol <- function(df, col_name) {
  # split rows by delimeters
  string_to_proc <- df %>% select(!!col_name) %>%
    unlist() %>% str_split(regex("\\, |\\,")) 
  # find unique entries
  unique_strings <- string_to_proc %>%
    unlist() %>% unique()
  # construct names of the new columns
  cols_names <- paste(col_name, unique_strings, sep = "_")
  # construct 0/1-content columns for each unique entry
  cols_content <- sapply(function(i) {
    as.integer(unlist(lapply(function(Z) any(Z %in% unique_strings[i]), 
                             X = string_to_proc)))
  }, X = seq_along(unique_strings))
  res <- data.frame(cols_content)
  names(res) <- cols_names
  return(res)
}

total_df <- acquired_report %>%
  left_join(mut_quant[c("gyrA","gyrB","parC","parE","id")], by = "id") %>%
  left_join(mut_report[c("id","marR","marA","rpoB","soxR","robA")], by = "id") %>%
  left_join(isolate_data[c("id", "species")])

names(total_df) <- c("id","qepA4","qnrA1","qnrB19",
                     "qnrB6","qnrB60","qnrS1","qnrS2",
                     "qnrS4","GyrA","GyrB","ParC",
                     "ParE","MarR","MarA","RpoB",
                     "SoxR","RobA","species")

col_proc <- c("GyrA","GyrB","ParC","ParE")

id <- total_df$id

cols_to_add <- sapply(function(i) {AddCol(df = total_df, col_name = col_proc[i])},
                      X = seq_along(col_proc)) %>%
  bind_cols() %>%
  cbind(id) %>%
  gather(key, value, -id) %>%
  filter(grepl("(.*?)_0", key) == FALSE) %>%
  spread(key, value)

total_df2 <- total_df %>%
  select(-c(GyrA, GyrB, ParC, ParE)) %>%
  left_join(cols_to_add, by = "id") %>%
  select(id, species, everything()) %>%
  mutate_at(vars(qnrA1,qnrB19,qnrB6,qnrB60,qnrS1,qnrS2,qnrS4),
            funs(as.numeric(.))) %>%
  mutate(qnr = ifelse(rowSums(.[,4:10]) == 0, 0, 1)) %>%
  select(-c(qnrA1,qnrB19,qnrB6,qnrB60,qnrS1,qnrS2,qnrS4)) %>%
  left_join(mic_data[c("id","CIP")], by = "id") %>%
  select(-id) %>%
  mutate_at(vars(-species, CIP),
            funs(as.numeric(as.character(.)))) %>%
  group_by(species, CIP) %>%
  summarise_all(funs(sum, n = n())) %>%
  dplyr::rename("n" = MarR_n) %>%
  select(-contains("_n")) %>%
  gather(gene, value, -c(species, CIP, n)) %>%
  mutate(gene2 = sub("(.*?)_(.*?)_sum", "\\1", gene),
         gene2 = sub("(.*?)_sum", "\\1", gene2),
         gene = sub("(.*?)_(.*?)_sum", "\\2", gene),
         gene = sub("(.*?)_sum", "\\1", gene),
         type = case_when(gene2 == "GyrA" ~ 1,
                          gene2 == "GyrB" ~ 2,
                          gene2 == "ParC" ~ 3,
                          gene2 == "ParE" ~ 4,
                          gene2 == "MarR" ~ 5,
                          gene2 == "MarA" ~ 6,
                          gene2 == "RpoB" ~ 7,
                          gene2 == "SoxR" ~ 8,
                          gene2 == "RobA" ~ 9,
                          gene2 == "qnr" ~ 10,
                          gene2 == "qepA4" ~ 11),
         perc = value/n*100)



p1 <- ggplot(total_df2, aes(reorder(gene, type), factor(CIP), fill = perc))+
  geom_tile()+
  scale_fill_viridis(option = "A",
                     begin = 0.1,
                     guide = guide_colorbar(barheight = 10,
                                            barwidth = 0.8,
                                            label.position = "left",
                                            label.hjust = 0.5,
                                            label.vjust = 0.3,
                                            ticks = TRUE))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_rect(colour = "black", fill = "grey70"),
        strip.text = element_text(size = 12),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  labs(y = "CIP MIC")+
  facet_wrap(~species, ncol = 1)+
  coord_fixed()



nplot_data <- total_df %>%
  left_join(mic_data[c("id","CIP")], by = "id") %>%
  group_by(species, CIP) %>%
  dplyr::count() %>%
  ungroup()

n_names <- c("Broiler" = "n",
             "Pig" = "n",
             "Red fox" = "n",
             "Wild bird" = "n")

p2 <- ggplot(nplot_data, aes(factor(CIP), n))+
  geom_col()+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        strip.background = element_rect(colour = "black", fill = "grey70"),
        plot.margin = unit(c(0,0,-10,0), "cm"))+
  facet_wrap(~species, ncol = 1, labeller = as_labeller(n_names))+
  coord_flip()

p1g <- ggplotGrob(p1)
p2g <- ggplotGrob(p2)

g <- cbind(p1g, p2g, size = "first")

ggsave(
  "C:/Users/vi1511/OneDrive - Veterinærinstituttet/Artikler/qrec_wgs/figures/mechanisms_per_species_mic.tiff",
  grid.arrange(g),
  device = "tiff",
  units = "cm",
  dpi = 600,
  height = 25,
  width = 30)
