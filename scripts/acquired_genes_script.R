library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Import data

func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

ex_samples <- c(
  "2016-17-565-1-S",
  "2016-02-428-2-S",
  "2016-02-486-2-S",
  "2016-02-732-2-S",
  "2014-01-1741-1-S"
)

isolate_data <- read.table("data/isolate_data.txt",
                           sep = "\t",
                           header = T,
                           stringsAsFactors = F) %>%
  filter(!id %in% ex_samples)

id_data <- read.table("data/id.txt",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)

mic_data <- read.table("data/mic_data.txt",
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F) %>%
  filter(!id %in% ex_samples) %>%
  select(-antres) %>%
  mutate(SMX_R = if_else(SMX > 64, 1, 0),
         TMP_R = if_else(TMP > 2, 1, 0),
         CIP_R = if_else(CIP > 0.06, 1, 0),
         TET_R = if_else(TET > 8, 1, 0),
         NAL_R = if_else(NAL > 16, 1, 0),
         CTX_R = if_else(CTX > 0.25, 1, 0),
         CHL_R = if_else(CHL > 16, 1, 0),
         AMP_R = if_else(AMP > 8, 1, 0),
         GEN_R = if_else(GEN > 2, 1, 0)) %>%
  select(id, contains("_R"))

names(mic_data) <- sub("_R", "", names(mic_data))

acquired_data <- read_delim("data/VAMPIR/amr_ac/acquired_gene_report.tsv",
                            delim = "\t") %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref)



## Create groups from acquired gene data

aac3II_names <- grep("aac3", names(acquired_data), value = TRUE, ignore.case = TRUE)
aadA_names <- grep("aadA", names(acquired_data), value = TRUE, ignore.case = TRUE)
aph_names <- grep("aph", names(acquired_data), value = TRUE, ignore.case = TRUE)
blaCTXM_names <- grep("ctx", names(acquired_data), value = TRUE, ignore.case = TRUE)
blaSHV_names <- grep("shv", names(acquired_data), value = TRUE, ignore.case = TRUE)
blaTEM_names <- grep("tem", names(acquired_data), value = TRUE, ignore.case = TRUE)
dfrA_names <- grep("dfr", names(acquired_data), value = TRUE, ignore.case = TRUE)
qnr_names <- grep("qnr", names(acquired_data), value = TRUE, ignore.case = TRUE)
sul_names <- grep("sul", names(acquired_data), value = TRUE, ignore.case = TRUE)
tet_names <- grep("tet", names(acquired_data), value = TRUE, ignore.case = TRUE)

excluded_genes <- c("ermB",
                    "lnuF",
                    "mphA",
                    "mphB",
                    "mphE",
                    "msrE")

single_genes <-
  names(acquired_data)[!names(acquired_data) %in% c(
    aac3II_names,
    aadA_names,
    aph_names,
    blaCTXM_names,
    blaSHV_names,
    blaTEM_names,
    dfrA_names,
    qnr_names,
    sul_names,
    tet_names,
    excluded_genes,
    "id",
    "mdfA"
  )]

selected_cols <- c(
  "id",
  "aac3II_group",
  "aadA_group",
  "aph_group",
  "blaCTXM_group",
  "blaSHV_group",
  "blaTEM_group",
  "dfrA_group",
  "qnr_group",
  "sul_group",
  "tet_group",
  single_genes)


# Import intrinsic data
intrinsic_genes <- read_delim("data/ariba/mut_report.txt",
                              delim = "\t") %>%
  select(id, gyrA, parC, parE)


gene_order <- c("gyrA",
                "parC",
                "parE",
                "qnr",
                "qepA4",
                "blaCMY2",
                "blaCTXM",
                "blaSHV",
                "blaTEM",
                "aac3II",
                "aadA",
                "aph",
                "catA1",
                "cmlA1",
                "floR",
                "sul",
                "tet",
                "dfrA")

group_df <-
  acquired_data[!names(acquired_data) %in% excluded_genes] %>%
  mutate(
    aac3II_group = if_else(rowSums(.[, aac3II_names]) > 0, 1, 0),
    aadA_group = if_else(rowSums(.[, aadA_names]) > 0, 1, 0),
    aph_group = if_else(rowSums(.[, aph_names]) > 0, 1, 0),
    blaCTXM_group = if_else(rowSums(.[, blaCTXM_names]) > 0, 1, 0),
    blaSHV_group = if_else(rowSums(.[, blaSHV_names]) > 0, 1, 0),
    blaTEM_group = if_else(rowSums(.[, blaTEM_names]) > 0, 1, 0),
    dfrA_group = if_else(rowSums(.[, dfrA_names]) > 0, 1, 0),
    qnr_group = if_else(rowSums(.[, qnr_names]) > 0, 1, 0),
    sul_group = if_else(rowSums(.[, sul_names]) > 0, 1, 0),
    tet_group = if_else(rowSums(.[, tet_names]) > 0, 1, 0)
  ) %>%
  select(selected_cols) %>%
  left_join(intrinsic_genes, by = "id") %>%
  group_by_at(vars(-id)) %>%
  mutate(group = group_indices())


# Data wrangling

## Sorting data frame

sort_df <- mic_data %>%
  mutate(antres = rowSums(.[,2:10])) %>%
  select(id, antres) %>%
  left_join(group_df[,c("id","group")], by = "id") %>%
  group_by(group) %>%
  mutate(mean_antres = round(mean(antres), 1)) %>%
  select(group, mean_antres) %>%
  group_by(group) %>%
  summarise_all(list(func_paste)) %>%
  mutate(mean_antres = as.numeric(mean_antres))


# barplot data
barplot_df <- group_df[,c("id","group")] %>%
  left_join(isolate_data[,c("id","species")], by = "id") %>%
  group_by(group, species) %>%
  count() %>%
  group_by(group) %>%
  mutate(sum_n = sum(n)) %>%
  left_join(sort_df, by = "group")



palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43",
             "Wild bird" = "#fdae61")


# create gene presence/absence plot

##correct x-axis to same scale as plot above
scale_df <- barplot_df %>%
  select(group, sum_n) %>%
  group_by(group) %>%
  summarise_all(list(func_paste)) %>%
  mutate(sum_n = as.integer(sum_n))
  


gene_plot_df <- group_df %>%
  gather(gene, value, -c(id, group)) %>%
  mutate(gene = sub("_group", "", gene),
         gene = factor(gene, levels = gene_order),
         type = case_when(gene %in% c("gyrA","parC","parE","qnr","qepA4") ~ "Quinolones",
                          gene %in% c("aac3II","aadA","aph") ~ "Aminoglycosides",
                          gene %in% c("blaCTXM","blaCMY2","blaSHV","blaTEM") ~ "Beta lactamases",
                          gene %in% c("catA1","cmlA1","floR") ~ "Phenicols",
                          gene == "sul" ~ "Sulfonamides",
                          gene == "tet" ~ "Tetracyclines",
                          gene == "dfrA" ~ "Diaminopyrimidines")) %>%
  filter(value == 1) %>%
  left_join(scale_df, by = "group") %>%
  left_join(sort_df, by = "group")



shape_palette <- c("Intrinsic" = 23,
                   "Acquired" = 21)


type_palette <- c("Quinolones" = "#8dd3c7",
                  "Aminoglycosides" = "#ffffb3",
                  "Beta lactamases" = "#bebada",
                  "Phenicols" = "#fb8072",
                  "Sulfonamides" = "#80b1d3",
                  "Tetracyclines" = "#fdb462",
                  "Diaminopyrimidines" = "#b3de69")


## Resistance plot

AM_order <- c("CIP",
              "NAL",
              "AMP",
              "CTX",
              "GEN",
              "CHL",
              "SMX",
              "TET",
              "TMP")


res_plot_df <- mic_data %>%
  left_join(group_df[,c("id","group")], by = "id") %>%
  mutate(antres = rowSums(.[,2:10])) %>%
  select(group, everything(), -id) %>%
  gather(AM, value, -c(group, antres)) %>%
  filter(value == 1) %>%
  { . ->> res_groups } %>%
  group_by(AM, group) %>%
  count() %>%
  mutate(type = case_when(AM %in% c("CIP","NAL") ~ "Quinolones",
                          AM %in% c("GEN") ~ "Aminoglycosides",
                          AM %in% c("CTX","AMP") ~ "Beta lactamases",
                          AM %in% c("CHL") ~ "Phenicols",
                          AM == "SMX" ~ "Sulfonamides",
                          AM == "TET" ~ "Tetracyclines",
                          AM == "TMP" ~ "Diaminopyrimidines")) %>%
  left_join(scale_df, by = "group") %>%
  left_join(res_groups, by = c("group","AM")) %>%
  left_join(sort_df, by = "group") %>%
  ungroup() %>%
  mutate(AM = factor(AM, levels = AM_order))


type_palette <- c("Quinolones" = "#8dd3c7",
                  "Aminoglycosides" = "#ffffb3",
                  "Beta lactamases" = "#bebada",
                  "Phenicols" = "#fb8072",
                  "Sulfonamides" = "#80b1d3",
                  "Tetracyclines" = "#fdb462",
                  "Diaminopyrimidines" = "#b3de69")

# Plotting

bar_plot <- ggplot(barplot_df, aes(reorder(factor(group), mean_antres), n, fill = species)) +
  geom_col() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.background = element_blank()) +
  guides(fill = FALSE) +
  labs(y = "No. of isolates") +
  scale_fill_manual(values = palette)

gene_labels <- c(expression(italic("gyrA")),
                 expression(italic("parC")),
                 expression(italic("parE")),
                 expression(italic("qnr")),
                 expression(italic("qepA4")),
                 expression(italic("bla")["CMY-2"]),
                 expression(italic("bla")["CTX-M"]),
                 expression(italic("bla")["SHV"]),
                 expression(italic("bla")["TEM"]),
                 expression(italic("AAC(3)-II")),
                 expression(italic("aadA")),
                 expression(italic("aph")),
                 expression(italic("catA1")),
                 expression(italic("cmlA1")),
                 expression(italic("floR")),
                 expression(italic("sul")),
                 expression(italic("tet")),
                 expression(italic("dfrA")))


gene_plot <- ggplot(gene_plot_df, aes(reorder(factor(group), mean_antres), gene, fill = type)) +
  geom_line(aes(group = group),
            size = 0.4,
            color = "grey50",
            alpha = 0.7) +
  annotate("rect",
           ymin = c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5),
           ymax = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5),
           xmin = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
           xmax = c(Inf, Inf, Inf,Inf, Inf, Inf,Inf, Inf, Inf),
           fill = "grey80",
           alpha = 0.3) +
  geom_point(size = 3.1,
             color = "grey50",
             pch = 21) +
  geom_hline(yintercept = 3.5,
             size = 0.3) +
  guides(fill = FALSE,
         shape = FALSE) +
  labs(y = "Genes") +
  scale_fill_manual(values = type_palette) +
  scale_y_discrete(labels = gene_labels) +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = "italic"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

res_plot <- ggplot(res_plot_df, aes(reorder(factor(group), mean_antres), AM, fill = type)) +
  geom_line(aes(group = group),
            size = 0.4,
            color = "grey50",
            alpha = 0.7) +
  annotate("rect",
           ymin = c(0.5, 2.5, 4.5, 6.5, 8.5),
           ymax = c(1.5, 3.5, 5.5, 7.5, 9.5),
           xmin = c(-Inf, -Inf, -Inf, -Inf, -Inf),
           xmax = c(Inf, Inf, Inf, Inf, Inf),
           fill = "grey80",
           alpha = 0.3) +
  geom_point(color = "grey50",
             pch = 23,
             size = 3.1) +
  guides(fill = FALSE) +
  labs(y = "Antimicrobials") +
  scale_fill_manual(values = type_palette) +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



p <- plot_grid(
  bar_plot,
  gene_plot,
  res_plot,
  ncol = 1,
  align = "v",
  rel_heights = c(0.2, 0.6, 0.2)
)

ggsave("test2.png",
       p,
       device = "png",
       dpi = 600,
       height = 12,
       width = 12)


names(group_df) <- sub("_group", "", names(group_df))

library(tibble)

cor_df <- group_df %>%
  column_to_rownames("id")


test <- cor.test(cor_df$qnr, cor_df$blaCTXM)

tests <- list()

for (i in 1:length(cor_df)) {
  res <- cor.test(cor_df[,8], cor_df[,i])
}






merged_acquired_report <- read.table("merged_acquired_report.txt",
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE)

testdf <- merged_acquired_report[, colSums(merged_acquired_report != 0) > 0] %>%
  column_to_rownames("id")



gene_names <- names(testdf)

tests <- list()

for (i in gene_names) {
  cor_data <- cor.test(testdf[["qnr"]],
                       testdf[[i]],
                       alternative = "two.sided",
                       method = "pearson",
                       conf.level = 0.95,
                       exact = TRUE)
  
  df <- data.frame(id = cor_data$data.name,
             cor_val = cor_data$estimate,
             p_val = cor_data$p.value,
             lwr = cor_data$conf.int[1],
             upr = cor_data$conf.int[2])
  
  tests[[i]] <- df
}

test2 <- bind_rows(tests, .id = "ref") %>%
  filter(ref != "qnr") %>%
  mutate(sign = ifelse(p_val < 0.05, "Significant", "Non-significant"))


p <- ggplot(test2, aes(reorder(ref, -cor_val), cor_val, fill = sign)) +
  geom_col() +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                size = 0.5,
                width = 0.4) +
  theme_classic() +
  theme(axis.title = element_blank()) +
  coord_flip() +
  ggtitle("Rho-correlation to qnr-genes")

ggsave("test.png",
       p,
       device = "png",
       units = "cm",
       dpi = 600,
       height = 20,
       width = 20)


ggplot(correlations, aes(Gene, Cor, fill = AB)) +
  geom_col(position = position_dodge())


mic_genes <- mic_data %>%
  left_join(testdf, by = "id") %>%
  column_to_rownames("id")


cor.test(mic_genes$NAL, mic_genes$qepA)

cordata <- read_excel("data/correlation_AB_gene.xlsx") %>%
  mutate(sign = ifelse(`p-value` > 0.05, 0, 1)) %>%
  mutate(AB = factor(AB, levels = c("AMP","CTX","CIP","NAL","GEN","SMX","TMP","CHL","TET")))


p<- ggplot(cordata, aes(reorder(Gene, Cor), Cor, fill = AB)) +
  geom_col(position = position_dodge(0.9),
           color = "black",
           size = 0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(0.9),
                width = 0.4,
                size = 0.8) +
  scale_fill_manual(values = type_palette) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.title = element_blank()) +
  annotate("text",
           x = c(14.9,13.9, 12.9, 11.9, 10.9, 9.9, 9.1, 8.65, 7.9, 6.9, 6.1, 5.65, 5.1, 4.65, 3.9, 1.1),
           y = c(1.04, 1.02, 1, 0.97, 0.95, 0.86, 0.85, 0.36, 0.53, 0.45, 0.44, 0.28, 0.44, 0.28, 0.32, -0.93), label = "*",
           size = 8) +
  coord_flip()

type_palette <- c("CIP" = "#8dd3c7",
                  "NAL" = "#ccebc5",
                  "GEN" = "#ffffb3",
                  "AMP" = "#bebada",
                  "CTX" = "#bc80bd",
                  "CHL" = "#fb8072",
                  "SMX" = "#80b1d3",
                  "TET" = "#fdb462",
                  "TMP" = "#b3de69")




do_cor_test <- function(df, list1, list2) {
  tests <- list()
  
  for (i in list1) {
    for (j in list2) {
      cor_data <- cor.test(df[[j]],
                           df[[i]],
                           alternative = "two.sided",
                           method = "pearson",
                           conf.level = 0.95,
                           exact = FALSE)
      
      df2 <- data.frame(id = cor_data$data.name,
                        cor_val = cor_data$estimate,
                        p_val = cor_data$p.value,
                        lwr = cor_data$conf.int[1],
                        upr = cor_data$conf.int[2])
      
      tests[[i]] <- df2
    }
  }
  
  final_df <- bind_rows(tests, .id = "ref")
  
  return(final_df)
}




