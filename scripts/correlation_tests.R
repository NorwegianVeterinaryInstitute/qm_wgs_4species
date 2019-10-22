# Correlation test script

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)

ex_samples <- c(
  "2016-17-565-1-S",
  "2016-02-428-2-S",
  "2016-02-486-2-S",
  "2016-02-732-2-S",
  "2014-01-1741-1-S"
)

id_data <- read.table("data/id.txt",
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = F)

isolate_data <- read.table("data/isolate_data.txt",
                           sep = "\t",
                           header = T,
                           stringsAsFactors = F) %>%
  filter(!id %in% ex_samples)

acquired_data <- read.table("data/VAMPIR/amr_ac/acquired_gene_report.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref)

merged_acquired_data <- read.table("merged_acquired_report.txt",
                                   sep = ",",
                                   header = TRUE,
                                   stringsAsFactors = FALSE)

intrinsic_genes <- read.table("data/ariba/mut_report.txt",
                             sep = "\t",
                             header = TRUE,
                             stringsAsFactors = FALSE)



all_mechanisms <- acquired_data %>%
  left_join(intrinsic_genes, by = "id") %>%
  select(-id)

merged_mechanisms <- merged_acquired_data %>%
  left_join(intrinsic_genes, by = "id") %>%
  select(-id)

results <- list()
results_merged <- list()

for (i in names(all_mechanisms)) {
  results[[i]] <- apply(all_mechanisms, 2, cor.test, all_mechanisms[[i]], method="pearson", exact = FALSE)
  
}

for (i in names(merged_mechanisms)) {
  results_merged[[i]] <- apply(merged_mechanisms, 2, cor.test, merged_mechanisms[[i]], method="pearson", exact = FALSE)
  
}

library(broom)

test <- lapply(results, function(x) lapply(x, tidy))


test2 <- lapply(test, function(x) bind_rows(x, .id = "ref2"))


test3 <- bind_rows(test2, .id = "ref1") %>%
  mutate(sign = ifelse(p.value < 0.05, 1, 0)) %>%
  filter(ref1 != ref2)

results_merged_clean <- lapply(results_merged, function(x) lapply(x, tidy))
results_merged_clean2 <- lapply(results_merged_clean, function(x) bind_rows(x, .id = "ref2"))

gene_order <- c("qnr",
                "qepA",
                "blaCMY",
                "blaCTXM",
                "blaSHV",
                "blaTEM",
                "aac3II",
                "aadA",
                "aph",
                "catA",
                "cmlA",
                "floR",
                "sul",
                "tet",
                "dfrA")



results_merged_clean3 <- bind_rows(results_merged_clean2, .id = "ref1") %>%
  mutate(sign = ifelse(p.value < 0.05, "Significant", "Nonsignificant"),
         sign2 = ifelse(p.value < 0.05, "*", NA),
         dir = ifelse(estimate < 0, "Neg", "Pos")) %>%
  filter(ref1 != ref2) %>%
  filter(ref2 %in% c("gyrA", "parC", "parE", "rpoB", "soxR", "marA", "marR"),
         !ref1 %in% c("gyrA", "gyrB", "parC", 
                      "parE", "rpoB", "soxR",
                      "marA", "marR", "robA", 
                      "fosA", "blaACT", "msrE", 
                      "ermB", "mph", "lnuF")) %>%
  mutate(ref2 = paste0(toupper(substr(ref2, 1,1)),
                       substr(ref2, 2,4)),
         ref1 = factor(ref1, levels = gene_order),
         ref2 = factor(ref2, levels = c("GyrA", "ParC", "ParE", "MarA", "MarR", "SoxR", "RpoB")))


palette <- c("Pos" = "#45ADA8",
             "Neg" = "#594F4F")

gene_labels <- c(expression(italic("qnr")),
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



p <- ggplot(results_merged_clean3, aes(ref1, estimate, fill = dir)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.4) +
  geom_text(aes(label = sign2),
            vjust = 0.75,
            hjust = ifelse(results_merged_clean3$estimate < 0, 1.5, -0.5),
            size = 4,
            color = "red") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(labels = gene_labels) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        strip.background = element_rect(fill = NA, color = "black"),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 6)) +
  guides(fill = FALSE) +
  coord_flip() +
  facet_wrap(~ref2)

ggsave("figures/corr_chrom_vs_plasmid.png",
       p,
       device = "png",
       units = "cm",
       dpi = 600,
       height = 15,
       width = 15)



gyrA_data <- results_merged_clean3 %>%
  filter(ref2 == "gyrA")


p1_gyrA <- ggplot(gyrA_data, aes(ref1, estimate, fill = dir)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.4) +
  geom_text(aes(label = sign2),
            vjust = 0.75,
            hjust = ifelse(gyrA_data$estimate < 0, 1.5, -0.5),
            size = 4,
            color = "red") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = palette) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = FALSE) +
  coord_flip()


n_ac <- merged_acquired_data %>%
  select(-id) %>%
  gather(gene, value) %>%
  group_by(gene, value) %>%
  count() %>%
  ungroup() %>%
  mutate(value = ifelse(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  mutate(Total = Present + Absent,
         Percent = round(Present / Total * 100, 1)) %>%
  filter(Present != 0)

p2 <- ggplot(n_ac, aes(gene, Percent)) +
  geom_col(color = "black") +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  coord_flip()

library(cowplot)

plot_grid(p1_gyrA, p2, ncol = 2, nrow = 1, align = "h")

count(merged_acquired_data, floR)















