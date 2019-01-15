library(stringr)
library(seqinr)
library(tidyr)
library(Biostrings)

gyrA <- readAAStringSet("data/Prokka/DNA_gyrase_GYRA_aa_sequences.fa")
parC <- readAAStringSet("data/Prokka/Topoisomerase_IV_subunit_A_aa_sequences.fa")

mismatches <- function(query, ref) {
  pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM50",
                    gapOpening = 3, gapExtension = 1) %>%
    mismatchTable() %>%
    mutate(ID=names(query),
           Pos=PatternStart,
           Reference_AA=as.character(PatternSubstring),
           Sample_AA=as.character(SubjectSubstring)) %>%
    select(ID, Reference_AA, Sample_AA, Pos) %>%
    filter(Sample_AA != "X")
}  

paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

ariba_flags_gyra <- mut_flags %>%
  select(-flag) %>%
  filter(gene == "gyrA") %>%
  mutate(gyrA_result = ref_ctg_change %>%
           str_extract_all("\\d+") %>% 
           map(as.integer) %>%
           map_lgl(~ any(.x >= 67L & .x <= 106L))) %>%
  mutate(type = case_when(gyrA_result == TRUE ~ "QRDR",
                          gyrA_result == FALSE & ref_ctg_change != "." ~ "non_QRDR",
                          gyrA_result == FALSE & ref_ctg_change == "." ~ "wild type")) %>%
  select(-gyrA_result) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup()

mut <- mut_quant %>%
  select(id, gyrA) %>%
  mutate(entries = sapply(strsplit(.$gyrA, ","), FUN = function(x) {length(x)})) %>%
  separate(gyrA, into = as.character(c(1:max(.$entries)))) %>%
  select(-entries) %>%
  gather(key, value, `1`,`2`) %>%
  filter(is.na(value) == FALSE) %>%
  select(-key) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup()

results <- bind_rows(lapply(seq_along(gyrA[-1]), function(i) mismatches(gyrA[i + 1], gyrA[1]))) %>%
  mutate(mutation = paste0(Reference_AA, Pos, Sample_AA)) %>%
  select(ID, mutation) %>%
  mutate(ID = sub(".faa", "", ID)) %>%
  dplyr::rename("ref" = ID,
                "gyrA_prokka" = mutation) %>%
  left_join(id_data, by = "ref") %>%
  select(id, gyrA_prokka) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup() %>%
  left_join(mut, by = c("id","test")) %>%
  left_join(ariba_flags_gyra, by = c("id","test")) %>%
  select(-c(test,ref_ctg_change,gene)) %>%
  dplyr::rename("gyrA_ariba" = value) %>%
  filter(id %not_in% c("2016-17-565-1-S",
                       "2016-02-428-2-S",
                       "2016-02-486-2-S",
                       "2016-02-732-2-S",
                       "2014-01-1741-1-S"))

non_seq <- id_data %>%
  mutate(test = id %in% results$id) %>%
  filter(test == FALSE) %>%
  filter(id %not_in% c("2016-17-565-1-S",
                       "2016-02-428-2-S",
                       "2016-02-486-2-S",
                       "2016-02-732-2-S",
                       "2014-01-1741-1-S")) %>%
  select(id,test) %>%
  left_join(ariba_flags_gyra, by = "id") %>%
  filter(type != "non_QRDR") %>%
  select(-c(test.x, test.y, gene)) %>%
  dplyr::rename("gyrA_ariba" = ref_ctg_change) %>%
  mutate(gyrA_prokka = NA) %>%
  select(id, gyrA_prokka, gyrA_ariba, flag_result, type)

total_results <- rbind(results, non_seq)

ids_gyrA <- unique(results$id)

test <- mut_report %>%
  filter(gyrA == 1)

test$id %in% ids_gyrA

# parC

ariba_flags_parc <- mut_flags %>%
  select(-flag) %>%
  filter(gene == "parC") %>%
  mutate(parC_result = ref_ctg_change %>%
           str_extract_all("\\d+") %>% 
           map(as.integer) %>%
           map_lgl(~ any(.x >= 51L & .x <= 170L))) %>%
  mutate(type = case_when(parC_result == TRUE ~ "QRDR",
                          parC_result == FALSE & ref_ctg_change != "." ~ "non_QRDR",
                          parC_result == FALSE & ref_ctg_change == "." ~ "wild type")) %>%
  select(-parC_result) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup()

mut_parC <- mut_quant %>%
  select(id, parC) %>%
  mutate(entries = sapply(strsplit(.$parC, ","), FUN = function(x) {length(x)})) %>%
  separate(parC, into = as.character(c(1:max(.$entries)))) %>%
  select(-entries) %>%
  gather(key, value, `1`,`2`) %>%
  filter(is.na(value) == FALSE) %>%
  select(-key) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup()

results_parC <- bind_rows(lapply(seq_along(parC[-1]), function(i) mismatches(parC[i + 1], parC[1]))) %>%
  mutate(Pos = Pos + 22) %>%
  mutate(mutation = paste0(Reference_AA, Pos, Sample_AA)) %>%
  select(ID, mutation) %>%
  mutate(ID = sub(".faa", "", ID)) %>%
  dplyr::rename("ref" = ID,
                "parC_prokka" = mutation) %>%
  left_join(id_data, by = "ref") %>%
  select(id, parC_prokka) %>%
  group_by(id) %>%
  mutate(test = 1:n()) %>%
  ungroup() %>%
  left_join(mut_parC, by = c("id","test")) %>%
  left_join(ariba_flags_parc, by = c("id","test")) %>%
  select(-c(test,ref_ctg_change,gene)) %>%
  dplyr::rename("parC_ariba" = value) %>%
  filter(id %not_in% c("2016-17-565-1-S",
                       "2016-02-428-2-S",
                       "2016-02-486-2-S",
                       "2016-02-732-2-S",
                       "2014-01-1741-1-S")) %>%
  filter(type == "QRDR")



ids <- unique(results_parC$id)

test$id %in% ids

test <- mut_report %>%
  filter(parC == 1)
