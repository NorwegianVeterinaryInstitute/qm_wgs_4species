library(ape)
library(Biostrings)
library(dplyr)
library(seqRFLP)

gyrA <- readDNAStringSet("data/Sequences/DNA_gyrase_GyrA_nt_sequences.fa")
gyrB <- readDNAStringSet("data/Sequences/DNA_gyrase_GyrB_nt_sequences.fa")
parC <- readDNAStringSet("data/Sequences/DNA_topoisomerase_IV_subunit_A_nt_sequences.fa")
parE <- readDNAStringSet("data/Sequences/DNA_topoisomerase_IV_subunit_B_nt_sequences.fa")
marA <- readDNAStringSet("data/Sequences/multiple_antibiotic_resistance_transcriptional_regulator_marA_nt_sequences.fa")
marR <- readDNAStringSet("data/Sequences/marR_nt_sequences.fa")
soxR <- readDNAStringSet("data/Sequences/soxR_nt_sequences.fa")
rpoB <- readDNAStringSet("data/Sequences/RNA_polymerase_beta_subunit_nt_sequences.fa")

gyrA_names <- names(gyrA)
gyrA_seq <- paste(gyrA)

gyrB_names <- names(gyrB)
gyrB_seq <- paste(gyrB)

parC_names <- names(parC)
parC_seq <- paste(parC)

parE_names <- names(parE)
parE_seq <- paste(parE)

marR_names <- names(marR)
marR_names <- sub("\\+.*","", marR_names)
marR_seq <- paste(marR)

marA_names <- names(marA)
marA_seq <- paste(marA)

soxR_names <- names(soxR)
soxR_names <- sub("\\+.*","", soxR_names)
soxR_seq <- paste(soxR)

rpoB_names <- names(rpoB)
rpoB_seq <- paste(rpoB)

df_gyrA <- data.frame(gyrA_names, gyrA_seq)
df_gyrB <- data.frame(gyrB_names, gyrB_seq)
df_parC <- data.frame(parC_names, parC_seq)
df_parE <- data.frame(parE_names, parE_seq)
df_marR <- data.frame(marR_names, marR_seq)
df_marA <- data.frame(marA_names, marA_seq)
df_soxR <- data.frame(soxR_names, soxR_seq)
df_rpoB <- data.frame(rpoB_names, rpoB_seq)

total_data <- df_gyrA %>%
  rename("names" = gyrA_names) %>%
  left_join(df_gyrB, by = c("names" = "gyrB_names")) %>%
  left_join(df_parC, by = c("names" = "parC_names")) %>%
  left_join(df_parE, by = c("names" = "parE_names")) %>%
  #left_join(df_marR, by = c("names" = "marR_names")) %>%
  left_join(df_marA, by = c("names" = "marA_names")) %>%
  left_join(df_soxR, by = c("names" = "soxR_names")) %>%
  left_join(df_rpoB, by = c("names" = "rpoB_names")) %>%
  mutate_at(vars(contains("_seq")),
            list(as.character)) %>%
  mutate(total_seq = paste0(gyrA_seq, gyrB_seq, parC_seq, parE_seq, marA_seq, soxR_seq, rpoB_seq)) %>%
  select(names, total_seq)


fasta <- dataframe2fas(total_data)

write.fasta(fasta, file = "test.fasta")

