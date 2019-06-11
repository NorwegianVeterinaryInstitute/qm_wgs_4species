## ENA submission file creation script

library(dplyr)
library(tidyr)
library(tibble)

func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ",")

isolate_data <- read.table("data/isolate_data.txt",
                           sep = "\t",
                           header = T,
                           stringsAsFactors = F)

manifest_file <- read.table("data/ENA_submission_files/name_files/filenames_hiseq_125bp.txt",
                           sep = "\t",
                           header = F,
                           stringsAsFactors = F) %>%
  mutate(STUDY = "QRECMAP",
         NAME = "TODO",
         INSTRUMENT = "Illumina HiSeq 2000",
         LIBRARY_NAME = "Nextera FLEX",
         LIBRARY_SOURCE = "GENOMIC",
         LIBRARY_SELECTION = "RANDOM",
         DESCRIPTION = "QREC from four different animal species") %>%
  mutate(ref = sub("_L00._R._001.fastq.gz", "", V1),
         SAMPLE = ref) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  separate(V1,
           into = c("FASTQ1",
                    "FASTQ2"),
           sep = ",") %>%
  select(ref,
         STUDY,
         SAMPLE,
         NAME,
         INSTRUMENT,
         LIBRARY_NAME,
         LIBRARY_SOURCE,
         LIBRARY_SELECTION,
         FASTQ1,
         FASTQ2) %>%
  split(f = .$ref)


fix_samples <- function(df) {
  df <- df %>%
    select(-ref) %>%
    gather(field_name, field_value) %>%
    mutate(field_name = ifelse(grepl("fastq",
                                     field_name,
                                     ignore.case = T) == T,
                               "FASTQ",
                               field_name))
  return(df)
}


file_list <- lapply(manifest_file, function(x) fix_samples(x))


for (file in 1:length(file_list)) {
  write.table(file_list[[file]],
              paste0("data/ENA_submission_files/manifest_files/hiseq_xt_125bp/", names(file_list)[file], ".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}



