# Libraries
library(dplyr)
library(tidyr)
library(Biostrings)

# Functions

## Mismatch analysis function
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

## Function for identifying mutations in sequence relative to reference
run_mutation_analysis <- function(path, output = "data frame") {
  # identify files in input path
  files <- list.files(path, pattern = ".fa")
  
  # Run function across all files identified
  total_results <- c()
  
  for (file in files) {
    seq_file <- readAAStringSet(paste(path, file, sep = "/"))
    
    results <- bind_rows(lapply(seq_along(seq_file[-1]),
                                function(i) mismatches(seq_file[i + 1],
                                                       seq_file[1])))
    
    total_results[[file]] <- results
  }
  
  # output results as list
  if (output == "list") {
    return(total_results)
  }
  
  # output results as a data frame
  if (output == "data frame") {
    df <- bind_rows(total_results, .id = "ref") %>%
      mutate(ref = sub(".fa", "", ref))
    return(df)
  }
  
}

# Run functions
path <- "../data/Prokka"

test <- run_mutation_analysis(path)

test2 <- test %>%
  mutate(ref = sub("_[0-9]$", "", ref),
         mut = paste0(Reference_AA, Pos, Sample_AA),
         ID = sub(".faa", "", ID)) %>%
  left_join(id_data, by = c("ID" = "ref")) %>%
  select(ref, id, mut) %>%
  mutate(n = 1:n()) %>%
  spread(ref, mut) %>%
  group_by(id) %>%
  summarise_all(funs(func_paste)) %>%
  select(-n) %>%
  na_if("")

test3 <- id_data %>%
  left_join(test2, by = "id") %>%
  select(-ref)


test4 <- mut_filtered %>%
  select(id, gene, mut) %>%
  spread(gene, mut)

complete_test <- test3 %>%
  left_join(test4, by = "id") %>%
  select(id, GyrA, gyrA, ParC, parC, ParE, parE, MarA, marA, MarR, marR, RpoB, rpoB, SoxR, soxR) %>%
  na_if("0") %>%
  filter(id %not_in% ex_samples) %>%
  mutate(MarR = sub("V1M, ", "", MarR),
         MarR = ifelse(MarR == "V1M", NA, MarR))
