# Testing import of gyrB sequences from .gz files created by ARIBA
library(Biostrings)
library(msa)
library(DECIPHER)

# File locations
path <- "data/ariba/megares_full/"

## Functions
# Identifies filenames in input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath,
                      pattern = "assembled_genes.fa.gz",
                      recursive = TRUE)
  return(files)
}

# Import data from assembled_genes.fa.gz files created by ariba
get_seq_data <- function(filepath) {
  files <- file_names(filepath)
  
  data_list <- lapply(files,
                      FUN = function(file) {
                        readDNAStringSet(paste0(filepath, "/", file))
                      })
  names(data_list) <- files
  return(data_list)
}

# Filter out specific gene of interest
filter_gyrB <- function(fasta_file, gene) {
  gyrB_pos <- which(grepl(gene, names(fasta_file), ignore.case = TRUE) == TRUE)
  gyrB_seq <- fasta_file[gyrB_pos]
  return(gyrB_seq)
}

# Append sequences to fasta file
append_fasta_file <- function(list) {
  i <- 1
  for (x in list) {
    write.fasta(x, names(list[i]), "data/sequences.fasta", open = "a")
    i <- i + 1
  }
}

## Run functions
sequences <- get_seq_data(path)
gyrB_sequences <- lapply(sequences, function(x) filter_gyrB(x, "gyrb"))
gyrB_ref <- readDNAStringSet("data/gyrB_ref.fasta")
write.fasta(gyrB_ref, "gyrB_ref", "data/sequences.fasta", open = "w")
append_fasta_file(gyrB_sequences)

sequences_gyrB <- readDNAStringSet("data/sequences.fasta")
trans_seqs <- Biostrings::translate(sequences_gyrB)
alignment <- AlignSeqs(trans_seqs, iterations = 0, refinements = 0)
writeXStringSet(alignment, "data/alignment.fasta")