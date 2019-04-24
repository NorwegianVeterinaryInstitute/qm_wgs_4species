# Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)

# Report directories
megares_report_loc <- "../data/ariba/megares_results"
resfinder_report_loc <- "../data/ariba/resfinder_results"

# Genes of interest. These genes are partially matched to 
# the generated gene names from the "fix_gene_names" function
ac_genes <- c("qnr","oqx","qep","aac6")
mut_genes <- c("gyr","par","marR","marA",
               "sox","rpoB","rob")

# -------------------------- Functions ----------------------------

# Paste unique cell values together, separate with ","
func_paste2 <- function(x) paste(unique(x[!is.na(x)]),
                                 collapse = ",")


# Confidence interval function
get_binCI <- function(x, n) as.numeric(
  setNames(
    binom.test(x,n)$conf.int*100, 
    c("lwr", "upr")
  )
)


# Identifies filenames in input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath, 
                      pattern = "amr_report.tsv")
  return(files)
}


# Import ariba data from report.tsv identified with file_names
get_ariba_data <- function(filepath) {
  # Identify file names in filepath
  files <- file_names(filepath) 
  
  # Import data
  data_list <- lapply(files,
                      FUN = function(file) {
                        read.delim(
                          paste0(filepath, "/", file),
                          stringsAsFactors = F,
                          header = TRUE,
                          sep = "\t"
                        )
                      })
  
  # Set file names in list
  names(data_list) <- files
  
  # Bind to data frame and convert columns to character
  data <- bind_rows(
    lapply(
      data_list,
      function(x) map(x, as.character)
    ),
    .id = "ref")
  
  return(data)
}


# Corrects gene names 
fix_gene_names <- function(df) {
  # Identify unique entries in ref_name column
  genes <- unique(df$ref_name)
  
  # Strip away excess characters
  new_names <- gsub("^(.*?)\\..*", "\\1", genes)
  new_names <- gsub("_", "", new_names, fixed = T)
  new_names <- gsub("-", "", new_names, fixed = T)
  
  # make all genes lower case, with last letter 
  # upper case if four characters long
  gene_names <- c()
  
  for (i in new_names) {
    p <- paste(
      tolower(
        substring(i, 1,3)
      ),
      substring(i, 4),
      sep = "",
      collapse = " ")
    gene_names <- c(gene_names,p)
  }
  
  # Create data frame for matching names below
  df2 <- data.frame(genes,gene_names) %>%
    mutate(genes = as.character(genes)) %>%
    dplyr::rename(ref_name = genes)
  
  # Left join to match with ref_name column
  # and strip ref column for file suffix
  df <- df %>%
    left_join(df2, by = "ref_name") %>%
    mutate(gene_names = as.character(gene_names),
           ref = gsub("(.*?)_amr_report.tsv", "\\1", ref)) %>%
    mutate(gene_names = ifelse(gene_names == "soxS", "soxR", gene_names))
  
  return(df)
}


# Function for selecting genes of interest and filtering the 
# columns in the data frame on the resulting vector
select_genes <- function(df, gene_string) {
  names_df <- names(df)
  x <- "ref"
  for (i in names_df) {
    for (j in gene_string) {
      if (str_detect(i, regex(j, ignore_case = T)) == TRUE) {
        x <- c(x, i)
      }
    }
  }
  df <- df %>%
    select(one_of(x))
  return(df)
}


# Allowed flags from the ARIBA report
flag_selection <- c("19","27","147","155","403",
                    "411","915","923","787","795",
                    "531","539","659","667","787","795")  


# Function that returns info on flag selection
check_flags <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change, ref_nt, ctg_nt) %>%
    mutate(flag_result = flag %in% flag_selection,
           flag_result = as.integer(flag_result),
           nt_change = paste(ref_nt, ctg_nt, sep = "-")) %>%
    mutate(nt_change = ifelse(nt_change == ".-.", "0", nt_change)) %>%
    dplyr::rename("gene" = gene_names) %>%
    select(-c(ref_nt, ctg_nt))
  return(df)
}


# Function that handles megares data and returns a data 
# frame with information on whether a gene is mutated 
# or not
create_mut_table <- function(df) {
  df <- df %>%
    mutate(ref_ctg_change = ifelse(ref_ctg_change == ".", "", ref_ctg_change),
           ref_nt = ifelse(ref_nt == ".", "", ref_nt),
           ctg_nt = ifelse(ctg_nt == ".", "", ctg_nt),
           ref_ctg_nt_change = paste(ref_nt, ctg_nt, sep = "-"),
           mutation = paste(ref_ctg_change, ref_ctg_nt_change, sep = "/"),
           mutation = ifelse(mutation == "/-", "0", mutation)) %>%
    select(ref, gene_names, flag, mutation) %>%
    filter(flag %in% flag_selection) %>%
    mutate(id = 1:n()) %>%
    spread(gene_names, mutation) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id, flag)) %>%
    gather(gene, mut, -ref) %>%
    mutate(mut = ifelse(mut == "" | mut == "." | is.na(mut) == TRUE, 0, mut),
           result = ifelse(mut != 0, 1, 0),
           result = as.integer(result),
           type = "mut",
           nt_change = gsub("[^,]+?/([^,]+?)", "\\1", mut),
           mut = gsub("([^/]+)/[^,]+", "\\1", mut))
  return(df)
}                    


# Corrects for mutations in QRDR for gyrA, gyrB, parC and parE genes.
# The function returns columns of 1/0 values for whether or not the 
# mutations reported by ARIBA is within the QRDR in the gene:
# gyrA: AA 67 - 106
# gyrB: AA 333 - 481
# ParC: AA 51 - 170
# parE: AA 366 - 523
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1142146/pdf/cjvr68pg229.pdf
fix_gyr_par_results <- function(df) {
  df <- df %>%
    # control for mutation within QRDR for gyrA
    mutate(gyrA_result = mut %>% 
             str_extract_all("\\d+") %>% # from reprex package
             map(as.integer) %>% # converts all to integer
             # returns TRUE/FALSE whether value is within range or not
             map_lgl(~ any(.x >= 67L & .x <= 106L)), 
           gyrA_result = if_else(gene != "gyrA", NA, gyrA_result),
           gyrA_result = as.integer(gyrA_result),
           # control for mutation within QRDR for gyrB
           gyrB_result = mut %>%
             str_extract_all("\\d+") %>%
             map(as.integer) %>%
             map_lgl(~ any(.x >= 333L & .x <= 481L)),
           gyrB_result = if_else(gene != "gyrB", NA, gyrB_result),
           gyrB_result = as.integer(gyrB_result),
           # control for mutation within QRDR for parC
           parC_result = mut %>% 
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 51L & .x <= 170L)),
           parC_result = if_else(gene != "parC", NA, parC_result),
           parC_result = as.integer(parC_result),
           # control for mutation within QRDR for parE
           parE_result = mut %>% 
             str_extract_all("\\d+") %>% 
             map(as.integer) %>% 
             map_lgl(~ any(.x >= 366L & .x <= 523L)),
           parE_result = if_else(gene != "parE", NA, parE_result),
           parE_result = as.integer(parE_result)) %>%
    # mergeresults from all four genes into one column
    mutate(result_gyr_par = case_when(gene == "gyrA" ~ gyrA_result,
                                      gene == "gyrB" ~ gyrB_result,
                                      gene == "parC" ~ parC_result,
                                      gene == "parE" ~ parE_result)) %>%
    # merge results from other genes with the four genes
    mutate(result_total = ifelse(gene %in% c("gyrA","gyrB","parC","parE"),
                                 result_gyr_par, result)) %>%
    select(-c(gyrA_result,
              gyrB_result,
              parC_result,
              parE_result,
              result_gyr_par,
              result))
  return(df)
}


# Function that handles resfinder data and returns a data frame with 
# presence/absence for acquired genes                           
create_acquired_table <- function(df) {
  df <- df %>%
    select(ref, gene_names, flag, ref_ctg_change) %>%
    mutate(id = 1:n()) %>%
    # Filter flags
    filter(flag %in% flag_selection) %>%
    spread(gene_names, ref_ctg_change) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate(ref = gsub("^(.*?)_amr_report.tsv$", "\\1", ref)) %>%
    select(-c(id, flag)) %>%
    gather(gene, result,-ref) %>%
    mutate(result_total = ifelse(result == "", 0, 1),
           result_total = as.character(result_total),
           type = "gene") %>%
    select(-result)
  return(df)
}


# Function that returns a filtered dataframe based 
# on the vector "mut_genes"
filter_mut_table <- function(df) {
  # Identify unique genes in column
  report_genes <- unique(df$gene)
  grep_genes <- c()
  
  # match names from name vector to column entries
  for (gene in mut_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  
  # Filter data frame based on values in grep_genes
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = as.character(gene))
  return(df)
}


# Function that returns a filtered dataframe based
# on the string "ac_genes"
filter_acquired_table <- function(df) {
  # Identify unique genes in column
  report_genes <- unique(df$gene)
  grep_genes <- c()
  
  # match names from vector to column entries
  for (gene in ac_genes) {
    for (g in report_genes) {
      if (grepl(gene, g, ignore.case = T) == TRUE) {
        grep_genes <- c(grep_genes, g)
      }
    }
  }
  
  # Filter data frame based on values in grep_genes
  df <- df %>%
    filter(gene %in% grep_genes) %>%
    mutate(gene = gsub("_", "", gene),
           gene = as.character(gene))
  return(df)
}

# calculates percentage of present/absent mutations and genes
calc_stats <- function(df) {
  df <- df %>%
    # Count Presence/Absence of each gene
    group_by(gene, result_total) %>%
    dplyr::count() %>%
    ungroup() %>%
    # Convert binary to "Present" and "Absent"
    mutate(result_total = if_else(result_total == 1, 
                                  "Present", "Absent")) %>%
    # Spread to separate Present and Absent columns
    spread(result_total, n, fill = 0) %>%
    rowwise() %>%
    # Calculate percentages
    mutate(Total = Present + Absent,
           Percent = round(Present/Total*100, 1))
  return(df)
}

# calculates how many mutations are present in the genes
# gyrA, gyrB, parC and parE
calc_no_of_mut <- function(df) {
  df1 <- df %>%
    # Filter out all other genes
    filter(gene %in% c("gyrA","parC","gyrB","parE")) %>%
    select(-nt_change) %>%
    # Count how many mutations in each gene per row
    mutate(entries = sapply(strsplit(.$mut, ","),
                            FUN = function(x) {length(x)})) %>%
    # Separate each mutation to own column
    separate(mut, into = as.character(c(1:max(.$entries)))) %>%
    select(-entries) %>%
    # Gather all mutations together into one column,
    # to get one row per mutation
    gather(id, mut, -c(ref, gene, result, type))
  
  # Correct for QRDR mutations in the four genes
  df2 <- fix_gyr_par_results(df1)
  
  # Generate final table
  df3 <- df2 %>%
    # Filter out mutations outide QRDR
    mutate(test = ifelse(mut != "0" & result_total == 0, 0, 1)) %>%
    filter(test == 1) %>%
    # Create own column for each gene
    spread(gene, mut) %>%
    group_by(ref) %>%
    # Summarise all columns based on unique values
    summarise_all(funs(func_paste)) %>%
    # Remove "0" values for correct counts below
    mutate_at(.vars = vars(c("gyrA","gyrB","parC","parE")),
              .funs = function(x) ifelse(x == "0", "", x)) %>%
    # Count mutations in each column entry
    mutate(mut_gyrA = sapply(strsplit(.$gyrA, ","),
                             FUN = function(x) {length(x)}),
           mut_gyrB = sapply(strsplit(.$gyrB, ","),
                             FUN = function(x) {length(x)}),
           mut_parC = sapply(strsplit(.$parC, ","),
                             FUN = function(x) {length(x)}),
           mut_parE = sapply(strsplit(.$parE, ","),
                             FUN = function(x) {length(x)})) %>%
    select(-c(type, id, result_total, test)) %>%
    # Reinsert "0" values 
    mutate_at(.vars = vars(c("gyrA","gyrB","parC","parE")),
              .funs = function(x) ifelse(x == "", "0", x))
  return(df3)
}

# creates a data frame with one row per sample, and 1/0 results 
# for mutations in respective genes in columns
create_mut_report <- function(df) {
  df <- df %>%
    select(-mut) %>%
    group_by(id) %>%
    mutate(id2 = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id2, type))
  return(df)
}

# creates a data frame with one row per sample, and 1/0 results
# for acquired genes in columns
create_acquired_report <- function(df) {
  df <- df %>%
    group_by(id) %>%
    mutate(id2 = 1:n()) %>%
    spread(gene, result_total) %>%
    summarise_all(funs(func_paste)) %>%
    select(-c(id2, type))
  return(df)
}

# -------------------------- Analysis ---------------------------
# Import data
mut_data <- get_ariba_data(megares_report_loc)
acquired_data <- get_ariba_data(resfinder_report_loc)

# Clean data
clean_mut_data <- fix_gene_names(mut_data)
clean_acquired_data <- fix_gene_names(acquired_data)

# Check flags
mut_flags <- check_flags(clean_mut_data) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref) %>%
  # Filters out unwanted samples
  filter(id %not_in% ex_samples)

acquired_flags <- check_flags(clean_acquired_data) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref) %>%
  filter(id %not_in% ex_samples)

mut_table <- create_mut_table(clean_mut_data)

mut_table_fixed <- fix_gyr_par_results(mut_table) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref) %>%
  filter(id %not_in% ex_samples)

mut_filtered <- filter_mut_table(mut_table_fixed)

acquired_table <- create_acquired_table(clean_acquired_data) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref) %>%
  filter(id %not_in% ex_samples)

acquired_filtered <- filter_acquired_table(acquired_table)

mut_stats <- calc_stats(mut_filtered)
mut_report <- create_mut_report(mut_filtered)

acquired_stats <- calc_stats(acquired_filtered)
acquired_report <- create_acquired_report(acquired_filtered)

mut_quant <- suppressWarnings(calc_no_of_mut(mut_table)) %>%
  # Correct for erroneous entry from ARIBA
  mutate(parC = ifelse(parC == "S58I", "S80I", parC)) %>%
  left_join(id_data, by = "ref") %>%
  select(id, everything(), -ref) %>%
  filter(id %not_in% ex_samples)

# number of mutations table
mut_no <- mut_quant %>%
  group_by(mut_gyrA, mut_gyrB, mut_parC, mut_parE) %>%
  dplyr::count()

# number of isolates per mutation combination
mut_comb <- mut_quant %>%
  group_by(gyrA, gyrB, parC, parE) %>%
  dplyr::count()