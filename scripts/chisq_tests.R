# Chi-squared tests

# function for doing chi-squared tests on multiple entries in a list
do_chisq_test <- function(df) {
  df <- df %>%
    # create columns with all possible comparisons
    dplyr::slice(expand.grid(1:length(unique(species)),
                             1:length(unique(species))) %>%
                   rev %>% 
                   dplyr::filter(Var2 < Var1) %>% t) %>%
    # create test ID
    mutate(test = rep(1:(n()/2), each = 2)) %>%
    group_by(test) %>%
    # specify values for chi-squred tests
    dplyr::do(data_frame(gene = .$key,
                         test = dplyr::first(.$test),
                         species1 = dplyr::first(.$species),
                         species2 = dplyr::last(.$species),
                         data = list(matrix(c(.$Absent,
                                              .$Present),
                                            ncol = 2)))) %>%
    # do chi-squared tests
    mutate(chi_test = map(data, chisq.test,
                          correct = FALSE)) %>%
    mutate(p.value = round(map_dbl(chi_test, "p.value"),3),
           df = map_dbl(chi_test, "parameter"),
           x_sq = round(map_dbl(chi_test, "statistic"),2),
           dupl = duplicated(test)) %>%
    # filter out duplicated tests
    dplyr::filter(dupl == FALSE) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(n = sum(unlist(data))) %>%
    select(gene, species1, species2, n, df, x_sq, p.value) %>%
    # round off p-value column
    mutate(p.value = round(p.value, 5)) %>%
    dplyr::filter(p.value <= 0.05)
  return(df)
}

# create total mechanism table
total_df <- acquired_report %>%
  left_join(mut_report, by = "id") %>%
  left_join(isolate_data[c("id", "species")], by = "id") %>%
  select(-nt_change)

# create list of genes with total presence/absence values
# from each animal species
mechanism_per_species <- total_df %>%
  select(id, species, everything()) %>%
  select(-id) %>%
  gather(key, value, -species) %>%
  group_by_all() %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(value = if_else(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  select(key, species, Present, Absent) %>%
  split(f = .$key)

# Get intrinsic gene names
intrinsic_names <- names(mut_report)
intrinsic_names <- intrinsic_names[-c(1,2)]

# create percentage table over intrinsic genes
intrinsic_gene_perc <- mut_report %>%
  left_join(isolate_data, by = "id") %>%
  gather(key, value, intrinsic_names) %>%
  group_by(key, value, species) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(value = if_else(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  rowwise() %>%
  mutate(Total = Present + Absent,
         Percent = round(Present / Total * 100, 1)) %>%
  select(key, species, Percent) %>%
  dplyr::rename("gene" = key)

# Get acquired gene names
acquired_names <- names(acquired_report)
acquired_names <- acquired_names[-1]

# create percentage table over acquired genes
acquired_genes_perc <- acquired_report %>%
  left_join(isolate_data, by = "id") %>%
  gather(key, value, acquired_names) %>%
  group_by(key, value, species) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(value = if_else(value == 1, "Present", "Absent")) %>%
  spread(value, n, fill = 0) %>%
  rowwise() %>%
  mutate(Total = Present + Absent,
         Percent = round(Present / Total * 100, 1)) %>%
  select(key, species, Percent) %>%
  dplyr::rename("gene" = key)

# create total percentage table
total_perc <- rbind(intrinsic_gene_perc, acquired_genes_perc)

# run chi-squared tests
p_values_total <- lapply(mechanism_per_species, function(x) do_chisq_test(x)) %>%
  # bind results
  bind_rows() %>%
  # specify percentages for species 1
  left_join(total_perc, by = c("gene", c("species1" = "species"))) %>%
  dplyr::rename("Percent species 1" = Percent) %>%
  # specify percentages for species 2
  left_join(total_perc, by = c("gene", c("species2" = "species"))) %>%
  dplyr::rename("Percent species 2" = Percent) %>%
  dplyr::select(gene,
                species1,
                `Percent species 1`,
                species2,
                `Percent species 2`,
                n,
                df,
                x_sq,
                p.value)

p_values_total %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover"))