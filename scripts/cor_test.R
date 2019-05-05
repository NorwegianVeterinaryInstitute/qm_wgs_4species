AddCol <- function(df, col_name) {
  # split rows by delimeters
  string_to_proc <- df %>% select(!!col_name) %>%
    base::unlist() %>% str_split(regex("\\, |\\,")) 
  # find unique entries
  unique_strings <- string_to_proc %>%
    base::unlist() %>% base::unique()
  # construct names of the new columns
  cols_names <- paste(col_name, unique_strings, sep = "_")
  # construct 0/1-content columns for each unique entry
  cols_content <- sapply(function(i) {
    as.integer(base::unlist(lapply(function(Z) any(Z %in% unique_strings[i]), 
                                   X = string_to_proc)))
  }, X = seq_along(unique_strings))
  res <- data.frame(cols_content)
  names(res) <- cols_names
  return(res)
}

col_proc <- c("gyrA","gyrB","parC","parE")

tot_mut_report <- sapply(function(i) {AddCol(df = mut_quant, col_name = col_proc[i])},
                         X = seq_along(col_proc)) %>%
  bind_cols() %>%
  select(-c(gyrA_0, parC_0, gyrB_0, parE_0)) %>%
  mutate_all(funs(as.numeric))


test <- mut_quant %>%
  cbind(tot_mut_report) %>%
  select(-contains("mut")) %>%
  select(-c(gyrA, gyrB, parC, parE)) %>%
  left_join(acquired_report, by = "id") %>%
  mutate_at(vars(-id),
            funs(as.numeric)) %>%
  mutate(qnr = ifelse(rowSums(.[,21:27]) == 0, 0, 1)) %>%
  left_join(mut_report[,c("id","gyrA")], by = "id") %>%
  mutate(gyrA = as.numeric(gyrA)) %>%
  select(-id)

cor.test(test$gyrA_S83L, test$qnr)
cor.test(test$gyrA, test$qnr)