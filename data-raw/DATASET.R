## code to prepare `DATASET` dataset goes here


# generate mock data
# Load required library
library(dplyr)

# Function to generate realistic CDR3 amino acid sequences
generate_cdr3_aa <- function(length = sample(12:16, 1)) {
  paste0(sample(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], length, replace = TRUE), collapse = "")
}

# Function to generate realistic nucleotide sequences
generate_cdr3_nt <- function(length = sample(36:45, 1)) {
  paste0(sample(strsplit("ATGC", "")[[1]], length, replace = TRUE), collapse = "")
}

# Function to generate a realistic V gene name
generate_v_name <- function() {
  paste0("TRBV", sample(1:30, 1))
}

# Function to generate a realistic J gene name
generate_j_name <- function() {
  paste0("TRBJ", sample(1:10, 1))
}
# Function to generate a realistic J gene name
generate_d_name <- function() {
  paste0("TRBD", sample(1:2, 1))
}


# Function to generate VDJtools formatted data
generate_vdjtools_df <- function(expanded_clonotype = NULL,
                                 expanded_count = 1000,
                                 base_count = 10,
                                 n_clonotypes = 50) {
#  total_count <- expanded_count + base_count * (n_clonotypes - 1)
  df <- data.frame(
    count = rep(base_count, n_clonotypes),
    freq = NA,
    cdr3nt = replicate(n_clonotypes, generate_cdr3_nt()),
    cdr3aa = replicate(n_clonotypes, generate_cdr3_aa()),
    v = replicate(n_clonotypes, generate_v_name()),
    d = replicate(n_clonotypes, generate_d_name()),
    j = replicate(n_clonotypes, generate_j_name()),
    VEnd	= sample(5:13, 1, size = n_clonotypes),
    DStart = NA,
    DEnd	= sample(1:5, 1, size = n_clonotypes),
    JStart = sample(16:35, 1, size = n_clonotypes),

    stringsAsFactors = FALSE
  )

  if (!is.null(expanded_clonotype)) {
    # get random index for a clone to be expanded
    i = sample(1:n_clonotypes, 1)
    df$cdr3aa[i] <- expanded_clonotype
    df$cdr3nt[i] <- generate_cdr3_nt()
    df$v[i] <- generate_v_name()
    df$j[i] <- generate_j_name()
    df$d[i] <- generate_d_name()
    df$count[i] <- expanded_count
  }

  df$freq <- round(df$count / sum(df$count), 6)
  return(df)
}

# Generate and save peptide-stimulated samples
for (i in 1:5) {
  clonotype <- generate_cdr3_aa()
  df <- generate_vdjtools_df(expanded_clonotype = clonotype)
  write.table(df, file = paste0("sample1_peptide", i, ".tsv"), sep = '\t',
              row.names = FALSE, quote = FALSE)
}

# Generate and save control sample
df_control <- generate_vdjtools_df()
write.table(df_control, file = "sample1_control.tsv",  sep = '\t',
            row.names = FALSE, quote = FALSE)



# Generate data with replicates
rep = c("A","B","C")
for (i in 1:5) {
  clonotype <- generate_cdr3_aa()
  for(j in rep){
    df <- generate_vdjtools_df(expanded_clonotype = clonotype)
    write.table(df, file = paste0("replicates/sample1_peptide", i, "_",j,".tsv"), sep = '\t',
                row.names = FALSE, quote = FALSE)
  }
}

# Generate and save reolicates of control sample
for(j in rep){
  df_control <- generate_vdjtools_df()
write.table(df_control, file = paste0("replicates/sample1_control", "_",j,".tsv"),  sep = '\t',
            row.names = FALSE, quote = FALSE)
}

#usethis::use_data(DATASET, overwrite = TRUE)
