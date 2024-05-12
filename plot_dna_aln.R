plot_DNA_aln <- function(DNAStringSet){
  library(tidyverse)
  protein_names <- DNAStringSet@ranges@NAMES
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(DNAStringSet), "")), nrow = length(as.character(DNAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  ggplot() +
    geom_tile(data = df_aa, aes(x = column, y = row, fill = aa)) +
    scale_fill_manual(values = c("G" = "palegreen", "T" = "skyblue", "A" = "gold", "-" = "white", "C" = "tomato")) +
    # geom_text(data = df_aa, aes(x = column, y = row, label = aa), color = "black", size = 3) +
    theme_minimal() +
    labs(x = "Position", y = "DNA", fill = "Structure")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}