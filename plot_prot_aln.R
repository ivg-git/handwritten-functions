plot_prot_aln <- function(AAStringSet){
  protein_names <- AAStringSet@ranges@NAMES
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(AAStringSet), "")), nrow = length(as.character(AAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  
  ggplot(df_aa) +
    geom_tile(aes(x = column, y = row, fill = aa)) +
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "AA") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values = viridis::turbo(21))
}
