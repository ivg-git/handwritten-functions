plot_DNA_aln2 <- function(DNAStringSet, highlight_consensus = FALSE, consensus_percent = 0.5){
  library(tidyverse)
  library(DECIPHER)
  
  protein_names <- DNAStringSet@ranges@NAMES
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(DNAStringSet), "")), 
                                nrow = length(as.character(DNAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  
  if(highlight_consensus){
    cons <- ConsensusSequence(DNAStringSet, 
                              threshold = consensus_percent 
                              )[[1]]
    cons_chars <- str_split(as.character(cons), "")[[1]]
    
    df_aa$consensus <- cons_chars[df_aa$column]
    df_aa$matches_consensus <- df_aa$aa == df_aa$consensus & df_aa$aa != "-" & df_aa$consensus != "-"
    
    p <- ggplot() +
      geom_tile(data = df_aa, aes(x = column, y = row, fill = aa, alpha = matches_consensus)) +
      scale_fill_manual(values = c("G" = "palegreen", "T" = "skyblue", "A" = "gold", "-" = "white", "C" = "tomato")) +
      scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1), guide = "none") +
      theme_minimal() +
      labs(x = "Position", y = "DNA", fill = "Structure") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  } else {
    p <- ggplot() +
      geom_tile(data = df_aa, aes(x = column, y = row, fill = aa)) +
      scale_fill_manual(values = c("G" = "palegreen", "T" = "skyblue", "A" = "gold", "-" = "white", "C" = "tomato")) +
      theme_minimal() +
      labs(x = "Position", y = "DNA", fill = "Structure") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  return(p)
}
