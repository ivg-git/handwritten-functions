plot_HEC <- function(HEC, protein_names){
  library(tidyverse)
  matrix_data <- matrix(unlist(str_split(HEC, "")), nrow = length(HEC), byrow = TRUE)
  df <- as.data.frame(matrix_data)
  df$row <- protein_names
  df <- tidyr::pivot_longer(df, cols = -row, names_to = "column", values_to = "structure")
  df$column <- as.numeric(str_split_fixed(df$column, "V",2)[,2])
  
  ggplot(df, aes(x = column, y = row, fill = structure)) +
    geom_tile() +
    scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold")) +
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "Structure") +
    ggtitle("Secondary Structure of Proteins")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
plot_HEC2 <- function(HEC, AAStringSet){
  library(tidyverse)
  matrix_data <- matrix(unlist(str_split(HEC, "")), nrow = length(HEC), byrow = TRUE)
  df <- as.data.frame(matrix_data)
  protein_names <- AAStringSet@ranges@NAMES
  df$row <- protein_names
  df <- tidyr::pivot_longer(df, cols = -row, names_to = "column", values_to = "structure")
  df$column <- as.numeric(str_split_fixed(df$column, "V",2)[,2])
  
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(AAStringSet), "")), nrow = length(as.character(AAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  
  ggplot() +
    geom_tile(data = df, aes(x = column, y = row, fill = structure)) +
    scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold")) +
    geom_text(data = df_aa, aes(x = column, y = row, label = aa), color = "black", size = 3) +
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "Structure") +
    ggtitle("Secondary Structure of Proteins with Amino Acids")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlim(0,max(nchar(HEC))+1)}

plot_HEC3 <- function(AAStringSet){
  library(DECIPHER)
  library(tidyverse)
  HEC <- PredictHEC(AAStringSet, type = "states")
  matrix_data <- matrix(unlist(str_split(HEC, "")), nrow = length(HEC), byrow = TRUE)
  df <- as.data.frame(matrix_data)
  protein_names <- AAStringSet@ranges@NAMES
  df$row <- protein_names
  df <- tidyr::pivot_longer(df, cols = -row, names_to = "column", values_to = "structure")
  df$column <- as.numeric(str_split_fixed(df$column, "V",2)[,2])
  
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(AAStringSet), "")), nrow = length(as.character(AAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  
  ggplot() +
    geom_tile(data = df, aes(x = column, y = row, fill = structure)) +
    scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold")) +
    geom_text(data = df_aa, aes(x = column, y = row, label = aa), color = "black", size = 3) +
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "Structure") +
    ggtitle("Secondary Structure of Proteins with Amino Acids")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlim(0,max(nchar(HEC))+1)}

plot_protein_features <- function(x){
  df <- data.frame(name = str_split_fixed(x@ranges@NAMES, " ",2)[,1], seq = as.character(x))
  df <- ampir::calculate_features(df)
  df <- df[,c(1:6)]
  df <- pivot_longer(df, cols = 2:6, names_to = "parameter", values_to = "value")
  ggplot(df)+
    geom_line(aes(x = seq_name, y = value, col = parameter, group = 1), linewidth = 1)+
    facet_wrap(~parameter, scales = "free_y", ncol = 5)+theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
}

plot_HEC4 <- function(AAStringSet){
  library(DECIPHER)
  library(tidyverse)
  HEC <- PredictHEC(AAStringSet, type = "states")
  gapPositions <- start(vmatchPattern("-", AAStringSet))
  for (i in seq_along(HEC)) {
    currentGapPositions <- gapPositions[[i]]
    for (pos in currentGapPositions) {
      beforeGap <- substring(HEC[[i]], 1, pos - 1)
      afterGap <- substring(HEC[[i]], pos, nchar(HEC[[i]]))
      HEC[[i]] <- paste0(beforeGap, "-", afterGap)
    }
  }
  matrix_data <- matrix(unlist(str_split(HEC, "")), nrow = length(HEC), byrow = TRUE)
  df <- as.data.frame(matrix_data)
  protein_names <- AAStringSet@ranges@NAMES
  df$row <- protein_names
  df <- tidyr::pivot_longer(df, cols = -row, names_to = "column", values_to = "structure")
  df$column <- as.numeric(str_split_fixed(df$column, "V",2)[,2])
  
  df_aa <- as.data.frame(matrix(unlist(str_split(as.character(AAStringSet), "")), nrow = length(as.character(AAStringSet)), byrow = TRUE))
  df_aa$row <- protein_names
  df_aa <- tidyr::pivot_longer(df_aa, cols = -row, names_to = "column", values_to = "aa")
  df_aa$column <- as.numeric(str_split_fixed(df_aa$column, "V",2)[,2])
  
  ggplot() +
    geom_tile(data = df, aes(x = column, y = row, fill = structure)) +
    scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold", "-" = "white")) +
    geom_text(data = df_aa, aes(x = column, y = row, label = aa), color = "black", size = 3) +
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "Structure") +
    ggtitle("Secondary Structure of Proteins with Amino Acids")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlim(0,max(nchar(HEC))+1)}

