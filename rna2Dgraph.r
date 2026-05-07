library(tidyverse)
library(igraph)

ct2coord_graph <- function(ct) {
  
  ct <- ct %>% arrange(pos2)
  n <- nrow(ct)
  
  edges <- data.frame(from = character(), to = character())
  
  for(i in 1:(n-1)) {
    edges <- rbind(edges, data.frame(from = as.character(ct$pos2[i]), 
                                     to = as.character(ct$pos2[i+1])))
  }
  
  for(i in 1:n) {
    if(ct$bound[i] > ct$pos2[i]) {
      edges <- rbind(edges, data.frame(from = as.character(ct$pos2[i]), 
                                       to = as.character(ct$bound[i])))
    }
  }
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  set.seed(42)
  layout <- layout_with_kk(g)
  
  layout_df <- data.frame(
    pos2 = as.numeric(V(g)$name),
    x = layout[, 1],
    y = layout[, 2]
  )
  
  result <- ct %>%
    left_join(layout_df, by = "pos2") %>%
    mutate(group = 1)
  
  return(result)
}

plotRNA2D_graph <- function(input_rna){
  dbn <- PredictDBN(input_rna, type = "states", processors = NULL, pseudoknots = 0, useFreeEnergy = T)
  ct <- RRNA::makeCt(dbn,as.character(input_rna))
  dat <- ct2coord_graph(ct)
  covalent_bonds <- dat %>%
    arrange(pos2) %>%
    mutate(next_x = lead(x), #lead shifts a vector by 1
           next_y = lead(y),
           next_pos2 = lead(pos2)) %>%
    filter(!is.na(next_pos2))
  wc_bonds <- dat %>%
    filter(bound > 0, pos2 < bound) %>%  
    dplyr::rename("x_from" = "x", "y_from" = "y", "seq_from" = "seq", "pos2_from" = "pos2") %>%
    left_join(dat %>% dplyr::rename("x_to" = "x", "y_to" = "y", "seq_to" = "seq", "pos2_to" = "pos2"),
              by = c("bound" = "pos2_to")) %>%
    select("pos2" = "pos2_from", "x_from", "y_from", "bound", "x_to", "y_to", "seq_from", "seq_to")
  
  p <-  ggplot(dat) +
    geom_segment(data = covalent_bonds,
                 aes(x = x, y = y, xend = next_x, yend = next_y),
                 color = "dodgerblue", linewidth = 0.8, alpha = 0.6) +
    geom_segment(data = wc_bonds,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 color = "red", linewidth = 0.5, alpha = 0.6) +
    geom_text(aes(x = x, y = y, label = seq), size = 3) +
    theme_void() +
    coord_fixed() 
  return(p)
}


plotRNA2D_graph(RNAStringSet(nad6cds)[1])
