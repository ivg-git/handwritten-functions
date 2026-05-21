plot_prot_aln2 <- function(alignment, highlight_mismatches = TRUE, consensus_percent = 0.5){
  
  colors <- c(`-`="#000000", `A`="#BDB1E8", `R`="#EFA2C5", `N`="#F6602F",
              `D`="#FD5559", `C`="#12C7FE", `Q`="#DDACB4", `E`="#FEA097", `G`="#F46802",
              `H`="#FCA708", `I`="#369BD9", `L`="#2E95EC", `K`="#CF7690", `M`="#4B8EFE",
              `F`="#76997D", `P`="#FD2AE3", `S`="#A08A9A", `T`="#9A84D5", `W`="#74C80D",
              `Y`="#9BB896", `V`="#89B9F9")
  
  consensus <- DECIPHER::ConsensusSequence(alignment, threshold = 1 - consensus_percent)%>% 
    as.character() %>% str_split("") %>% unlist()
writeXStringSet(alignment, "temp.aln")
ph <- as.phyDat(ape::read.FASTA("temp.aln", type = "AA"))
tip_order <- get_taxa_name(ggtree(phangorn::upgma(phangorn::dist.ml(ph))))
  
  
  aln_mat <- as.matrix(alignment)
  aln_mat <- aln_mat[tip_order,]
  
  aln_long <- as.data.frame(aln_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "position", values_to = "character") %>%
    mutate(
      taxon = factor(taxon, levels = rev(tip_order)),
      position = as.numeric(gsub("V", "", position)),
      character = factor(character, levels = names(colors))
    )
  
  if(highlight_mismatches) {
    aln_long <- aln_long %>%
      mutate(
        consensus_char = consensus[position],
        is_match = character == consensus_char | character == "-",
        alpha_val = ifelse(is_match, 0.2, 1)
      )
    
    p2 <- ggplot(aln_long)+
      geom_tile(aes(x = position,
                    y = taxon,
                    fill = character,
                    alpha = alpha_val),
                height = 0.8) +
      scale_fill_manual(values = colors) +
      scale_alpha_identity() +  #using the alpha directly
      theme_light() + 
      theme(legend.position = "none")
    
  } else {
    p2 <- ggplot(aln_long)+
      geom_tile(aes(x = position,
                    y = taxon,
                    fill = character),
                height = 0.8) +
      scale_fill_manual(values = colors) +
      theme_light() + 
      theme(legend.position = "none")
  }
  return(p2)
}
