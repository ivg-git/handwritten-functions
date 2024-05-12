complex_graph <- function(path_to_genome, path_to_gff){
  
  Cmvd <- readDNAStringSet(path_to_genome)
  com1 <- sequence_complexity(Cmvd, window.size = 100, method = "WoottonFederhen", return.granges = T)
  com2 <- sequence_complexity(Cmvd, window.size = 100, method = "Trifonov", return.granges = T)
  com3 <- sequence_complexity(Cmvd, window.size = 100, method = "DUST", return.granges = T)
  
  calc_skew <- function(window_size, dna){
    
    calculate_gc_skew <- function(window_size, dna) {
      library(tidyverse)
      num_windows <- dna@ranges@width - window_size + 1
      gc_skew <- numeric(num_windows)
      for (i in seq(1, dna@ranges@width - window_size + 1, by = window_size)) {
        window <- subseq(dna, i, i + window_size - 1)
        c_count <- vcountPattern("C", window)
        g_count <- vcountPattern("G", window)
        gc_skew[i] <- (c_count - g_count) / (c_count + g_count)
      }
      gc_skew <- as_tibble(gc_skew)
      gc_skew$len <- c(1:nrow(gc_skew))
      gc_skew <- gc_skew[seq(1, dna@ranges@width - window_size + 1, by = window_size),]
      return(gc_skew)
    }
    
    calculate_at_skew <- function(window_size, dna) {
      library(tidyverse)
      num_windows <- dna@ranges@width - window_size + 1
      gc_skew <- numeric(num_windows)
      for (i in seq(1, dna@ranges@width - window_size + 1, by = window_size)) {
        window <- subseq(dna, i, i + window_size - 1)
        c_count <- vcountPattern("A", window)
        g_count <- vcountPattern("T", window)
        gc_skew[i] <- (c_count - g_count) / (c_count + g_count)
      }
      gc_skew <- as_tibble(gc_skew)
      gc_skew$len <- c(1:nrow(gc_skew))
      gc_skew <- gc_skew[seq(1, dna@ranges@width - window_size + 1, by = window_size),]
      return(gc_skew)
    }
    at <- calculate_at_skew(window_size, dna)
    at$skew <- "AT-skew"
    gc <- calculate_gc_skew(window_size, dna)
    gc$skew <- "GC-skew"
    skew <- rbind(at,gc)
    return(skew)}
  
  library(ggbio)
  p1 <- autoplot(com1, aes(y = complexity), geom = "line", size = 1, col = "blue", alpha = 0.8)+ylab("")
  p2 <- autoplot(com2, aes(y = complexity), geom = "line", size = 1, col = "blue", alpha = 0.8)+ylab("")
  p3 <- autoplot(com3, aes(y = complexity), geom = "line", size = 1, col = "blue", alpha = 0.8)+ylab("")
  
  skew <- calc_skew(10, Cmvd)
  
  ps <- ggplot(skew)+
    geom_line(aes(x = len, y = value, col = skew, group = skew), alpha = 0.4 )+
    stat_smooth(aes(x = len, y = value, col = skew, group = skew), method = "loess", se = F, span = 0.05)+
    geom_hline(yintercept = 0, col = "black", linetype = "dashed", linewidth = 1)+
    scale_color_manual(values = c("AT-skew" = "dodgerblue", "GC-skew" = "tomato"))+ylab("")
  
  GCcon <- letterFrequencyInSlidingView(DNAString(as.character(Cmvd)), letters = "GC", view.width = 10)
  GCcon <- as_tibble(GCcon)
  GCcon$len <- c(1:nrow(GCcon))
  colnames(GCcon)[1] <- "value"
  GCcon$value <- GCcon$value*10
  
  pgc <- ggplot(GCcon)+
    geom_line(aes(x = len, y = value), alpha = 0.4, col = "tomato" )+
    stat_smooth(aes(x = len, y = value), method = "loess", se = F, span = 0.05, col = "tomato")+
    geom_hline(yintercept = 50, col = "black", linetype = "dashed")+
    ylab("")
  
  cmv <- read_gff3(path_to_gff)
  p4 <- autoplot(cmv[cmv$type %in% c("CDS")], geom = "arrowrect", col = "black", aes(fill = product))
  #pdf 7 x 10  
  tracks(GC = pgc,Skew = ps,WF = p1, Trif = p2, DUST = p3, PCG = p4, 
         main = "Sequence complexity, GC-content and skew", heights = c(1,1,1,1,1,2)) + theme_bw()
}