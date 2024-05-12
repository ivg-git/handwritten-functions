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