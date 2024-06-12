calculate_sirc_mst_stats <- function(x){
  eg <- ego(x, nodes = V(x), order = 1)
  drws <- lapply(eg, FUN = function(x){as.numeric(x$DRw)})
  spas <- lapply(eg, FUN = function(x){as.numeric(x$spacers)})
  lens <- lapply(eg, FUN = function(x){as.numeric(x$len)})
  meandiff <- function(x){
    for(i in seq(2,length(x))){
      proxy[[i-1]] <- x[i]-x[1]
    }
    rnn <- Reduce(mean,proxy)
    return(rnn)
  }
  
  mst.stats.fullseq <- tibble(DRw.diff = abs(sapply(drws, FUN = meandiff)),
                              spacers.diff =abs(sapply(spas, FUN = meandiff)),
                              len.diff = abs(sapply(lens, FUN = meandiff)))
  stats.fs2 <- mst.stats.fullseq %>% mutate(EN = c(1:nrow(mst.stats.fullseq))) %>% pivot_longer(cols = 1:3, names_to = "Parameter", values_to = "Value")
  return(stats.fs2)
}