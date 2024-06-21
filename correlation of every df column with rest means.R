iscor <- function(x){
  rlist <- list()
  for(i in 1:ncol(x)){
    cr <- tibble(res.cor =cor(x[,i], rowMeans(x[,-i],na.rm = T), use = "pairwise.complete.obs", method = "spearman"),Series = colnames(x)[i])
    rlist[[i]] <- cr}
  
  rlist <- do.call(rbind,rlist)
  return(rlist)
}