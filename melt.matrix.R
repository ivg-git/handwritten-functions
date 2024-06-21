melt_matrix <- function(x){
  res <- pivot_longer(as_tibble(as.data.frame(x) %>% rownames_to_column("Var2")), cols = 2:(ncol(x)+1), names_to = "Var1")
  return(res)
}