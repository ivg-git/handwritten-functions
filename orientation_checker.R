checker <- function(i){
  pid1 <- pid(pairwiseAlignment(pattern = as.character(dataset[names(dataset)==ref]), subject = 
                                  as.character(dataset[i])))
  pid2 <- pid(pairwiseAlignment(pattern = as.character(dataset[names(dataset)==ref]), subject = 
                                  as.character(reverseComplement(dataset[i]))))
  res <- tibble(name = names(dataset)[i], forvard = pid1, reverse = pid2) %>% mutate(
    orientation = ifelse(forvard > reverse, T, F)
  )
  return(res)
}
orientation_checker <- function(dataset, ref){
  library(msa)
  library(parabar)
  dataset <- dataset
  ref <- ref
    backend <- start_backend(cores = parallel::detectCores()-1, cluster_type = "psock", backend_type = "async")
  evaluate(backend, expression = {library(msa); library(tidyverse)})
  export(backend, variables = c("checker", "dataset", "ref"), environment = environment())
  configure_bar(type = "modern", format = "[:bar] :percent eta: :eta")
  
  result <- par_lapply(backend = backend, x = seq_along(names(dataset)), fun = checker)
result <- do.call(rbind, result)
x1 <- dataset[names(dataset) %in% result$name[result$orientation == T]]
x2 <- dataset[names(dataset) %in% result$name[result$orientation == F]]
x2 <- reverseComplement(x2)
x3 <- c(x1,x2)
return(x3)
}
