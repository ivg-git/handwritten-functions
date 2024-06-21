get_list_tibble <- function(ls){
  temp <- lapply(names(ls), function(name) {
    df <- ls[[name]]
    if(is.data.frame(df) && ncol(df) > 0) {
      data.frame(name = name, colnames = colnames(df), stringsAsFactors = FALSE)
    } else {
      data.frame(name = name, colnames = NA_character_, stringsAsFactors = FALSE)
    }
  }) %>% 
    bind_rows()%>% as_tibble()
  return(temp)}