normal_st_read <- function(path){
  layers <- st_layers(path)[[1]]
  res <- list()
  for(i in seq_along(layers)){
    res[[i]] <- st_read(dsn = path, layer = layers[i])
    res[[i]]$layer <- layers[i]
  }
  res <- do.call(rbind,res)
  return(res)
}
