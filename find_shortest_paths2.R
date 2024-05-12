find_shortest_paths2 <- function(graph, line, togo, weight_){
  to <- unlist(str_split(ego$geneID[ego$ID == togo & ego$line == line], "/"))
  mst.probes <- dme2$ProbeName[which(dme2[[line]] !=0)]
  mst.tair <- mods$sseqid[mods$qseqid %in% mst.probes]
  test <- all_shortest_paths(graph, from = which(V(graph)$name %in% mst.tair),
                             to = which(V(graph)$name %in% to),mode = "all", weights = weight_)
  tl <- list()
  for(i in seq_along(test$res)){
    tl[[i]] <- make_empty_graph()+ test$res[[i]]$name +igraph::path(test$res[[i]]$name)  
  }
  sp <- do.call(igraph::union, tl)
  V(sp)$type <- "path"
  V(sp)$type[which(V(sp)$name %in% mst.tair)] <- "root"
  V(sp)$type[which(V(sp)$name %in% to)] <- "ends"
  V(sp)$OEM15 <- dme$M15[match(V(sp)$name, dme$locus)]
  V(sp)$OEM15[is.na(V(sp)$OEM15)==T] <- 0
  V(sp)$OEP12 <- dme$P12[match(V(sp)$name, dme$locus)]
  V(sp)$OEP12[is.na(V(sp)$OEP12)==T] <- 0
  V(sp)$Tmp <- dme$Tmp[match(V(sp)$name, dme$locus)]
  V(sp)$Tmp[is.na(V(sp)$Tmp)==T] <- 0
  
  V(sp)$symbol <- str_split_fixed(mods$symbols[match(V(sp)$name, mods$sseqid)], ";",2)[,1]
  V(sp)$symbol[is.na(V(sp)$symbol)==T] <- unis$`Gene.Names.(primary)`[match(V(sp)$name[is.na(V(sp)$symbol)==T], unis$TAIR)]
  V(sp)$symbol[is.na(V(sp)$symbol)==T] <- V(sp)$name[is.na(V(sp)$symbol)==T]
  V(sp)$symbol2 <- NA
  V(sp)$symbol2[V(sp)$type %in% c("root", "ends")] <- V(sp)$symbol[V(sp)$type %in% c("root", "ends")]
  # topdegree <- names(table(degree(sp)))[(length(table(degree(sp)))-5):length(table(degree(sp)))]
  # V(sp)$symbol2[degree(sp) %in% topdegree] <- V(sp)$symbol[degree(sp) %in% topdegree]
  E(sp)$weight <- E(gppi)$lsMR[match(E(sp), E(gppi))]
  E(sp)$color[E(sp)$weight<0] <- "royalblue" 
  E(sp)$color[E(sp)$weight>0] <- "firebrick1"
  E(sp)$color[E(sp)$weight==0] <- "grey"
  E(sp)$weight <- abs(E(sp)$weight)
  E(sp)$weight <- E(sp)$weight / max(E(sp)$weight)
  E(sp)$weight <- E(sp)$weight+0.01
  return(sp)
}