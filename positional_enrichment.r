positional_enrichment <- function(query_gff, annotation_gff, minoverlap = 30) {
  
  annotation_list <- split(annotation_gff, annotation_gff$type)
  
  intersections <- map_dfr(annotation_list, function(anno) {
    tibble(
      type = unique(anno$type),
      n_intersect = filter_by_overlaps(query_gff, anno, minoverlap = minoverlap) %>% length(),
      n_annotation = length(anno)
    )
  })
  
  results <- intersections %>%
    mutate(
      n_query = length(query_gff),
      n_universe = length(annotation_gff),
      
      m = n_annotation,                    
      n = n_universe - n_annotation,      
      k = n_query,                        
      
      expected = k * (m / (m + n)),
      
      variance = k * (m / (m + n)) * (n / (m + n)) * ((m + n - k) / (m + n - 1)),
      
      sd = sqrt(variance),
      
      z_score = (n_intersect - expected) / sd,
      
      enrichment_log10 = log10((n_intersect + 1e-9) / (expected + 1e-9)),
      
      fold_enrichment = (n_intersect + 1e-9) / (expected + 1e-9),
      
      p_value = phyper(
        q = n_intersect - 1,
        m = m,
        n = n,
        k = k,
        lower.tail = FALSE
      ),
      
      ci_lower = qhyper(0.025, m = m, n = n, k = k),
      ci_upper = qhyper(0.975, m = m, n = n, k = k)
    ) %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      p_significance = -log10(p_adjusted + 1e-9),
      significance = p_adjusted < 0.05
    ) %>%
    arrange(p_adjusted) %>%
    mutate(type = factor(type, levels = type))
  
  p1 <- ggplot(results) +
    geom_point(
      aes(x = enrichment_log10, y = type, 
          size = p_significance, col = z_score, shape = significance),
      alpha = 0.8
    ) +
    geom_vline(xintercept = 0, col = "gray50", linetype = "dotted") +
    geom_vline(xintercept = 1, col = "firebrick1", linetype = "dashed") +
    scale_color_gradientn( colours = c("darkblue", "firebrick1")
                           , name = "Z-score"
    ) +
    scale_shape_manual(values = setNames(nm = c(TRUE, FALSE), c(16, 4))) +
    labs(
      title = "Enrichment of annotations in query regions",
      x = "Enrichment (log10)",
      y = "Annotation type",
      size = "-log10(p_adjusted)"
    ) +
    theme_bw()
  
  print(p1)
  
  return(results)
}
