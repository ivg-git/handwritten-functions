#plots beautiful morphology graphs using groups
morph_plots <- function(plant_indices){
  library(rstatix)
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  plants <- levels(as.factor(flo2$Species))
  
  fl1 <- flo2[,-c(1,3)] %>% filter(Species %in% plants[plant_indices])
  fr1 <- fru2[,-c(1,3)] %>% filter(Species %in% plants[plant_indices])
  
  pcr1 <- PCA(X = fl1,quali.sup =  1)
  p0 <- fviz_pca_biplot(pcr1, habillage = "Species", means = F, addEllipses = T, 
                        geom = "point", repel = T, title = "PCA - Biplot. Flowering parameters", palette = "lancet")
  
  pcr2 <- PCA(X = fr1,quali.sup =  1)
  
  p1 <- fviz_pca_biplot(pcr2, habillage = "Species", means = F, addEllipses = T, 
                        geom = "point", repel = T, title = "PCA - Biplot. Fruiting parameters", palette = "lancet")
  p_pca <- ggpubr::ggarrange(p0,p1, ncol = 1, labels = c("A","B"), font.label = list(size = 12))
  
  if(length(plant_indices) > 2){
  
  library(rstatix)
  wl1 <-  fl1 %>% pivot_longer(cols = 2:ncol(fl1), names_to = "Par", values_to = "Val") %>% group_by(Par) %>% wilcox_test(Val~Species)%>%
    mutate(p.adj.significance = -log10(p.adj))
  
  wr1 <-  fr1 %>% pivot_longer(cols = 2:ncol(fr1), names_to = "Par", values_to = "Val") %>% group_by(Par) %>% wilcox_test(Val~Species)%>%
    mutate(p.adj.significance = -log10(p.adj))
  
  kl1 <-  fl1 %>% pivot_longer(cols = 2:18, names_to = "Par", values_to = "Val") %>% group_by(Par) %>% kruskal_test(Val~Species)
  kr1 <-  fr1 %>% pivot_longer(cols = 2:17, names_to = "Par", values_to = "Val") %>% group_by(Par) %>% kruskal_test(Val~Species)
  kl1 <- kl1 %>% mutate(p.sig = ifelse(p<0.05, "sig", "insig"))
  kr1 <- kr1 %>% mutate(p.sig = ifelse(p<0.05, "sig", "insig"))
  p5 <- ggplot(kl1)+
    geom_point(aes(y = Par, x = -log10(p), shape = p.sig), size = 3)+
    scale_shape_manual(values = setNames(c(19,13), c("sig", "insig")))+
    geom_vline(xintercept = -log10(0.05), col = "firebrick1", linetype = "dashed")+theme_light()+
    ylab("Flowering parameters") + xlab("-log10(Kruskal-Wallis P)")+
    theme(axis.text.y = element_text(face = "bold", size = 15))
  p6 <- ggplot(kr1)+
    geom_point(aes(y = Par, x = -log10(p), shape = p.sig), size = 3)+
    scale_shape_manual(values = setNames(c(19,13), c("sig", "insig")))+
    geom_vline(xintercept = -log10(0.05), col = "firebrick1", linetype = "dashed")+theme_light()+
    ylab("Fruiting parameters") + xlab("-log10(Kruskal-Wallis P)")+
    theme(axis.text.y = element_text(face = "bold", size = 15))
  
  p_kw <- ggpubr::ggarrange(p5,p6, ncol = 1, nrow = 2, common.legend = T, legend = "bottom", labels = c("C", "D"))
  
  p3 <- ggplot(wl1)+
    geom_tile(aes(x = group1, y = group2, fill = p.adj.significance))+
    geom_text(aes(x = group1, y = group2, label = p.adj.signif))+
    facet_grid(rows = vars(Par))+scale_fill_gradientn(colours = c("white", "firebrick1"))+theme_light()+
    theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
          axis.text.y = element_text(face = "italic"),
          strip.background = element_rect(fill = "white", colour = "grey"),
          strip.text = element_text(colour = "black", face = "bold"))+xlab("")+ylab("")+
    ggtitle(paste("Flowering"))
  
  p4 <- ggplot(wr1)+
    geom_tile(aes(x = group1, y = group2, fill = p.adj.significance))+
    geom_text(aes(x = group1, y = group2, label = p.adj.signif))+
    facet_grid(rows = vars(Par))+scale_fill_gradientn(colours = c("white", "firebrick1"))+theme_light()+
    theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
          axis.text.y = element_text(face = "italic"),
          strip.background = element_rect(fill = "white", colour = "grey"),
          strip.text = element_text(colour = "black", face = "bold"))+xlab("")+ylab("")+
    ggtitle(paste("Fruiting"))
  
  p_w <- ggpubr::ggarrange(p3,p4, ncol = 2, common.legend = T, legend = "bottom", labels = c("E", "F"))
  
  p_ex <- ggpubr::ggarrange(p_pca,p_kw, p_w, ncol = 3, labels = c("",""), widths = c(1, 0.3, 0.8))
  pdf(paste0("group", plant_indices, "group.pdf"), 18, 12)
  print(p_ex)
  dev.off()
  
  } else {library(rstatix)
    wl1 <-  fl1 %>% pivot_longer(cols = 2:ncol(fl1), names_to = "Par", values_to = "Val") %>% group_by(Par) %>% wilcox_test(Val~Species) %>%
      mutate(p.significance = -log10(p))
    
    wr1 <-  fr1 %>% pivot_longer(cols = 2:ncol(fr1), names_to = "Par", values_to = "Val") %>% group_by(Par) %>% wilcox_test(Val~Species) %>%
      mutate(p.significance = -log10(p))
    
    
    wl1$p.signif = NA_character_
    
    for(i in 1:nrow(wl1)){
      if(wl1$p[i] >=0.05 ){
        wl1$p.signif[i] = "ns"  
      } else if (wl1$p[i] < 0.05 & wl1$p[i] >= 0.01) {wl1$p.signif[i] = "*"
      } else if (wl1$p[i] < 0.01 & wl1$p[i] >= 0.001) {wl1$p.signif[i] = "**"
      } else if (wl1$p[i] < 0.001 & wl1$p[i] >= 0.0001) {wl1$p.signif[i] = "***"
      } else if (wl1$p[i] < 0.0001) {wl1$p.signif[i] = "****"}}
    
    wr1$p.signif <- NA_character_
    
    for(i in 1:nrow(wr1)){
      if(wr1$p[i] >=0.05 ){
        wr1$p.signif[i] = "ns"  
      } else if (wr1$p[i] < 0.05 & wr1$p[i] >= 0.01) {wr1$p.signif[i] = "*"
      } else if (wr1$p[i] < 0.01 & wr1$p[i] >= 0.001) {wr1$p.signif[i] = "**"
      } else if (wr1$p[i] < 0.001 & wr1$p[i] >= 0.0001) {wr1$p.signif[i] = "***"
      } else if (wr1$p[i] < 0.0001) {wr1$p.signif[i] = "****"}}
    
    wl1 <- wl1[order(wl1$Par),]
    wr1 <- wr1[order(wr1$Par),]
    
    wl1$Par <- as.factor(wl1$Par)
    wr1$Par <- as.factor(wr1$Par)
    
    p3 <- ggplot(wl1)+
      geom_tile(aes(x = group1, y = group2, fill = p.significance))+
      geom_text(aes(x = group1, y = group2, label = p.signif))+
      facet_grid(rows = vars(Par))+scale_fill_gradientn(colours = c("white", "firebrick1"))+theme_light()+
      theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
            axis.text.y = element_text(face = "italic"),
            strip.background = element_rect(fill = "white", colour = "grey"),
            strip.text = element_text(colour = "black", face = "bold"))+xlab("")+ylab("")+
      ggtitle(paste("Flowering"))
    
    p4 <- ggplot(wr1)+
      geom_tile(aes(x = group1, y = group2, fill = p.significance))+
      geom_text(aes(x = group1, y = group2, label = p.signif))+
      facet_grid(rows = vars(Par))+scale_fill_gradientn(colours = c("white", "firebrick1"))+theme_light()+
      theme(axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
            axis.text.y = element_text(face = "italic"),
            strip.background = element_rect(fill = "white", colour = "grey"),
            strip.text = element_text(colour = "black", face = "bold"))+xlab("")+ylab("")+
      ggtitle(paste("Fruiting"))
    
    p_w <- ggpubr::ggarrange(p3,p4, ncol = 2, common.legend = T, legend = "bottom", labels = c("C", "D"))
    
    p_ex <- ggpubr::ggarrange(p_pca, p_w, ncol = 3, labels = c("",""), widths = c(1, 0.8))
    pdf(paste0("group", plant_indices, "group.pdf"), 17, 10)
    print(p_ex)
    dev.off()}
  
}
