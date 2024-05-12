predict_nonb <- function(x){
  library(msa)
  library(gquad)
  pb <- txtProgressBar(min = 0, max = length(x), style = 3)
  results <- list()
  for (i in seq_along(x)){
    gq <- gquad(x[i])
    gqas <- gquad(reverseComplement(x[i]))
    gq$strand = "+"
    gqas$strand = "-"
    gq <- rbind(gq,gqas)
    gq$type <-  "G-quadruplex"
    gq$source <- names(x[i])
    
    sli <- slipped(x[i])
    slias <- slipped(reverseComplement(x[i]))
    sli$strand = "+"
    slias$strand = "-"
    sli <- rbind(sli,slias)
    sli$type <- "slipped"
    sli$source <- names(x[i])
    
    aph <- aphased(x[i])
    aphas <- aphased(reverseComplement(x[i]))
    aph$strand = "+"
    aphas$strand = "-"
    aph <- rbind(aph,aphas)
    aph$type <- "A-phased"
    aph$source <- names(x[i])
    aph$likeliness = ""
    
    str <- str(x[i])
    str$strand = "*"
    str$type <- "Short Tandem Repeat"
    str$source <- names(x[i])
    str$likeliness = ""
    
    tfo <- tfo(x[i])
    tfoas <- tfo(reverseComplement(x[i]))
    tfo$strand = "+"
    tfoas$strand = "-"
    tfo <- rbind(tfo,tfoas)
    tfo$type <- "Triplex-forming oligo"
    tfo$source <- names(x[i])
    tfo$likeliness = ""
    
    zda <- zdna(x[i])
    zdaas <- zdna(reverseComplement(x[i]))
    zda$strand = "+"
    zdaas$strand = "-"
    zda <- rbind(zda,zdaas)
    zda$type <- "Z-DNA"
    zda$source <- names(x[i])
    zda$likeliness = ""
    results[[i]] <- rbind(gq,sli,aph,str,tfo,zda)
    # names(results)[i] <- str_split_fixed(names(x[i]), " ",2)[,1]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  results <- do.call("rbind",results)
  results <- results[-which(results$sequence_position == "-"),]
  results[,c(1,3)] <- lapply(results[,c(1,3)], "as.numeric")
  results$end <- results$sequence_position+results$sequence_length-1
  colnames(results)[c(1,8,7)] <- c("start", "end", "seqnames")
  return(results)
}