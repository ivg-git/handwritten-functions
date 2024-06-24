find_LTR <- function(i, mism = 80,LTR_length = 100){
  X <- copd[i]

wha <- vmatchPattern(as.character(subseq(X, start = 1, end = LTR_length)), as.character(X), max.mismatch = mism)
while(length(unlist(wha@ends))>=2){
  LTR_length = LTR_length +1
  wha <- vmatchPattern(as.character(subseq(X, start = 1, end = LTR_length)), as.character(X), max.mismatch = mism)
  cat(paste("LTR_length of ", LTR_length),"\r")}
  
LTR_length = LTR_length - 1 


wha <- vmatchPattern(as.character(subseq(X, start = 1, end = LTR_length)), as.character(X), max.mismatch = mism)

forg <- tibble(seqnames = names(copd)[i], 
               end = unlist(wha@ends),
               start = end - LTR_length, ID = c("LTR_L", "LTR_R"))
if(forg$end[2] > width(copd)[i]){
  dif <- forg$end[2] - width(copd)[i]
  forg$end <- forg$end - dif
LTR_length = LTR_length - dif }
message(paste("LTR_length is ", LTR_length))
res <- as_granges(forg)
return(res)
}
