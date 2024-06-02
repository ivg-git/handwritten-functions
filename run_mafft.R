run_mafft <- function(x){
  writeXStringSet(x, "temp.fa")
  system("E:/coursework/SIRC/mafft-win/mafft.bat --thread -1  --out out.fa --localpair temp.fa", invisible = T, ignore.stderr = T)
  readDNAStringSet("out.fa")}