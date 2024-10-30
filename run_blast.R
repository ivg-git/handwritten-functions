makeblastdb <- function(input_path, db_name = "bdb"){
  system(paste("makeblastdb -in ", input_path, "-parse_seqids -blastdb_version 5 -title", db_name,
               "-dbtype prot -out bdb"))
  return(print("All clear!"))
}
run_blastp <- function(query, db_name = "bdb", num_threads = parallel::detectCores()){
  writeXStringSet(query, "query.fa")
  system(paste("blastp -query query.fa -db", db_name, "-out tempoutput.tab -outfmt 6 -num_threads", num_threads))
  temp <- read.delim("tempoutput.tab", header = F)
  colnames(temp) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart" ,"send" ,"evalue" ,"bitscore")
  return(temp)
}
