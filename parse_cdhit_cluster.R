parse_cdhit_clstr <- function(clstr_file) {
  # Читаем файл построчно
  lines <- readLines(clstr_file)
  
  # Инициализируем пустые векторы
  clusters <- c()
  seq_ids <- c()
  
  current_cluster <- NA
  
  # Проходим по каждой строке
  for (line in lines) {
    # Если строка начинается с ">Cluster" — это начало нового кластера
    if (grepl("^>Cluster", line)) {
      # Извлекаем номер кластера (число после "Cluster")
      current_cluster <- as.numeric(gsub(">Cluster (\\d+)", "\\1", line))
    } else {
      # Это строка с последовательностью
      # Формат: "0\t92nt, >NC_082286.1_134... *"
      # или "1\t60nt, >NC_083147.1_1561... at 1:60:33:92/+/98.33%"
      
      # Извлекаем ID последовательности — всё, что между ">" и первым пробелом или "..."
      # Паттерн: ищем ">" за которым идет ID до пробела или "..." или конца строки
      seq_id <- gsub(".*>([^ ]+).*", "\\1", line)
      
      # Если ID содержит "..." — обрезаем его
      seq_id <- gsub("\\.\\.\\.$", "", seq_id)
      
      # Добавляем в векторы
      clusters <- c(clusters, current_cluster)
      seq_ids <- c(seq_ids, seq_id)
    }
  }
  
  # Создаем data.frame
  result <- data.frame(
    cluster = clusters,
    seq_id = seq_ids,
    stringsAsFactors = FALSE
  )
  
  return(result)
}
