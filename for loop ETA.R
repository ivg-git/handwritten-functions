library(lubridate)
#result <- list()
n_iter <- 256
init = numeric(n_iter)
end = numeric(n_iter)
for (i in 1:n_iter){
  init[i] <- Sys.time()
  
  #result[[i]] <-  
  
  end[i] <- Sys.time()
  time <- round(seconds_to_period(sum(end - init)), 0)
  est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
  remainining <- round(seconds_to_period(est), 0)
  
  cat(paste("iter", i,"/",round(i/n_iter*100, digits = 2), "%",  " // Execution time:", time,
            " // Estimated time remaining:", remainining, "     "), "",sep="\r")
  flush.console() }