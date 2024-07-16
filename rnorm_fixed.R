rnorm_fixed = function(n, mu=0, sigma=1) { #mu = mean, sigma = sd
  x = rnorm(n)  
  x = sigma * x / sd(x) 
  x = x - mean(x) + mu 
  return(x)
}
