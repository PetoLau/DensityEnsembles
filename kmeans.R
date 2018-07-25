kmeanspp <- function(x, k, iter.max = 10, nstart = 1, ...) {
  
  res <- kmeans(x, k, iter.max = iter.max, nstart = nstart)
  
  res
}
