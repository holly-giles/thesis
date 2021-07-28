#indicate how to split test set and training set
testPartition05 <- function(y, ratio) {
  #balanced sampling of test set
  ord <- seq_along(y)
  testIdx <- lapply(unique(y),function(n) {
    subY <- ord[y == n]
    sample(subY, size = as.integer(length(subY)  * ratio)) 
  }) %>% do.call(c,.) %>% sort()
  return(testIdx)
}