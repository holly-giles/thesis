


#*Sum coefficients:* to drop all features from multinomial regression that don't meet specified cut off criteria, and gather coefficients for all remaining other features, for all repeats of the regression. Accepts a matrix `coefMat`, plus numeric values for `coefCut`, the minimum value of that the average coefficient should be, and `freqCut`, the minimum proportion of repeats that a coefficient should be significant.  

sumCoef <- function(coefMat, coefCut = 0, freqCut =1) {
  meanCoef <- rowMeans(abs(coefMat))
  freqCoef <- rowMeans(coefMat != 0)
  subMat <- coefMat[meanCoef > coefCut & freqCoef >= freqCut,,drop=FALSE]
  eachTab <- data.frame(subMat) %>%
     rownames_to_column("feature") %>% gather(key = "rep",value = "coef",-feature) %>%
     mutate(rep = gsub("X","",rep))
  return(eachTab)
} 
