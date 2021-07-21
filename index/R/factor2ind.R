
#*Factor2Ind*: To generate indicator matrix from a factor. Given a factor  `x`, create an indicator matrix of dimension `length(x)` multiplied by `nlevels(x)-1`, dropping the column corresponding to the baseline level (by default the first level is used as baseline).

factor2ind <- function(x, baseline)
{
  
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  X[is.na(x),] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  return(X[,-1,drop=FALSE])
}

