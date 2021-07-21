

#*Function for multinomial regression*: To perform multinomial regression to identify genetic features that are predictors of cluster assignment. Provide a feature matrix `X` and a response matrix `y`, and specify `method` (regression method), `repeats` (number of repeats of the regression) and `folds` (number of folds to split the data into for cross validation). 

runGlm_multiCluster <- 
  function(X, y, method = "ridge", repeats=20, folds = 3) {
    modelList <- list()
    lambdaList <- c()
    
    coefMat <- 
      lapply(unique(y), function(n) {
        mat <- matrix(NA, ncol(X), repeats)
        rownames(mat) <- colnames(X)
        mat
      })
    
    names(coefMat) <- unique(y)
    
    alpha = switch(method, lasso = 1, ridge = 0, stop("Please provide a valid method: lasso or ridge"))
    
    
    for (i in seq(repeats)) {
      
      #balanced sampling
      vecFold <- mltools::folds(y, nfolds = folds, stratified = TRUE, seed = i*1996)
      res <- cv.glmnet(X, y, type.measure = "class",
                       foldid = vecFold, alpha = alpha, standardize = FALSE,
                       intercept = TRUE, family = "multinomial")
      lambdaList <- c(lambdaList, res$lambda.min)
      modelList[[i]] <- res
      
      coefModel <- coef(res, s = "lambda.min")
      for (n in names(coefModel)) {
        coefMat[[n]][,i] <- coefModel[[n]][-1]
      }
    }
    list(modelList = modelList, lambdaList = lambdaList, coefMat = coefMat)
  }



