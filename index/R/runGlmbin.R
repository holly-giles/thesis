runGlm.bin <- function(X, y, method = "lasso", repeats=20, folds = 3, testRatio = NULL, lambda = "lambda.1se") {
  modelList <- list()
  lambdaList <- c()
  aucCV <- c()
  aucTest <- c()
  rocTest <- list()
  coefMat  <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)
  
  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }
  
  for (i in seq(repeats)) {
    if (!is.null(testRatio)) {
      testIdx <- testPartition(y, testRatio)
      X.test <- X[testIdx,]
      X.train <- X[-testIdx,]
      y.test <- y[testIdx]
      y.train <- y[-testIdx]
    } else {
      X.train <- X
      y.train <- y
    }
    
    #need to add a set seed for this function, but ensure that it changes with each loop otherwise all 10 repeats are the same 
    vecFold <- mltools::folds(y.train, nfolds = folds, stratified = TRUE, seed = i)
    
    #train model
    res <- cv.glmnet(X.train,y.train, type.measure = "auc",
                     foldid = vecFold, alpha = alpha, standardize = FALSE,
                     intercept = TRUE, family = "binomial")
    lambdaList <- c(lambdaList, res[[lambda]])
    #If you arent splitting the data into a training and testing dataset. you can use the AUC for the CV
    aucCV <- c(aucCV, res$cvm[res$lambda == res[[lambda]]])
    modelList[[i]] <- res
    coefMat[,i] <- coef(res, s = lambda)[-1]
    
    #test model if testRatio is speficied- can use this AUC to identify the best model 
    if(!is.null(testRatio)) {
      rocRes <- plotROC(res, X.test, y.test, lambda)
      aucTest <- c(aucTest, rocRes$auc)
      rocTest[[i]] <- rocRes$plot
    }
  }
  list(modelList = modelList, lambdaList = lambdaList, aucCV = aucCV, coefMat = coefMat,
       aucTest = aucTest, rocTest = rocTest)
}