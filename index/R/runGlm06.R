#Function for multi-variant regression
runGlm06 <- function(X, y, method = "lasso", repeats = 20, folds = 3) {
  #set up objects to store results
  modelList <- list()
  lambdaList <- c()
  varExplain <- c()
  coefMat <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)
  
  #set alpha
  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }
  
  #for the set number of repeats, run the regression
  for (i in seq(repeats)) {
    if (ncol(X) > 2) {
      #run cross validated generalised linear model with given parameters
      res <- cv.glmnet(X,y, type.measure = "mse", family="gaussian", 
                       nfolds = folds, alpha = alpha, standardize = FALSE)
      
      #store lamdas and min lambda value
      lambdaList <- c(lambdaList, res$lambda.min)
      
      #store result of cv.glmnet
      modelList[[i]] <- res
      
      #get coefficents with min lambda value 
      coefModel <- coef(res, s = "lambda.min")[-1] #remove intercept row
      
      #store coefficients for this repeat
      coefMat[,i] <- coefModel
      
      #calculate variance explained
      if(sum(coefModel !=0)){
        y.pred <- predict(res, s = "lambda.min", newx = X)
        #if there are no predictors, all y.pred will be the same so its not possible to calculate variance explained because the SD is 0
        varExp <- cor(as.vector(y),as.vector(y.pred))^2
      }else{ varExp <- NA}
      varExplain[i] <- ifelse(is.na(varExp), 0, varExp) 
      
      
      
    } else {
      fitlm<-lm(y~., data.frame(X))
      varExp <- summary(fitlm)$r.squared
      varExplain <- c(varExplain, varExp)
      
    }
    
  }
  #store all results 
  list(modelList = modelList, lambdaList = lambdaList, varExplain = varExplain, coefMat = coefMat)
}
