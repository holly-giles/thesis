#select lasso or ridge 
runGlm05 <- function(X, y, method = "lasso", repeats=30, folds = 3) {
  
  #set up objects
  modelList <- list()
  lambdaList <- c()
  varExplain <- c()
  
  #set up a matrix for values of coefficients, with a row for each feature, and a column for each repeat
  coefMat <- matrix(NA, ncol(X), repeats)
  
  #make row names = genetic features
  rownames(coefMat) <- colnames(X)
  #set alpha according to selected method
  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }
  
  #Run cv.glmnet for chosen number of repeats 
  for (i in seq(repeats)) {
    
    #if there are more than two features, fit a glm
    if (ncol(X) > 2) {
      
      res <- cv.glmnet(X,y, type.measure = "mse", family="gaussian", 
                       nfolds = folds, alpha = alpha, standardize = FALSE)
      
      #add lambda min from this repeat to the list of lambdas
      lambdaList <- c(lambdaList, res$lambda.min)
      
      #put the res object (with lambdas) into the list of models
      modelList[[i]] <- res
      
      #extract the coefficients for each feature, for lambda.min
      coefModel <- coef(res, s = "lambda.min")[-1] #remove intercept row
      
      #put these coefficients into  column of coefMatrix corresponding to the repeat
      coefMat[,i] <- coefModel
      
      #calculate variance explained
      y.pred <- predict(res, s = "lambda.min", newx = X)
      varExp <- cor(as.vector(y),as.vector(y.pred))^2
      varExplain[i] <- ifelse(is.na(varExp), 0, varExp) 
      
      
      
    } else {
      #if there are only two features, fit a linear model
      fitlm<-lm(y~., data.frame(X))
      varExp <- summary(fitlm)$r.squared
      varExplain <- c(varExplain, varExp)
      
    }
  }
  #gather all lists
  list(modelList = modelList, lambdaList = lambdaList, varExplain = varExplain, coefMat = coefMat)
}

