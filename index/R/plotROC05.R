#function to plot ROC curves
plotROC05 <- function(glmModel, X, y, lambda = "lambda.1se") {
  lambdaChose <- glmModel[[lambda]]
  glmPred <- prediction(predict(glmModel, type = "response", newx = X, s=lambdaChose), y)
  glmPerform <- performance(glmPred,"tpr","fpr")
  aucVal <- performance(glmPred, measure = "auc")@y.values[[1]]
  xname <- glmPerform@x.name
  yname <- glmPerform@y.name
  plotTab <- tibble(x= glmPerform@x.values[[1]],
                    y = glmPerform@y.values[[1]])
  p <- ggplot(plotTab, aes(x=x, y =y )) + geom_line(color = "red") +
    xlab(xname) + ylab(yname) + theme(panel.grid = element_blank())
  
  if (!is.null(aucVal)) {
    p <- p + annotate("text", x= 0.75, y = 0.25, label = sprintf("AUC: %1.2f", aucVal))
  }
  list(plot=p, auc = aucVal)
}