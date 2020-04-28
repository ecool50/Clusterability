run.SVM <- function(xdata, xlabel, ydata,kern="linear")
{
  f.g <- intersect(colnames(xdata),colnames(ydata))
  xdata <- xdata[,f.g,drop=F]
  ydata <- ydata[,f.g,drop=F]
  model <- e1071::svm(xdata, xlabel, kernel = kern)
  ylabel <- predict(model, newdata=ydata)
  names(ylabel) <- rownames(ydata)
  return(list("ylabel"=ylabel,"svm"=model))
}