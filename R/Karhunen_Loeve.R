Karhunen_Loeve <-
function(data,mKL,scale=TRUE)
{
  dataa <- data
  MEAN <- apply(data,2,mean)
  if(scale== TRUE){ data <- scale((data),center=TRUE,scale=FALSE)}
  else {data <- (data)}
  n <- nrow(data )
  M <- ncol(data )
  COVARIANCE <- t(data ) %*% (data )/n
  ei=eigen(COVARIANCE, symmetric = TRUE)
  q <- mKL
  EIGENVECTORS <- matrix(ei$vectors[, 1:q],ncol=q)
  EIGENVALUES <- ei$values[1:q]
  VAR <- sum(ei$values[1:q])/sum(ei$values)
  COMPONENT1 <- data %*% EIGENVECTORS
  return(list(Coeff=COMPONENT1,ExplainedVar=VAR,center=MEAN,eigenfunction=EIGENVECTORS))
}
