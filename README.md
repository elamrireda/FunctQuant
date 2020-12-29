---
title: "FunctQuant R Package"
author: "Reda El Amri"
date: "28/12/2020"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r echo=FALSE}
GFQ <- function(data,mKL,size,method,deepstart=TRUE)
{
  support <- FPCA(data,mKL,scale=TRUE)$Coeff
  if(method == 'L2')
  {
    DESIGN <- NULL
    DESIGN <- as.vector(which.min(apply(((support - colMeans(support))^2),1,sum)))
    n <- dim(support)[1]
    for(j in 2:n)
    {
      e2 <- rep(NA,n)
      for(i in 1:n)
      {
        Distance <- as.matrix(dist(rbind(support,support[c(DESIGN,i),])))
        if(j==1) e2[i] <- min((Distance[1:n,(n+1)])[-i])
        else e2[i] <- sum(apply(Distance[1:n,(n+1):(n+(j))],1,min))
      }
      DESIGN <- c(DESIGN,which.min(e2))
      if(j > size) break
    }
  }
  
  if(method=='maximin')
  {
    if(deepstart==FALSE) res <- sample(seq(1,dim(data)[1]),1)
    else res <- as.vector(which.min(apply(((support - colMeans(support))^2),1,sum)))
    
    x <- support[res,]
    dim <- dim(as.matrix(x))[1]
    n <- size+1
    N <- dim(support)[1]
    support <- cbind(1:N,support)
    dist <- matrix(NA*1:(2*N),ncol=2)
    dist[,1] = c(1:N)
    p <- x
    for(j in 1:(size-1))
    {
      for(i in 1:N)
      {
        Distance <- as.matrix(dist(rbind(p,support[i,c(2:(dim+1))])))
        diag(Distance) <- 1.0E30
        min1 <- apply(Distance,2,min)
        dist[i,2] <- min(min1)
      }
      ind <- dist[order(dist[,2],decreasing=TRUE),1][1]
      p <- rbind(p,support[ind,c(2:(dim+1))])
      res <- c(res,ind)
    }
    DESIGN <- res
  }
  
  ccenters <- data[DESIGN,]
  pairs <- data
  distances <- outer(
    1:nrow(pairs),
    1:nrow(ccenters),
    Vectorize( function(m,l) {
      sum( (pairs[m,] - ccenters[l,])^2 )
    } )
  )
  clusters <- apply(distances,1,which.min)
  weight <- as.vector(as.matrix(NA,1,size))
  for(q in 1:size) weight[q] <- length(which(clusters==q))/dim(data)[1]
  return(list(data=data,quantizer=data[DESIGN,],weights=weight,index=DESIGN))
}
FPCA <- function(data,ncomp,scale=FALSE)
{
  dataa <- data
  MEAN <- apply(data,2,mean)
  if(scale== TRUE){ data <- scale((data),center=TRUE,scale=FALSE)}
  else {data <- (data)}
  n <- nrow(data )
  M <- ncol(data )
  
  COVARIANCE <- t(data ) %*% (data )/n
  
  ei=eigen(COVARIANCE, symmetric = TRUE)
  q <- ncomp
  EIGENVECTORS <- matrix(ei$vectors[, 1:q],ncol=q)
  EIGENVALUES <- ei$values[1:q]
  VAR <- sum(ei$values[1:q])/sum(ei$values)
  COMPONENT1 <- data %*% EIGENVECTORS
  donnee <- matrix(0,ncol=M ,nrow=n)
  for(i in 1:n)
  {
    test <- matrix(0,ncol=q,nrow=M)
    for(j in 1:q)
    {
      test[,j] <- COMPONENT1[i,j] * (EIGENVECTORS[,j]  )
    }
    donnee[i,] <- apply(test,1,sum) + MEAN
  }
  
  
  return(list(Data = (dataa), Appro = donnee, eigenvectors = EIGENVECTORS, 
              Coeff = COMPONENT1, eigenvalues = EIGENVALUES,var=VAR,
              vartot=ei$values,center = MEAN,t=ei$values))
}
BM <- function(N=1000,M=1,x0=0,t0=0,T=1,Dt=NULL)
{
  Dt <- (T - t0)/N
  t <- seq(t0, T, by=Dt)
  res <- data.frame(sapply(1:M,function(i) c(0,cumsum(rnorm(N,mean  =0,sd=sqrt(Dt))))))
  names(res) <- paste("X",1:M,sep="")
  X <- ts(res, start = t0, deltat = Dt)
  return(X)
}
```

## Including Plots

You can also embed plots, for example:

```{r echo=TRUE}
data <- t(BM(N = 200 - 1, M = 300))
size <- 3
mKL <- 2
quant <- GFQ(data,mKL,size,method="maximin",deepstart=TRUE)
```
```{r echo=FALSE}
matplot(t(quant$data),type='l',col='grey',xlim=c(0,300),
        xlab="Time (s)",ylab="Functional output",main="Greedy Functional Quantization via maximin method")
for(i in 1:size) lines(quant$quantizer[i,],col=i,lwd=3)
for(i in 1:size)  points(200,quant$quantizer[i,200],pch=19,col=i)
for(i in 1:size)  text(250,quant$quantizer[i,200],labels=paste("weight[",i,"]=",round(quant$weights[i],3),sep=""),col=i)
```

```{r echo=TRUE}
quant <- GFQ(data,mKL,size,method="L2",deepstart=TRUE)
```

```{r echo=FALSE}
matplot(t(quant$data),type='l',col='grey',xlim=c(0,300),
        xlab="Time (s)",ylab="Functional output",main="Greedy Functional Quantization via L2 method")
for(i in 1:size) lines(quant$quantizer[i,],col=i,lwd=3)
for(i in 1:size)  points(200,quant$quantizer[i,200],pch=19,col=i)
for(i in 1:size)  text(250,quant$quantizer[i,200],labels=paste("weight[",i,"]=",round(quant$weights[i],3),sep=""),col=i)
```

