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