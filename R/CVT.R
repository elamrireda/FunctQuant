CVT <- function(data,size,iter=22)
{
  Nsim <- dim(data)[1]
  discretization <- dim(data)[2]
  id <- sample(seq(1,Nsim,by=1),size)
  center <- data[id,]
  for(j in 1:iter)
  {
    ccenters <- center
    pairs <- (data)
    distances <- outer(
      1:nrow(pairs),
      1:nrow(ccenters),
      Vectorize( function(m,l) {
        sqrt(sum( (pairs[m,] - ccenters[l,])^2 )/discretization)
      } )
    )
    clusters <- apply(distances,1,which.min)
    weight <- rep(0,size)
    for(q in 1:size) weight[q] <- length(which(clusters==q))/Nsim
    for(i in 1:size)
    {
      if( length(which(clusters==i)) > 1) center[i,] <- apply(data[which(clusters==i),],2,mean)
      if( length(which(clusters==i)) == 1) center[i,] <- data[which(clusters==i),]
      if( length(which(clusters==i)) == 0) center[i,] <- center[i,]
    }
    data <- data
  }
  return(list(data=data,quantizer=center,weights=weight))
}
