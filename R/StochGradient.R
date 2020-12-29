StochGradient <- function(data,mKL,size)
{
  choice.grid_Reda <- function(X, N, ng=1, p = 2) {
    if (!is.numeric(X))
      stop("X must be numeric")
    if (!is.vector(X) & !is.matrix(X))
      stop("X must be a matrix or a vector")
    if ((!(floor(N) == N)) | (N <= 0))
      stop("N must be entire and positive")
    if ((!(floor(ng) == ng)) | (ng <= 0))
      stop("B must be entire")
    if (p < 1)
      stop("p must be at least 1")
    if (is.vector(X)) {
      d <- 1
      n <- length(X)
      primeX <- matrix(sample(X, n * ng, replace = TRUE),
                       nrow = ng)
      hatX <- replicate(ng, sample(unique(X), N, replace = FALSE))
      hatX0 <- hatX  #save the initial grids
      gammaX <- array(0, dim = c(ng, n + 1))
      # initialisation of the step parameter gamma
      a_gamma <- 4 * N
      b_gamma <- pi^2 * N^(-2)
      BtestX = array(Inf, dim = c(ng, 1))
      # choice of gamma0X
      projXbootinit <- array(0, dim = c(n, ng))
      # index of the grid on which X is projected
      iminx <- array(0, dim = c(n, ng))
      for (i in 1:n) {
        RepX <- matrix(rep(X[i], N * ng), ncol = ng, byrow = TRUE)
        Ax <- sqrt((RepX - hatX0)^2)
        iminx[i, ] <- apply(Ax, 2, which.min)
        mx <- matrix(c(iminx[i, ], c(1:ng)), nrow = ng)
        projXbootinit[i, ] <- hatX0[mx]
      }
      RepX <- matrix(rep(X, ng), ncol = ng)
      distortion <- apply((RepX - projXbootinit)^2, 2, sum)/n
      temp_gammaX <- which(distortion > 1)
      if (length(temp_gammaX) > 0) {
        distortion[temp_gammaX] <- array(1, dim = c(length(temp_gammaX),
                                                    1))
      }
      gamma0X <- distortion
      if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
      gammaX = array(0, dim = c(ng, n + 1))
      for (i in 1:ng) {
        # calculation of the step parameter
        gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] *
                                               b_gamma * c(1:(n + 1) - 1))
      }
      gammaX[, 1] <- gamma0X
      iminX <- array(0, dim = c(n, ng))  #index that will change
      tildeX <- array(0, dim = c(N, ng))
      # updating of the grids, providing optimal grids
      for (i in 1:n) {
        for (j in 1:ng) {
          tildeX[, j] <- matrix(rep(primeX[j, i], N), nrow = N,
                                byrow = TRUE)
        }
        Ax <- sqrt((tildeX - hatX)^2)
        # calculation of each distance to determine the point of
        # the grid the closer of the stimuli
        iminX[i, ] <- apply(Ax, 2, which.min)
        mX <- matrix(c(iminX[i, ], c(1:ng)), nrow = ng)
        if(sqrt(sum((hatX[mX] - primeX[, i])^2))==0){
          hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] - primeX[, i])
        }else{
          hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] -
                                                      primeX[, i]) * (sqrt(sum((hatX[mX] - primeX[,
                                                                                                  i])^2)))^(p - 1)/sqrt(sum((hatX[mX] - primeX[,                                                                                                                          i])^2))
        }
      }
    } else {
      n <- ncol(X)
      d <- nrow(X)
      primeX <- array(X[, sample(c(1:n), n * ng,
                                 replace = T)], dim = c(d, n, ng))
      hatX <- replicate(ng, unique(X)[, sample(c(1:n),
                                               N, replace = F)])  # initial grids chosen randomly in the sample
      hatX <- array(hatX, dim = c(d, N, ng))
      hatX0 <- hatX  #save the initial grids
      # initialisation of the step parameter gamma
      a_gamma <- 4 * N^(1/d)
      b_gamma <- pi^2 * N^(-2/d)
      BtestX <- array(Inf, dim = c(1, ng))
      # choice of gamma0X
      for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          Bx <- array(0, dim = c(ng, 1))
          Bx <- sqrt(apply((hatX[, i, , drop = FALSE] -
                              hatX[, j, , drop = FALSE])^2, c(2, 3), sum))/2
          temp_gammaX <- which(Bx < BtestX)
          BtestX[temp_gammaX] <- Bx[temp_gammaX]
        }
      }
      temp_gammaX = which(BtestX > 1)
      if (length(temp_gammaX) > 0) {
        BtestX[temp_gammaX] <- array(1, dim = c(length(temp_gammaX),
                                                1))
      }
      gamma0X <- BtestX
      if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
      
      gammaX <- array(0, dim = c(ng, n + 1))
      for (i in 1:ng) {
        # calculation of the step parameter
        gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] *
                                               b_gamma * c(1:(n + 1) - 1))
      }
      gammaX[, 1] <- gamma0X
      iminX <- array(0, dim = c(n, ng))  #index that will change
      tildeX <- array(0, dim = c(d, N, ng))
      # updating of the grids, providing optimal grids
      for (i in 1:n) {
        for (j in 1:ng) {
          tildeX[, , j] <- matrix(rep(primeX[, i, j], N),
                                  nrow = N, byrow = FALSE)
        }
        Ax <- sqrt(apply((tildeX - hatX)^2, c(2, 3), sum))
        # calculation of each distance to determine the point of the grid
        #the closer of the stimuli
        iminX[i, ] <- apply(Ax, 2, which.min)
        for (k in 1:d) {
          m <- matrix(c(rep(k, ng), iminX[i, ], c(1:ng)), ncol = 3)
          if(p==2){
            hatX[m] = hatX[m] - gammaX[, i + 1] * (hatX[m] - primeX[k, i, ])
          }else{
            hatX[m] = hatX[m] - gammaX[,i + 1] * (hatX[m] - primeX[k, i, ]) *
              (sqrt(sum((hatX[m] - primeX[k, i, ])^2)))^(p - 1)/sqrt(sum((hatX[m] - primeX[k,
                                                                                           i, ])^2))
          }
        }
      }
    }
    output <- list(init_grid=hatX0,opti_grid=hatX)
    output
  }
  
  dimRed <- Karhunen_Loeve(data,mKL,scale=TRUE)
  support <- dimRed$Coeff
  MEAN <- dimRed$center
  
  l <- choice.grid_Reda(t(support),size,ng=1,p=2)
  
  newCoeff <- t(l$opti_grid[,,1])
  
  quantizer <- matrix(0,ncol=ncol(data),nrow=size)
  for(i in 1:size)
  {
    test <- matrix(0,ncol=mKL,nrow=ncol(data))
    for(j in 1:mKL) test[,j] <- newCoeff[i,j] * (dimRed$eigenfunction[,j]  )
    quantizer[i,] <- apply(test,1,sum) + MEAN
  }
  
  ccenters <- quantizer
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
  
  return(list(data=data,quantizer=quantizer,weights=weight))
}
