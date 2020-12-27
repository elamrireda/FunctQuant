\name{StochGradient}
\alias{StochGradient}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functional quantization with stochastic gradient method}
\description{Data-driven functional quantization with stochastic gradient method}

\usage{StochGradient(data,mKL,size)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{matrix that we want to quantize.}
\item{mKL}{truncation argument for the dimension reduction.}
\item{size}{size of the quantization grids.}
}

\value{
\item{data}{the input matrix.}
\item{quantizer}{the quantizer grid.}
\item{weights}{the associated weight of each curve (calculated using the input matrix).}
}

\seealso{
\code{\link{CVT}} and \code{\link{GFQ}}
}

\examples{
##### function to generate realizations of BM
BM <- function(N=1000,M=1,x0=0,t0=0,T=1,Dt=NULL)
{
  Dt <- (T - t0)/N
  t <- seq(t0, T, by=Dt)
  res <- data.frame(sapply(1:M,function(i) c(0,cumsum(rnorm(N,mean  =0,sd=sqrt(Dt))))))
  names(res) <- paste("X",1:M,sep="")
  X <- ts(res, start = t0, deltat = Dt)
  return(X)
}
data <- t(BM(N = 200 - 1, M = 200))
size <- 10
mKL <- 2
quant <- StochGradient(data,mKL,size)
}
