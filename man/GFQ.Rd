\name{GFQ}
\alias{GFQ}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Greedy Functional Quantization}
\description{Data-driven greedy functional quantization based on the distorsion error or the maximin (a space filling design criterion) }

\usage{GFQ(data,mKL,size,method,deepstart=TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{matrix that we want to quantize.}
\item{mKL}{truncation argument for the dimension reduction.}
\item{size}{size of the quantization grids.}
\item{method}{"L2" or "maximin".}
\item{deepstart}{maximin, if TRUE: the quantization is started by the central curve.}
}

\value{
\item{data}{the input matrix.}
\item{quantizer}{the quantizer grid (curves are chosen among the input data.}
\item{weights}{the associated weight of each curve (calculated using the input matrix).}
}

\seealso{
\code{\link{CVT}} and \code{\link{StochGradient}}
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

data <- t(BM(N=200-1,M=200))
mKL <- 2
size <- 10
method <- "maximin"
quant <- GFQ(data,mKL,size,method,deepstart=TRUE)
}
