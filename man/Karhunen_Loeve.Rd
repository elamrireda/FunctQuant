\name{Karhunen_Loeve}
\alias{Karhunen_Loeve}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Karhunen Loeve expansion}
\description{Dimension reduction : Karhunen Loeve expansion}


\usage{Karhunen_Loeve(data,mKL,scale=TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{matrix that we want to quantize.}
\item{mKL}{truncation argument for the dimension reduction.}
\item{scale}{scale or not the input data.}
}

\value{
\item{Coeff}{the coefficient in the reduced space.}
\item{ExplainedVar}{the explained variance.}
\item{center}{the mean curve.}
\item{eigenfunction}{the eigenfunctions.}
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
mKL <- 2
DimReduction <- Karhunen_Loeve(data,mKL,scale=TRUE)
}
