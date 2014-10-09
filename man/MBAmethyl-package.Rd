\name{MBAmethyl-package}
\alias{MBAmethyl-package}
\docType{package}
\title{
Model-based analysis of DNA methylation data
}
\description{
This package provides functions for reconstructing DNA methylation 
values from raw measurements. It utilize both the information from biological replicates and neighboring probes by explicitly modeling the probe-specific effect and encouraging the neighboring similarity by a group fused lasso penalty. 
}
\details{
\tabular{ll}{
Package: \tab MBAmethyl\cr
Type: \tab Package\cr
Version: \tab 0.99.0 \cr
Date: \tab 2014-08-24\cr
License: \tab  Artistic-2.0 \cr
}
}
\author{
Tao Wang, Mengjie Chen

Maintainer: Tao Wang <tao.wang.tw376@yale.edu>

}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\examples{
 p <- 80
  n <- 40
  K <- 2
  k <- K - 1
  cp <- numeric()
  L <- c(0, floor(p / K) * (1 : k), p)
  cp <- floor(p / K) * (1 : k) + 1

  ## phi0: probe effects; theta0: true methylation values; part: partition of probe indices
  phi0 <- runif(p, 0.5, 2.0)
  theta0 <- matrix(0, p, n)
  part <- list()

  for (s in 1 : K) {
    part[[s]] <- (L[s] + 1) : L[s + 1]
    phi0[part[[s]]] <- phi0[part[[s]]] / sqrt(mean(phi0[part[[s]]]^2))
  }

  theta0[part[[1]], ] <- rep(1, length(part[[1]])) %x% t(runif(n, 0.1, 0.6))
  theta0[part[[2]], ] <- rep(1, length(part[[2]])) %x% t(runif(n, 0.4, 0.9))
  
  error <- matrix(runif(p * n, 0, 0.1), p, n)
  Y <- theta0 * phi0  + error
  fit <- MBAmethyl(Y, steps = 10)
}
