\name{MBAmethyl}
\alias{MBAmethyl}
\title{
Model-based analysis of DNA methylation data
}
\description{
This function reconstructs DNA methylation values from raw measurements. It 
iteratively implements the group fused lars to smooth related-by-location 
methylation values and the constrained least squares to remove probe affinity 
effect across multiple sequences. It also contains a criterion-based method 
(AIC or BIC) for selecting the tuning paramter.  
}
\usage{
MBAmethyl(Y, wts = .defaultWeights(nrow(Y)), steps = min(dim(Y)) - 1)
}
\arguments{
  \item{Y}{An observed matrix (p x n) of methylation values (beta values); p is the number of probes and n is the number of samples;}
  \item{wts}{A pre-specified vector of weights. By default, we use the probe index-dependent weight scheme, $wts_i = sqrt(p / i / (p - i))$ for $i = 1, ... , p$; }
  \item{steps}{Limit the number of steps taken. One can use this option to 
perform early stopping.}
}

\value{
\item{ans.aic}{A list corresponds to the AIC, 
containing estimated beta values, estimated probed effects, estimated 
change-point locations, residual sum of squares, and degree of freedom. }

\item{ans.bic}{A list corresponds to the BIC, 
containing estimated beta values, estimated probed effects, estimated 
change-point locations, residual sum of squares, and degree of freedom. }
}
\references{
paper under review
}
\author{
Tao Wang, Mengjie Chen
}
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
\keyword{ methylation }

