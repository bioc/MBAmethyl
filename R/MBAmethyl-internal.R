.doGFLars <- function (Y, K, phi = rep(1, nrow(Y)), weights = .defaultWeights(nrow(Y)), epsilon = 1e-09, verbose = FALSE)
{
  if (is.null(dim(Y)) || is.data.frame(Y)) {
    if (verbose) {
      print("Coercing 'Y' to a matrix")
    }
    Y <- as.matrix(Y)
  }
  else if (!is.matrix(Y)) {
    stop("Argument 'Y' should be a matrix, vector or data.frame")
  }
  if (any(is.na(Y))) {
    stop("Missing values are not handled by the current implementation of group-fused LARS")
  }
  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]
  if (K >= n) {
    stop("Too many breakpoints are required")
  }
  if (is.null(weights)) {
    weights <- .defaultWeights(n)
  }
  res.meansignal <- colMeans(Y)
  res.lambda <- numeric(K)
  res.bkp <- numeric(K)
  res.value <- list()
  res.c <- matrix(NA, n - 1, K)
  Y <- Y - (phi / n) %*% (t(phi) %*% Y)
  AS <- numeric(0)
  c <- .leftMultiplyByXt(Y = Y, phi = phi, w = weights, verbose = verbose)
  for (ii in 1:K) {
    cNorm <- rowSums(c^2)
    res.c[, ii] <- cNorm
    bigcHat <- max(cNorm)
    if (verbose) {
      print(paste("optimize LARS : ", ii))
    }
    if (ii == 1) {
      AS <- which.max(cNorm)
      res.bkp[ii] <- AS
    }
    I <- order(AS)
    AS0 <- AS
    AS <- AS[I]
    w <- .leftMultiplyByInvXAtXA(n, AS, matrix(c[AS, ], ncol = p), phi,
                                weights, verbose = verbose)
    a <- .multiplyXtXBySparse(n = n, ind = AS, val = w, phi = phi, w = weights,
                             verbose = verbose)
    a1 <- bigcHat - rowSums(a^2)
    u <- a * c
    a2 <- bigcHat - rowSums(u)
    a3 <- bigcHat - cNorm
    gammaTemp = matrix(NA, n - 1, 2)
    subset <- which(a1 > epsilon)
    delta <- a2[subset]^2 - a1[subset] * a3[subset]
    delta.neg <- subset[which(delta < 0)]
    delta.pos <- subset[which(delta >= 0)]
    gammaTemp[delta.neg, 1] <- NA
    gammaTemp[delta.pos, 1] <- (a2[delta.pos] + sqrt(delta[which(delta >=
                                                                   0)]))/a1[delta.pos]
    gammaTemp[delta.neg, 2] <- NA
    gammaTemp[delta.pos, 2] <- (a2[delta.pos] - sqrt(delta[which(delta >=
                                                                   0)]))/a1[delta.pos]
    subset <- which((a1 <= epsilon) & (a2 > epsilon))
    gammaTemp[subset, ] = a3[subset]/(2 * a2[subset])
    maxg <- max(gammaTemp, na.rm = TRUE) + 1
    subset <- which((a1 <= epsilon) & (a2 <= epsilon))
    gammaTemp[subset, ] <- maxg
    gammaTemp[AS, ] <- maxg
    
    gammaTemp[gammaTemp <= 0] <- maxg
    
    gamma <- min(gammaTemp, na.rm = TRUE)
    idx <- which.min(gammaTemp)
    nexttoadd <- 1 + (idx - 1)%%(n - 1)
    res.lambda[ii] <- sqrt(bigcHat)
    res.value[[ii]] <- matrix(numeric(ii * p), ncol = p)
    res.value[[ii]][I, ] <- gamma * w
    if (ii > 1) {
      res.value[[ii]][1:(ii - 1), ] <- res.value[[ii]][1:(ii -
                                                            1), ] + res.value[[ii - 1]]
    }
    if (ii < K) {
      AS <- c(AS0, nexttoadd)
      res.bkp[ii + 1] <- nexttoadd
      c <- c - gamma * a
    }
  }
  return(list(bkp = res.bkp, lambda = res.lambda, mean = res.meansignal,
              value = res.value, c = res.c))
}

###
.defaultWeights <- function (n)
{
  a <- seq(length = n - 1)/n
  b <- a * (1 - a)
  1/sqrt(b * n)
}

###
.leftMultiplyByXt <- function (Y, phi = rep(1, nrow(Y)), w = .defaultWeights(nrow(Y)), verbose = FALSE)
{
  n <- as.numeric(dim(Y)[1])
  p <- ncol(Y)
  
  Y <- Y * phi
  
  u <- apply(Y, 2, cumsum)
  if (length(w) != (n - 1)) {
    stop("w needs to be of length nrow(Y)-1")
  }
  C <- apply(u, 2, function(x) {
    w * (cumsum(phi^2)[1:(n - 1)] * x[n]/n - x[1:(n - 1)])
  })
  dim(C) <- c(n - 1, p)
  return(C)
}

###
.leftMultiplyByInvXAtXA <- function (n, ind, val, phi = rep(1, n), w = .defaultWeights(n), verbose = FALSE)
{
  a <- dim(val)[1]
  p <- dim(val)[2]
  o <- order(ind)
  ind <- ind[o]
  val <- val[o, , drop = FALSE]
  r <- matrix(numeric(a * p), nrow = a, ncol = p)
  if (length(w) != (n - 1)) {
    stop("'w' needs to be of length n-1")
  }
  if (a != 0) {
    v <- diff(c(0, cumsum(phi^2)[ind], n))
    d <- w[ind]
    R <- matrix(numeric((a + 2) * p), ncol = p)
    val <- apply(val, 2, function(x) {
      x/d
    })
    R[1, ] <- numeric(p)
    R[2:(a + 1), ] <- val
    R[(a + 2), ] <- numeric(p)
    gamma <- apply(R, 2, diff)
    delta <- gamma/v
    r <- -diff(delta)
    r <- r/d
  }
  return(r)
}

###
.multiplyXtXBySparse <- function (n, ind, val, phi = rep(1, n), w = .defaultWeights(n), verbose = FALSE)
{
  if (length(w) != (n - 1)) {
    stop("Argument 'w' has to be of length nrow(Y)-1")
  }
  a <- nrow(val)
  p <- ncol(val)
  indrev <- (n - 1):1
  if (a != 0) {
    o <- order(ind)
    ind <- ind[o]
    val <- val[o, , drop = FALSE]
    r <- val * w[ind]
    s <- cumsum(phi^2)[ind] * r
    s <- colSums(s)/n
    T <- matrix(numeric((n - 1) * p), nrow = (n - 1), ncol = p)
    T[indrev[ind], ] <- r[, , drop = FALSE]
    T <- apply(T, 2, cumsum)
    T <- T[indrev, , drop = FALSE]
    u <- sweep(T, 2, s)
    
    u <- u * (phi[-n])^2
    
    U <- apply(u, 2, cumsum)
    C <- U * w
  }
  return(C)
}
