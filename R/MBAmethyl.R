###################################
# Y: observed matrix (p \times n)
# theta: unobserved matrix (p \times n) of true values
# phi: unknown vector of p probe effects
# delta1 = theta[1, ]
# delta2[i, ] = (theta[i + 1, ] - theta[i, ]) / wts[i] for i = 1 : (p - 1)

MBAmethyl <- function(Y, wts = .defaultWeights(nrow(Y)), steps = min(dim(Y)) - 1) {
  
  p <- dim(Y)[1]; n <- dim(Y)[2]
  
  ############################################################
  phi <- rep(1, p); phi.old <- phi + 1e+10; err.phi <- 1e+10
  loop <- 0; theta1 <- NULL
  while ((loop < 20) & (err.phi > 1e-5)) {
    
    theta1 <- t(phi / p) %*% Y; dim(theta1) <- NULL
    
    phi <- Y %*% (theta1 / sum(theta1^2))
    phi <- phi / sqrt(mean(phi^2)); dim(phi) <- NULL
    
    err.phi <- sqrt(sum((phi.old - phi)^2) / p)
    if (err.phi > 1e-5) {
      phi.old <- phi
      rss.0 <- sum((Y - phi %*% t(theta1))^2) / n / p / 2
      loop <- loop + 1
    }
  }
  
  df.0 <- n + p - 1
  
  output <- list()
  output[[1]] <- list(theta = rep(1, p) %x% t(theta1), phi = phi, change.p = numeric(0),
                      rss = rss.0, df.aic = df.0, df.bic = df.0)
  
  #########################################
  for (k in 1 : (steps - 1)) {
    
    phi.old <- phi + 1e+10; err.phi <- 1e+10; loop <- 0
    
    while ((loop < 10) & (err.phi > 1e-5)) {
      
      grp.reg <- .doGFLars(Y, k, phi = phi)
      
      bkp <- grp.reg$bkp
      order.bkp <- order(bkp)
      bkp <- bkp[order.bkp]
      change.p <- bkp + 1
      
      delta2 <- matrix(0, p - 1, n)
      delta2[bkp, ] <- grp.reg$value[[k]][order.bkp, ]
      
      temp <- cumsum(phi[p : 1]^2)[p : 1][-1]
      delta2.up <- delta2 * wts
      delta1 <- t(phi) %*% Y / p - t(temp) %*% delta2.up / p
      
      delta <- rbind(delta1, delta2.up)
      theta <- delta
      for (i in 2 : p) theta[i, ] <- theta[i - 1, ] + delta[i, ]
      
      change.p.set <- list()
      count <- 1; change.p.set[[1]] <- 1
      for (i in 2 : p) {
        if (is.element(i, change.p)) {
          count <- count + 1
          change.p.set[[count]] <- i
        }
        if (!is.element(i, change.p)) change.p.set[[count]] <- c(change.p.set[[count]], i)
      }
      
      for (g in 1 : length(change.p.set)) {
        cp.set.g <- change.p.set[[g]]
        l.g <- length(cp.set.g)
        if (l.g == 1) phi[cp.set.g] <- 1
        if (l.g > 1) {
          theta.g <- theta[cp.set.g[1], ]
          phi[cp.set.g] <- Y[cp.set.g, ] %*% (theta.g / sum(theta.g^2))
          phi[cp.set.g] <- phi[cp.set.g] / sqrt(mean(phi[cp.set.g]^2))
        }
      }
      
      df.aic   <- length(change.p.set) * n + p - length(change.p.set) + 3 * length(change.p)
      df.bic <- length(change.p.set) * n + p - length(change.p.set) + 2 * length(change.p)
      
      rss <- sum((Y - theta * phi)^2) / n / p / 2
      
      err.phi <- sqrt(sum((phi.old - phi)^2) / p)
      if (err.phi > 1e-5) {
        phi.old <- phi
        loop <- loop + 1
        
        output[[k + 1]] <- list(theta = theta, phi = phi, change.p = change.p,
                                rss = rss, df.aic = df.aic, df.bic = df.bic)
      }
    }
  }
  
  rss.path <- df.path.aic <- df.path.bic <- NULL
  for (k in 1 : length(output)) {
    rss.path <- c(rss.path, output[[k]]$rss)
    df.path.aic   <- c(df.path.aic, output[[k]]$df.aic)
    df.path.bic <- c(df.path.bic, output[[k]]$df.bic)
  }
  
  BIC.path.aic   <- log(rss.path) + 2 * df.path.aic / p / n
  BIC.path.bic <- log(rss.path) + log(p) * df.path.bic / p / n
  
  ans.aic   <- output[[which.min(BIC.path.aic)]]; ans.aic$df.bic <- NULL
  ans.bic <- output[[which.min(BIC.path.bic)]]; ans.bic$df.aic <- NULL
  
  return(list(ans.aic = ans.aic, ans.bic = ans.bic))
}
