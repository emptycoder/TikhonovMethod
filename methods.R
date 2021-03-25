TihonovMethod <- function(MTGmatrix, MTVcoefficients, alpha1, alpha2, delta, h, eps) {
  HalfDivisionNeeded = FALSE
  E <- diag(1,ncol(MTGmatrix))
  # Approximation to solve an ill-posed problem
  psi <- matrix(c(0,0),nrow = ncol(MTGmatrix), ncol = 1)
  
  x1 <- solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi)
  x2 <- solve(MTGmatrix + alpha2 * E, MTVcoefficients + alpha2 * psi)
  
  if (det(MTGmatrix) == 0) {
    stop("det of matrix has zero value")
  } else {
    x <- solve(MTGmatrix, MTVcoefficients)
    myu <- norm((matrix %*% x) - vector)
    delta1 <- 2 * norm(x) * h + 2 * delta + myu
  }
  
  phia1 <- norm(matrix %*% x1 - vector) - delta1
  phia2 <- norm(matrix %*% x2 - vector) - delta1
  
  delta11 <- h*norm(x1) + delta1
  delta12 <- h*norm(x2) + delta1

  while (TRUE)
  {
    if (phia1 * phia2 > 0 && phia1 * phia2 < eps / 100) {
      break
    }
    else if (phia1 * phia2 > 0) {
      alpha1 <- alpha1 / 2
      alpha2 <- alpha2 * 2
    }
    else if (phia1 * phia2 < 0) {
      HalfDivisionNeeded = TRUE
      break
    }
    
    x1 <- solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi)
    delta11 <- h * norm(x1) + delta11
    x2 <- solve(MTGmatrix + alpha2 * E, MTVcoefficients + alpha2 * psi)
    delta12 <- h * norm(x2) + delta12
    
    phia1 <- norm(matrix %*% x1 - vector) - delta11
    phia2 <- norm(matrix %*% x2 - vector) - delta12
  }
  
  while (HalfDivisionNeeded)
  {
    if (abs(alpha1 - alpha2) / 2 <= eps) {
      break
    }
    
    m <- (alpha1 + alpha2) / 2
    xm <- solve(MTGmatrix + m * E, MTVcoefficients + m * psi)
    deltam <- h * norm(xm) + delta
    phiam <- norm(matrix %*% xm - vector) - deltam
    if (phia1 * phiam < 0) {
      alpha2 = m
    } else {
      alpha1 = m
    }
  }
  
  return(solve(MTGmatrix + alpha1 * E, MTVcoefficients + alpha1 * psi))
}
