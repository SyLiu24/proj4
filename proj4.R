#XXXXXXXXXXXXXXXXXXXXXXXX proj2.r XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Group member
# Biwei Zhu, s2325784
# Guanhao Su, s2301705
# Shuying Liu, s2436365

# Github repo: 
# https://github.com/SyLiu24/proj4.git

#XXXXXXXXXXXXXXXXXXXXXXXXXXX Contribution XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX





newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
    
  n <- length(theta)
  
  # If hess is not given, then estimate it by finite differencing of grad
  if (is.null(hess)) {
    hess <- function(theta,...) {
      # n <- length(theta)
      g <- grad(theta,...)
      hfd <- matrix(0,n,n)
      for (i in 1:n) {
        thi <- theta
        thi[i] <- thi[i] + eps
        gi <- grad(thi,...)
        hfd[i,] <- (gi - g)/eps
      }
      # Make sure hessian is symmetric
      (t(hfd) + hfd)/2
    }
  }


  th0 <- theta
  f0 <- func(th0,...)
  g0 <- grad(th0,...)
  h0 <- hess(th0,...)

  # Check if the initial theta is suitable for optimization
  # Check if objective function value is -inf
  if (is.infinite(f0) & f0 < 0)
    stop("Function is unbounded below.")
  # Check if objective or derivatives are not finite
  else if (is.infinite(f0) | any(is.infinite(g0)) | any(is.infinite(h0)))
    stop("Objective or derivatives are not finite at the initial theta.")

  it <- 0
  while (any(abs(g0) >= tol*abs(f0+fscale))) {
  
  
    if (it == maxit)
      stop("Maxit is reached without convergence.")
  
  
    # If h0 is not positive definite then perturb it to be
    # by adding a multiple of identity matrix
    if (inherits(try(chol(h0),silent=TRUE),"try-error")) {
      m <- (1e-6)*norm(h0) # Multiplier of identity matrix
    
      while (inherits(try(chol(h0),silent=TRUE),"try-error")) {
        diag(h0) <- diag(h0) + m
        m <- 10*m
        cat(m)
      }
    }
  }
}

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}
newt(c(0,0),rb,gb,hb,k=2)