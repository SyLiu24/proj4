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