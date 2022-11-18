#XXXXXXXXXXXXXXXXXXXXXXXX proj2.r XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Group member
# Biwei Zhu, s2325784
# Guanhao Su, s2301705
# Shuying Liu, s2436365

# Github repo: 
# https://github.com/SyLiu24/proj4

#XXXXXXXXXXXXXXXXXXXXXXXXXXX Contribution XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Biwei Zhu completed the approximation to Hessian
# Guanhao Su completed the perturbation of Hessian
# Biwei Zhu and Guanhao Su performed tests to check the algorithm
# Shuying Liu complted the rest of the newt function

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# The aim of this project is to write an R function, newt, implementing 
# Newton’s method for minimization of functions.

# Suppose the objective function D(theta) to minimize is smooth with 2 bounded 
# derivatives - gradient g and Hessian H. Newton's method attempts to start
# from an initial guess theta, then produce a sequence of iterates
# minimizing successive quadratic approximations to D, to converge towards 
# a minimizer theta^hat. 
# In detail, the estimation quadratic to produce the iterates is the 
# second-order Taylor approximation 
# D(theta+delta) = D(theta) + delta^Tg(theta) + 1/2*delta^TH(theta)delta.
# Provided H(theta) is positive definite, the quadratic has minimizer
# delta = -(H(theta))^{-1}g(theta), note this also true for any positive
# definite matrix in place of H. Hence the next iterate theta+delta.
# To ensure convergence, first if H(theta) is not positive definite then we just
# perturb it to be by adding a multiple of the identity matrix large enough to
# force positive definiteness (tested by Cholesky decomposition); Secondly to
# avoid overshoot, we halve the step until it reduces the objective.
# We would repeat minimizing the quadratic to find improved guesses until g is
# close enough to zero.

newt <- function(theta,func,grad,hess=NULL,...,
                 tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6) {
  # Use Newton's method to find minimum of func
  
  # Input:
  #   theta - vector of initial values for the optimization parameters
  #   func - the objective function to minimize; Its first argument is the 
  #   vector of optimization parameters. Remaining arguments will be passed 
  #   from newt using ‘...
  #   grad - the gradient function. It has the same arguments as func.
  #   hess - the Hessian matrix function. It has the same arguments as func.
  #   If not supplied then newt should obtain an approximation to the
  #   Hessian by finite differencing of gradient.
  #   ... - any arguments of func, grad and hess after the first are passed
  #   tol - the convergence tolerance
  #   fscale - a rough estimate of the magnitude of func near the optimum
  #   maxit - the maximum number of Newton iterations to try
  #   max.half the maximum number of times a step can be halved
  #   eps - the finite difference intervals
  
  # Returns a list containing
  #   f - the minimum of the objective
  #   theta - the value of the parameters at the minimum
  #   iter - the number of iterations taken to reach the minimum
  #   g - the gradient vector at the minimum
  #   Hi - the inverse of the Hessian matrix at the minimum
  
  n <- length(theta)
  
  # If hess is not given, then estimate it by finite differencing of grad
  if (is.null(hess)) {
    hess <- function(theta,...) {
      # Produces hess estimation function that has the same arguments as func
      g <- grad(theta,...)
      hfd <- matrix(0,n,n) # estimated Hessian
      
      # Generate finite defferencing for each element of theta
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

  # Generate initial guess and its derivatives
  th0 <- theta
  f0 <- func(th0,...)
  g0 <- grad(th0,...)
  h0 <- hess(th0,...)

  # Check if the initial theta is suitable for optimization
  # Check if objective function value is -inf
  if (is.infinite(f0) & f0 < 0)
    sprintf("Function is unbounded below at %f.",th0) |> stop()
  # Check if objective or derivatives are not finite
  else if (is.infinite(f0) | any(is.infinite(g0)) | any(is.infinite(h0)))
    stop("Objective or derivatives are not finite at the initial theta.")
  
  it <- 0 # number of iterations
  while (any(abs(g0) >= tol*abs(f0+fscale))) { # convergence condition
    # Check if max iteration time has been reached
    if (it == maxit)
      stop("Maxit is reached without convergence.")
    
    # If h0 is not positive definite then perturb it to be
    # by adding a multiple of identity matrix
    if (inherits(try(chol(h0),silent=TRUE),"try-error")) {
      h0 <- as.matrix(h0) # avoid norm and diag error if h0 is 1x1
      m <- (1e-6)*norm(h0) # Multiplier of identity matrix
      while (inherits(try(chol(h0),silent=TRUE),"try-error")) {
        diag(h0) <- diag(h0) + m
        m <- 10*m
      }
    }
    
    # Calculate the step
    R <- chol(h0)
    step <- -backsolve(R,forwardsolve(t(R),g0))
    
    # Step halving to make sure step reduce the objective
    for (it_half in 0:max.half) {
      th1 <- th0 + step # new iterate
      f1 <- func(th1,...)
      
      # Check if objective value is valid
      if (is.na(f1) | is.nan(f1))
        sprintf("No objective value at %f.",th1) |> stop()
      else if (is.infinite(f1) & f1 < 0)
        sprintf("Function is unbounded below at %f.",th1) |> stop()
      # Check if step halving is needed
      else if (f1 >= f0) step <- step/2
      else break
    }
    
    # Check if f1 reduces after max half step halvings
    if (f1 >= f0)
      stop("Step fails to reduce the objective after max.half step halvings")
    
    # Update iterates
    th0 <- th1
    # cat(th0,'\n')
    f0 <- f1
    g0 <- grad(th1,...)
    h0 <- hess(th1,...)
    it <- it + 1
  }
  # Check if Hessian is positive definite at convergence
  if (inherits(try(chol(h0),silent=TRUE),"try-error")) {
    warning("Hessian is not positive definite at convergence, 
                its inverse is not returned.")
    return(list(f=f0,theta=th0,iter=it,g=g0))
  }
  else {
    return(list(f=f0,theta=th0,iter=it,g=g0,
                Hi=chol2inv(chol(h0))))
  }
}
