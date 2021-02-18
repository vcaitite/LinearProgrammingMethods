#  ---------------------------------------------------------------------------
#  A Fundamental Interior Point Algorithm for Solving LP Problem
#  Vannelli.(1993). "Teaching Large-Scale Optimization
#  by an Interior Point Approach"; IEEE Trans. on Education. (36)1:204-209
#  Solve the linear programming problem
#  min  c'x
#  s.t. Ax < b
#  and returns the path x to the solution with respective 
#  objective function values fx.
#  ---------------------------------------------------------------------------

# A - An m2 by n matrix of coefficients of constraints.
# b -	A vector of length m2 giving the right hand side of constraints. 
#	c - A vector of length n which gives the coefficients of the objective function.
# gama - 
# x0 - initial point 
# maxIter - maximum number of iterations
# epslon - 
InteriorPoint <- function(A, b, c, gama, x0, maxIter, epslon, plot){
  function_value <- rep(0, maxIter);
  x <- matrix(ncol = nrow(c), nrow = maxIter)
  function_value[1] <- t(c)%*%x0
  x[1, ] <- t(x0) 
  iterations <- 0
  stop_condition <- FALSE
  while(stop_condition == FALSE && iterations <= maxIter){
    # Step 2:
    vk <- b - A%*%x0
    print(vk)
    
    # Step 3: 
    Dk = diag(rep(1,length(vk)))
    for (k in 1:length(vk)) {
      Dk[k,k]<-1/vk[k]
    }
    print(format(Dk, digits = 4))
    
    # Step 4:
    dx<-pseudoinverse(t(A) %*% Dk %*% Dk %*% A) %*% c
    print(format(dx, digits = 4))
    
    # Step 5:
    dv <- (-A) %*% dx
    print(format(dv, digits = 4))
    
    # Step 6: 
    step = matrix(rep(0,length(vk)), ncol = 1)
    for (k in 1:length(step)) {
      step[k]<-vk[k]/dv[k]
      if (step[k] > 0){
        step[k] <- -1000000
      }
    }
    print(format(step, digits = 4))
    alfa <- max(step)*gama
    print(alfa)
    
    # Step 7:
    x1 <- x0 - alfa * dx
    print(x1)
    iterations <- iterations + 1
    print(function_value[iterations])
    function_value[iterations+1] <- t(c) %*% x1
    
    
    # Condição de parada
    yk <- -(Dk %*% Dk) %*% dv
    f <- abs(t(b)%*%yk - t(c)%*%x1)/max(c(1, abs(t(c)%*%x1))) 
    print(f)
    if (f < 10^(-6)){
      stop_condition <- TRUE;
    }
    x[iterations+1, ] <- t(x1[,1]) 
    x0 <- x1
  }
  result <- rep(0, iterations+1)
  xresult <- matrix(ncol = ncol(x), nrow = iterations+1)
  print(dim(xresult))
  print(dim(x))
  print(iterations)
  for (index in 1:(iterations+1)) {
    result[index] <- function_value[index]
    xresult[index, ] <- x[index, ]
  }
  if(plot == TRUE){
    plot(1:(iterations+1), result, type="b")
  }
  return (list(result, xresult, iterations))
}